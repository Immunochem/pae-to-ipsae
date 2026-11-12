#!/usr/bin/env python3
"""pae_vers_ipsae.py — donner a ipSAE une PAE qu il sait lire, ou la refuser.

    python3 tools/ipsae/pae_vers_ipsae.py <pae_in> <structure> <npz_out>
    python3 tools/ipsae/pae_vers_ipsae.py --autotest

IPSAE2 (DunbrackLab, MIT) accepte une PAE en `.json`, `.json.gz` ou `.npz`. Il
REFUSE le `.npy` nu :

    native_iptm.py:155   raise ValueError(f"PAE input must be JSON or NPZ: {only}")
    readers.py:57        raise ValueError(f"PAE file must be .json or .npz: {path}")

Or ESMFold2 ecrit precisement un `.npy` nu, sans metadonnee de chaine. Ce module
fait la conversion AVEC SON GARDE, parce qu une PAE sans carte de tokens est un blob
et qu un blob de la mauvaise taille produit un score silencieusement faux.

────────────────────────────────────────────────────────────────────────────────
🔴 RESERVE D ORIENTATION, NON LEVEE, A LIRE AVANT D UTILISER ipsae12/ipsae21.

L orientation lignes/colonnes du `.pae.npy` d ESMFold2 n est fixee NULLE PART :
ni l en-tete `.npy` (qui ne porte que dtype et shape), ni le depot ESMFold2, ni
nos protocoles. Mesure du 2026-11-10, consignee dans
measured 2026-11-10 by reading the producer and the `.npy` header; see README.

Consequence exacte, et elle n est pas la meme pour tous les scores :
  - `ipsae`  = max des deux directions          -> INSENSIBLE a l orientation. UTILISABLE.
  - `iptm`   = symetrique par construction      -> INSENSIBLE. UTILISABLE.
  - `ipsae12` / `ipsae21` / `iPAE12` / `iPAE21` -> ECHANGES si l orientation est
    inversee. NON FIABLES tant que la reserve tient.

Ce module ne devine pas : il ECRIT la reserve dans le `.npz` (cle `_reserve`) et
la rappelle a chaque appel. Pour la lever il faut lire le PRODUCTEUR (le skill
`cofold-esmfold2`), pas raisonner sur les valeurs.
────────────────────────────────────────────────────────────────────────────────

🔸 Les formats deja lisibles par ipSAE ne passent PAS par ici : OpenDDE et Mol-1
   ecrivent un `.json` avec `token_pair_pae`, accepte tel quel
   (`readers.py:29,45`). Convertir ce qui n a pas besoin de l etre, c est ajouter
   une source d erreur sans rien gagner.
"""
import json
import pathlib
import sys

RESERVE = (
    "orientation lignes/colonnes du .pae.npy NON ETABLIE "
    "see README : ipsae/iptm safe, "
    "ipsae12/ipsae21/iPAE12/iPAE21 NON FIABLES"
)


def compter_residus(chemin):
    """Nombre de residus de la structure, par CA pour un PDB, par unite pour un mmCIF.

    🔴 On compte les RESIDUS, pas les atomes : la PAE d ESMFold2 est par residu.
    Rend (n, detail_par_chaine) ; leve ValueError si le format est inconnu.
    """
    p = pathlib.Path(chemin)
    suf = "".join(p.suffixes).lower()
    par_chaine = {}
    if suf.endswith(".pdb") or suf.endswith(".pdb.gz"):
        ouvrir = __import__("gzip").open if suf.endswith(".gz") else open
        with ouvrir(p, "rt") as f:
            for l in f:
                if l.startswith("ATOM") and l[12:16].strip() == "CA":
                    par_chaine.setdefault(l[21:22], []).append(l[22:27])
    elif suf.endswith(".cif") or suf.endswith(".mmcif") or suf.endswith(".cif.gz"):
        ouvrir = __import__("gzip").open if suf.endswith(".gz") else open
        cols, dans = [], False
        with ouvrir(p, "rt") as f:
            for l in f:
                s = l.strip()
                if s.startswith("_atom_site."):
                    cols.append(s.split(".", 1)[1])
                    dans = True
                    continue
                if dans:
                    if s.startswith("ATOM") or s.startswith("HETATM"):
                        c = s.split()
                        i = {k: j for j, k in enumerate(cols)}
                        ka = i.get("label_asym_id", i.get("auth_asym_id"))
                        kn = i.get("label_seq_id", i.get("auth_seq_id"))
                        katm = i.get("label_atom_id", i.get("auth_atom_id"))
                        if ka is None or kn is None:
                            raise ValueError(f"colonnes de chaine/residu absentes dans {p}")
                        if katm is not None and c[katm].strip('"') != "CA":
                            continue
                        par_chaine.setdefault(c[ka], []).append(c[kn])
                    elif par_chaine:
                        break
    else:
        raise ValueError(f"extension de structure inconnue : {p.name} "
                         "(.pdb, .pdb.gz, .cif, .mmcif, .cif.gz attendus)")
    detail = {k: len(set(v)) for k, v in par_chaine.items()}
    return sum(detail.values()), detail


def convertir(pae_in, structure, npz_out, cle="pae"):
    """Convertit une PAE vers un .npz lisible par ipSAE, ou leve.

    🔴 GARDE BLOQUANT : n_tokens(PAE) == n_residus(structure). Une PAE dont la
       taille ne correspond pas a sa structure produit un score silencieusement
       faux — c est le cas que ce module existe pour attraper.
    Rend un dict de compte-rendu : il est IMPRIME par l appelant, a chaque appel.
    """
    import numpy as np
    p = pathlib.Path(pae_in)
    suf = "".join(p.suffixes).lower()
    if suf.endswith(".npy"):
        a = np.load(p)
    elif suf.endswith(".npz"):
        z = np.load(p)
        noms = list(z.files)
        k = next((n for n in (cle, "pae", "token_pair_pae", "predicted_aligned_error")
                  if n in noms), None)
        if k is None:
            raise ValueError(f"{p.name} : aucune cle de PAE parmi {noms}")
        a = z[k]
    elif suf.endswith(".json") or suf.endswith(".json.gz"):
        raise ValueError(f"{p.name} : deja lisible par ipSAE, ne pas convertir "
                         "(readers.py:29,45 acceptent token_pair_pae tel quel)")
    else:
        raise ValueError(f"extension de PAE inconnue : {p.name}")

    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"PAE non carree : shape {a.shape}")
    n_pae = int(a.shape[0])
    n_res, detail = compter_residus(structure)
    if n_pae != n_res:
        raise ValueError(
            f"🔴 GARDE : PAE {n_pae}x{n_pae} contre {n_res} residus dans "
            f"{pathlib.Path(structure).name} ({detail}). Une PAE sans carte de tokens "
            "est un blob : on refuse, on ne redimensionne pas.")

    out = pathlib.Path(npz_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, **{cle: a.astype("float32")},
             _reserve=np.array([RESERVE], dtype=object))
    return {"pae_in": str(p), "structure": str(structure), "npz_out": str(out),
            "n_tokens": n_pae, "n_residus": n_res, "chaines": detail,
            "reserve": RESERVE}


def autotest():
    """Trois cas, dont un qui DOIT echouer. Si le 3e passe, le garde ne marche pas."""
    import tempfile
    import numpy as np
    ok = True
    with tempfile.TemporaryDirectory() as d:
        d = pathlib.Path(d)
        # une structure jouet : 2 chaines, 5 + 3 residus = 8
        pdb = d / "jouet.pdb"
        lignes = []
        for i in range(1, 6):
            lignes.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                          f"{i*3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           C")
        for i in range(1, 4):
            lignes.append(f"ATOM  {i+10:5d}  CA  GLY B{i:4d}    "
                          f"{i*3.8:8.3f}{10.0:8.3f}{0.0:8.3f}  1.00 50.00           C")
        pdb.write_text("\n".join(lignes) + "\nEND\n")

        n, det = compter_residus(pdb)
        print(f"[1] comptage          : {n} residus {det}  (8 attendus)")
        ok &= (n == 8)

        bon = d / "bon.npy"
        np.save(bon, np.random.rand(8, 8).astype("float32") * 30)
        try:
            cr = convertir(bon, pdb, d / "bon.npz")
            print(f"[2] PAE 8x8 / 8 res   : ACCEPTE  n_tokens={cr['n_tokens']}")
            z = np.load(d / "bon.npz", allow_pickle=True)
            r = str(z["_reserve"][0])
            print(f"    reserve ecrite    : {'oui' if 'NON ETABLIE' in r else '🔴 NON'}")
            ok &= ("NON ETABLIE" in r)
        except Exception as e:
            print(f"[2] 🔴 refus inattendu : {e}")
            ok = False

        faux = d / "faux.npy"
        np.save(faux, np.random.rand(7, 7).astype("float32") * 30)
        try:
            convertir(faux, pdb, d / "faux.npz")
            print("[3] 🔴 PAE 7x7 / 8 res : ACCEPTE — LE GARDE NE MARCHE PAS")
            ok = False
        except ValueError as e:
            print(f"[3] PAE 7x7 / 8 res   : REFUSE  « {str(e).splitlines()[0][:72]}… »")

        try:
            convertir(d / "rien.json", pdb, d / "x.npz")
            print("[4] 🔴 .json : ACCEPTE — on convertit ce qui n en a pas besoin")
            ok = False
        except ValueError as e:
            print(f"[4] .json en entree   : REFUSE  « {str(e)[:72]}… »")

    print("\n" + ("✅ AUTOTEST VERT" if ok else "🔴 AUTOTEST ROUGE"))
    return 0 if ok else 1


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--autotest":
        return autotest()
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    try:
        cr = convertir(sys.argv[1], sys.argv[2], sys.argv[3])
    except Exception as e:
        print(f"🔴 REFUS : {e}")
        return 3
    print(json.dumps({k: v for k, v in cr.items() if k != "reserve"},
                     indent=1, ensure_ascii=False))
    print(f"🔸 {cr['reserve']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
