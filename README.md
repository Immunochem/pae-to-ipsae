# pae-to-ipsae

**Make ESMFold2 PAE files readable by ipSAE — and refuse them when they don't match.**

[ipSAE](https://github.com/DunbrackLab/IPSAE2) scores an interface from the predicted
aligned error, keeping only the residue pairs that are actually well-aligned. It reads a
PAE as `.json`, `.json.gz` or `.npz`. It **rejects a bare `.npy`**:

```
native_iptm.py:155   raise ValueError(f"PAE input must be JSON or NPZ: {only}")
readers.py:57        raise ValueError(f"PAE file must be .json or .npz: {path}")
```

ESMFold2 writes exactly that: a bare `.npy`, with no chain metadata. So ipSAE cannot be
applied to ESMFold2 output without a conversion step. This is that step — with the guard
that makes it safe.

## Why a converter needs a guard

A PAE matrix with no token map is a blob. A blob of the **wrong size** does not crash
anything: it produces a score that is silently wrong. So the conversion is refused unless

```
n_tokens(PAE) == n_residues(structure)
```

The tool never resizes, pads or truncates to make a mismatch fit. It stops and says so.

```bash
python3 pae_to_ipsae.py model.pae.npy model.cif model.npz
python3 pae_to_ipsae.py --autotest     # 4 cases, 2 of which MUST be refused
```

Then feed the result to ipSAE as usual:

```bash
python3 -m ipsae.cli model.npz model.cif --output-file out.tsv
```

## ⚠️ The orientation caveat — read this before using directional scores

The row/column orientation of an ESMFold2 `.pae.npy` is **not specified anywhere**: not in
the `.npy` header (which carries only dtype and shape), not in the ESMFold2 repository.
We looked. Until a producer-side source settles it, this matters differently per score:

| score | effect of a flipped orientation | safe to use? |
|---|---|---|
| `ipsae` (max over both directions) | none | ✅ yes |
| `iptm` (symmetric by construction) | none | ✅ yes |
| `ipsae12`, `ipsae21`, `iPAE12`, `iPAE21` | **swapped** | ❌ no |

The tool does not guess. It writes the caveat into the `.npz` (key `_reserve`) and prints
it on every run. To lift it, read the **producer**, not the values.

## What it deliberately does not do

**It does not convert formats ipSAE already reads.** AlphaFold3-style predictors — and
others writing `token_pair_pae` in JSON — are accepted as-is (`readers.py:29,45`).
Converting something that needs no conversion only adds a place for errors.

## Autotest

```
[1] residue count     : 8 residues {'A': 5, 'B': 3}  (8 expected)
[2] PAE 8x8 / 8 res   : ACCEPTED  n_tokens=8
    caveat written    : yes
[3] PAE 7x7 / 8 res   : REFUSED
[4] .json as input    : REFUSED
✅ AUTOTEST GREEN
```

Case [3] is the one that matters: **if it ever passes, the guard is broken.**

## Requirements

`numpy`. Python 3.9+. No other dependency — ipSAE itself is called separately.

## Credit

ipSAE is by **Roland L. Dunbrack, Jr.** ([IPSAE2](https://github.com/DunbrackLab/IPSAE2),
MIT). This repository contains **no code from ipSAE**; it only prepares input for it.

Released by [Immunochem](https://immunochem.be) under the MIT license. It came out of a
nanobody design campaign where ipSAE had to be applied to two co-folding engines at once,
one of which wrote `.npy`.
