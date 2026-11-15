# QPE cosmology: verification package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24247908.svg)](https://doi.org/10.5281/zenodo.24247908)

This repository contains the code that produces every number and table in

> Chen, D.-M. 2026, *Cosmological Tests of Scale-Dependent Quantum Gravity:
> Type Ia Supernovae, Full-Spectrum CMB, and Baryon Acoustic Oscillations*,
> MNRAS (submitted).

The aim is narrow and practical: **anyone can check the published fit quality
in a few seconds, and reproduce it from scratch in under twenty minutes.**

---

## 1. Quick start (about 1 second)

```bash
pip install -r requirements.txt
python get_data.py   # downloads the Pantheon+ catalogue (~0.6 MB)
python verify.py --fast
```

Actual output:

```
==============================================================
QPE cosmology: verification of Chen (2026)
==============================================================

[1/3] Type Ia supernovae (Pantheon+SH0ES, N = 1701)
  computed in 0.1 s
  chi2                     745.5586       paper: 745.6     OK
  chi2 / N                   0.4383       paper: 0.438     OK

[2/3] BAO (DESI DR2, 13 observables)
  computed in 0.0 s
  chi2                       8.0119       paper: 8         OK

[3/3] CMB (Planck 2018 TT, N = 836)
  SKIPPED (--fast).

==============================================================
All published values reproduced.
==============================================================
```

The SN step needs `data/Pantheon+SH0ES.dat` (see section 3); without it that
block is skipped with an instruction message, and the rest still runs.

## 2. Full verification, including the CMB (about 17 minutes)

```bash
python verify.py
```

Actual output of the third block:

```
[3/3] CMB (Planck 2018 TT, N = 836)
  building the spherical-Bessel grid (NU=900, NK=300) ...
  this is the expensive step: ~1 GB of RAM and 5-15 min.
  computed in 1029.0 s
  chi2                    1464.5716       paper: 1464.7    OK
  chi2 / N                   1.7519       paper: 1.752     OK
```

The CMB step builds a spherical-Bessel grid of size
`836 multipoles x 300 wavenumbers x 900 radii` (about 1 GB in float32) and then
evaluates the two-branch projection of Eqs. (149)-(150). This is the only
expensive part, and it is recomputed from scratch.

A cheaper sanity check, useful to confirm that the code runs:

```bash
python verify.py --nu 300 --nk 150
```

This takes about a minute and does **not** reproduce the published value; the
publication grid (NU=900, NK=300) is converged to about 2 per cent.

## 3. Data

| File | Size | Source |
|---|---|---|
| `data/Pantheon+SH0ES.dat` | 0.6 MB | Pantheon+ release (Scolnic et al. 2022); see below |
| `data/COM_PowerSpect_CMB-TT-full_R3.01.txt` | 0.2 MB | included (Planck Legacy Archive, public) |
| DESI DR2 BAO | 13 numbers | included in `qpe_model.py` as `DESI_DR2` |

The Pantheon+ catalogue is **not** redistributed, because its terms of use ask
users to obtain it from the collaboration. Download `Pantheon+SH0ES.dat` from
the Pantheon+ data release and place it in `data/`. The columns used are
`zHD`, `MU_SH0ES`, `MU_SH0ES_ERR_DIAG`; the script selects the 1701 supernovae
with `z > 0` and a positive uncertainty.

## 4. What is in the code

`qpe_model.py` is a single self-contained module, about 400 lines. The comments
give the equation numbers of the paper, so that the code can be read against the
text.

| Function | Paper | Notes |
|---|---|---|
| `gamma_of_z` | Eq. (86) | three constant segments, `tanh` joins |
| `x_of_z`, `z_of_x` | Eq. (75) | `1+z = 1/sqrt(1-2 sinh x)` |
| `d_prop_of_z` | Eq. (78) | light-travel distance |
| `C_of_z` | Eq. (84) | areal radius `C = (1+z)^(gamma-2) D_prop` |
| `D_L_of_z`, `mu_of_z` | Eqs. (82), (85) | luminosity distance, distance modulus |
| `K_gamma` | Eq. (153) | cumulative path integral for radial BAO |
| `sn_chi2` | Eq. (90) | SN Ia Hubble diagram |
| `bao_model`, `bao_chi2` | Eqs. (152)-(154) | `D_M/r_d`, `D_H/r_d`, `D_V/r_d` |
| `emission_kernel` | Eq. (142) | `W ~ u^(1-p/2) exp(-tau)` |
| `power_spectrum` | Eq. (124) | three-segment `n(k)`, `cos^2(k/k_w)` modulation |
| `CMBProjector` | Eqs. (149)-(150) | two-branch projection, `j_l(kC)` |

The published parameter values are collected at the top of the module in
`SN_PARAMS`, `BAO_PARAMS` and `CMB_PARAMS`.

Two points where the code is deliberately explicit:

* `CMB_PARAMS` holds the **full-precision** posterior medians of the archived
  chain, not the rounded values printed in Table 2 of the paper. Using the
  rounded values shifts the CMB chi^2 by about 6.
* `CMBProjector` fixes the behaviour of the two-branch change of variable for
  `u` below the value of `C` at the horizon (the `left=`/`right=` arguments
  of `np.interp`). This affects 9 of the 900 radial grid points and changes
  chi^2 by 0.4 per cent, but it is fixed here so that the published value is
  reproduced exactly.

## 5. Runtime

| Step | Time | Memory |
|---|---|---|
| `verify.py --fast` (SN + BAO) | ~1 s | < 100 MB |
| `verify.py` (adds CMB) | ~17 min | ~1.5 GB |
| `verify.py --nu 300 --nk 150` | ~1 min | ~200 MB |

Timings on a 2026 desktop CPU, single core.

## 6. Not included in this release

* The MCMC samplers (`emcee`) that produced the chains; the chains themselves
  are available from the author until the archival deposit is posted.
* The scripts that draw Figs. 2, 4, 5 and B1. Each is a thin wrapper around
  `qpe_model.py`: for Fig. 4, `CMBProjector(...).Dl_model(CMB_PARAMS)`
  reproduces the model curve.

## 7. Licence and citation

MIT licence. If you use this code, please cite the paper above. A versioned
archival copy with a DOI will be deposited on Zenodo; the DOI will be added
here and in the paper's data-availability statement.

## 8. Known limitations

* The BAO chi^2 uses the **diagonal** uncertainties reported with the DESI DR2
  measurements. The `D_M`-`D_H` correlations within a tracer are not included;
  the paper states this, and the same treatment is applied to the Lambda-CDM
  comparison.
* Only the best-fit points are verified, not the full posteriors.
* The CMB computation follows the exact grid and settings of the published
  analysis; changing `--nu` or `--nk` shifts the result.
