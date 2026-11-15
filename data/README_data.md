# Data files

## `COM_PowerSpect_CMB-TT-full_R3.01.txt`
Planck 2018 TT power spectrum, public release (Planck Legacy Archive).
Redistributed here without restriction.

## `Pantheon+SH0ES.dat`
**NOT redistributed.** The Pantheon+ collaboration asks users to obtain the
catalogue from the official release. The copy in this folder (if present) is a
local test copy and must be deleted before this repository is made public; it is
listed in `.gitignore`.

Columns used by `verify.py`: `zHD`, `MU_SH0ES`, `MU_SH0ES_ERR_DIAG`.

## DESI DR2 BAO
The thirteen measurements are reproduced verbatim in `qpe_model.py` as the
array `DESI_DR2` (Section 7.2 of the paper), so no external file is needed.
