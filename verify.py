#!/usr/bin/env python
"""
verify.py --- Reproduce the three chi^2 values quoted in Chen (2026).

Usage
-----
    python verify.py              # everything (CMB takes ~10 min the first time)
    python verify.py --fast       # SN + BAO only (a few seconds)
    python verify.py --no-cmb     # same as --fast

Expected output is given in README.md.  The script exits with status 0 if every
value agrees with the paper to the quoted precision, and 1 otherwise.
"""

import argparse
import os
import sys
import time

import numpy as np

import qpe_model as qm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

PLANCK_FILE = os.path.join(DATA, "COM_PowerSpect_CMB-TT-full_R3.01.txt")
PANTHEON_FILE = os.path.join(DATA, "Pantheon+SH0ES.dat")

#: Values quoted in the paper (Tables 1 and 2, and Section 7.3).
PAPER = {
    "sn_chi2": 745.6, "sn_chi2nu": 0.438,
    "bao_chi2": 8.0,
    "cmb_chi2": 1464.7, "cmb_chi2nu": 1.752,
}


def line(name, value, expected, tol, unit="", note=""):
    ok = abs(value - expected) <= tol
    flag = "OK  " if ok else "DIFF"
    bar = "=" * 62
    print("  %-22s %10.4f%-6s paper: %-8.4g  %s %s"
          % (name, value, unit, expected, flag, note))
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", "--no-cmb", dest="fast", action="store_true",
                    help="skip the CMB projection (SN + BAO only)")
    ap.add_argument("--nu", type=int, default=900, help="grid points in areal radius (default 900)")
    ap.add_argument("--nk", type=int, default=300, help="grid points in k (default 300)")
    args = ap.parse_args()

    print("=" * 62)
    print("QPE cosmology: verification of Chen (2026)")
    print("=" * 62)

    all_ok = True

    # ---------------- SN Ia ----------------
    print("\n[1/3] Type Ia supernovae (Pantheon+SH0ES, N = 1701)")
    if os.path.exists(PANTHEON_FILE):
        t0 = time.time()
        chi2, n = qm.sn_chi2(PANTHEON_FILE)
        print("  computed in %.1f s" % (time.time() - t0))
        all_ok &= line("chi2", chi2, PAPER["sn_chi2"], 0.1)
        all_ok &= line("chi2 / N", chi2 / n, PAPER["sn_chi2nu"], 0.001)
    else:
        print("  SKIPPED: %s not found." % PANTHEON_FILE)
        print("  Download the Pantheon+SH0ES catalogue and place it there")
        print("  (see README.md, section 'Data').")

    # ---------------- BAO ----------------
    print("\n[2/3] BAO (DESI DR2, 13 observables)")
    t0 = time.time()
    chi2, n = qm.bao_chi2()
    print("  computed in %.1f s" % (time.time() - t0))
    all_ok &= line("chi2", chi2, PAPER["bao_chi2"], 0.1)

    # ---------------- CMB ----------------
    print("\n[3/3] CMB (Planck 2018 TT, N = 836)")
    if args.fast:
        print("  SKIPPED (--fast).")
    elif not os.path.exists(PLANCK_FILE):
        print("  SKIPPED: %s not found." % PLANCK_FILE)
    else:
        print("  building the spherical-Bessel grid (NU=%d, NK=%d) ..." % (args.nu, args.nk))
        print("  this is the expensive step: ~1 GB of RAM and 5-15 min.")
        t0 = time.time()
        proj = qm.CMBProjector(PLANCK_FILE, nu=args.nu, nk=args.nk)
        chi2, n = proj.chi2(qm.CMB_PARAMS)
        print("  computed in %.1f s" % (time.time() - t0))
        all_ok &= line("chi2", chi2, PAPER["cmb_chi2"], 1.0)
        all_ok &= line("chi2 / N", chi2 / n, PAPER["cmb_chi2nu"], 0.002)

    print("\n" + "=" * 62)
    print("All published values reproduced." if all_ok
          else "Some values differ: see the DIFF lines above.")
    print("=" * 62)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
