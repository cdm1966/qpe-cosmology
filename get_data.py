#!/usr/bin/env python
"""
get_data.py --- Download the Pantheon+SH0ES catalogue used by verify.py.

The catalogue is not redistributed with this package: the Pantheon+ Data
Release carries no explicit licence, so the file is fetched from the official
repository each time instead of being copied.

Usage
-----
    python get_data.py

The file is written to data/Pantheon+SH0ES.dat (579 283 bytes).
"""

import os
import sys
import urllib.request

URL = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
       "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat")
TARGET = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "data", "Pantheon+SH0ES.dat")
EXPECTED = 579283          # bytes, as distributed


def main():
    os.makedirs(os.path.dirname(TARGET), exist_ok=True)
    if os.path.exists(TARGET) and os.path.getsize(TARGET) == EXPECTED:
        print("Already present: %s (%d bytes)" % (TARGET, EXPECTED))
        return 0

    print("Downloading the Pantheon+SH0ES catalogue from the official release")
    print("  %s" % URL)
    try:
        urllib.request.urlretrieve(URL, TARGET)
    except Exception as exc:                       # noqa: BLE001
        print("Download failed: %s" % exc)
        print("Download the file manually from")
        print("  https://github.com/PantheonPlusSH0ES/DataRelease")
        print("and save it as %s" % TARGET)
        return 1

    size = os.path.getsize(TARGET)
    if size != EXPECTED:
        print("Warning: got %d bytes, expected %d." % (size, EXPECTED))
        return 1
    print("Saved %s (%d bytes)" % (TARGET, size))
    return 0


if __name__ == "__main__":
    sys.exit(main())
