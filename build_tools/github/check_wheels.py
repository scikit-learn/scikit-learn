"""Checks that dist/* contains the number of wheels."""
from pathlib import Path
import sys

# Python 3.7, 3.8, 3.9 each has 7 wheels
# 1 osx wheel (x86_64)
# 4 linux wheel (i686 + x86_64) * (manylinux1 + manylinux2010)
# 2 windows wheel (win32 + wind_amd64)
n_wheels = 3 * 7

# NumPy on Python 3.10 only supports 64bit and is only avaliable with manylinux2014
# With macos and window support the number of wheels should go up to 3
n_wheels += 1

# plus one more for the sdist
n_wheels += 1

# aarch64 builds from travis has builds for Python 3.7, 3.8, 3.9
n_wheels += 3

dist_files = list(Path("dist").glob("**/*"))
n_dist_files = len(dist_files)

if n_dist_files != n_wheels:
    print(
        f"Expected {n_wheels} wheels in dist/* but "
        f"got {n_dist_files} artifacts instead."
    )
    sys.exit(1)

print(f"dist/* has the expected {n_wheels} wheels:")
print("\n".join(file.name for file in dist_files))
