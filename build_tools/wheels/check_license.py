"""Checks the bundled license is installed with the wheel."""

import platform
import site
from itertools import chain
from pathlib import Path

site_packages = site.getsitepackages()

site_packages_path = (Path(p) for p in site_packages)

try:
    distinfo_path = next(
        chain(
            s
            for site_package in site_packages_path
            for s in site_package.glob("scikit_learn-*.dist-info")
        )
    )
except StopIteration as e:
    raise RuntimeError("Unable to find scikit-learn's dist-info") from e

license_text = (distinfo_path / "COPYING").read_text()

assert "Copyright (c)" in license_text

assert (
    "This binary distribution of scikit-learn also bundles the following software"
    in license_text
), f"Unable to find bundled license for {platform.system()}"
