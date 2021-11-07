"""Checks that dist/* contains the number of wheels built from the
.github/workflows/wheels.yml config."""
import yaml
from pathlib import Path
import sys

gh_wheel_path = Path.cwd() / ".github" / "workflows" / "wheels.yml"
with gh_wheel_path.open("r") as f:
    wheel_config = yaml.safe_load(f)

build_matrix = wheel_config["jobs"]["build_wheels"]["strategy"]["matrix"]
n_python_versions = len(build_matrix["python"])

# Python 3.7, 3.8, and 3.9 each produces 7 wheels
# 1 osx wheel (x86_64)
# 4 linux wheel (i686 + x86_64) * (manylinux1 + manylinux2010)
# 2 windows wheel (win32 + wind_amd64)
n_wheels = 7 * n_python_versions

# NumPy on Python 3.10 only supports 64bit and is only avaliable with manylinux2014
# With macos and window support the number of wheels should go up to 3
n_wheels += 3

# plus one more for the sdist
n_wheels += 1

# aarch64 builds from travis
travis_config_path = Path.cwd() / ".travis.yml"
with travis_config_path.open("r") as f:
    travis_config = yaml.safe_load(f)

jobs = travis_config["jobs"]["include"]
travis_builds = [j for j in jobs if any("CIBW_BUILD" in env for env in j["env"])]
n_wheels += len(travis_builds)

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
