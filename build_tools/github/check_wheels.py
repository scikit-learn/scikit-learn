"""Checks that dist/* contains the number of wheels built from the
.github/workflows/wheels.yml config."""
import yaml
from pathlib import Path
import sys

gh_wheel_path = Path.cwd() / ".github" / "workflows" / "wheels.yml"
with gh_wheel_path.open("r") as f:
    wheel_config = yaml.safe_load(f)

build_matrix = wheel_config["jobs"]["build_wheels"]["strategy"]["matrix"]["include"]
n_wheels = len(build_matrix)

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
