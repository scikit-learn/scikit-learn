"""Checks that dist/* contains the number of wheels built from the
.github/workflows/wheels.yml config."""
import sys
from pathlib import Path

import yaml

gh_wheel_path = Path.cwd() / ".github" / "workflows" / "wheels.yml"
with gh_wheel_path.open("r") as f:
    wheel_config = yaml.safe_load(f)

build_matrix = wheel_config["jobs"]["build_wheels"]["strategy"]["matrix"]["include"]
n_wheels = len(build_matrix)

# plus one more for the sdist
n_wheels += 1

# arm64 builds from cirrus
cirrus_path = Path.cwd() / "build_tools" / "cirrus" / "arm_wheel.yml"
with cirrus_path.open("r") as f:
    cirrus_config = yaml.safe_load(f)

n_wheels += len(cirrus_config["macos_arm64_wheel_task"]["matrix"])
n_wheels += len(cirrus_config["linux_arm64_wheel_task"]["matrix"])

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
