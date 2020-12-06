"""Checks that dist/* contains the number of wheels built from the
.github/workflows/wheels.yml config."""
import yaml
from pathlib import Path
import sys

gh_wheel_path = Path.cwd() / ".github" / "workflows" / "wheels.yml"
with gh_wheel_path.open('r') as f:
    wheel_config = yaml.safe_load(f)

build_matrix = wheel_config['jobs']['build_wheels']['strategy']['matrix']
n_python_versions = len(build_matrix['python'])

# For each python version we have: 5 wheels
# 1 osx wheel (x86_64)
# 2 linux wheel (i686 + x86_64)
# 2 windows wheel (win32 + wind_amd64)
n_wheels = 5 * n_python_versions

# plus one more for the sdist
n_wheels += 1

dist_files = list(Path("dist").glob('**/*'))
n_dist_files = len(dist_files)

if n_dist_files != n_wheels:
    print(f"Expected {n_wheels} wheels in dist/* but "
          f"got {n_dist_files} artifacts instead.")
    sys.exit(1)

print(f"dist/* has the expected {n_wheels} wheels:")
print("\n".join(file.name for file in dist_files))
