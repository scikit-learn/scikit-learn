"""Update dependencies pinned environment variables.

Requires conda-lock and ruamel.yaml to run. ruamel.yaml is used because handle comments.

```bash
pip install conda-lock ruamel.yaml
```

This script assumes that the platform environment file is in
`build_tools/azure/env_files` and will place the pinned environment file into
`build_tools/azure/lock_files`.

Usage:

python build_tools/ci_update_dependencies.py linux-64 py38_conda_defaults_openblas.yml
"""
import os
from pathlib import Path
import tempfile
import subprocess
import argparse
import sys

from ruamel.yaml import YAML
from sklearn._min_dependencies import dependent_packages  # noqa


parser = argparse.ArgumentParser(description="Update environment yaml")
parser.add_argument("platform", help="platform to generate locked environment file")
parser.add_argument("env_file", help="conda environment file")

args = parser.parse_args()
env_file = args.env_file

env_path = Path("build_tools") / "azure" / "env_files" / env_file
if not env_path.exists():
    print(f"{env_file} does not exist")
    sys.exit(1)

yaml = YAML()
with env_path.open("r") as f:
    env_yaml = yaml.load(f)


# Check that min dependencies are correct
dependencies = env_yaml["dependencies"]
dependencies_comments = dependencies.ca.items

incorrect_mins = []
for i, dependency in enumerate(dependencies):
    comment = dependencies_comments.get(i)
    if comment is not None and "min" in comment[0].value:
        package, version = dependency.split("=")
        min_version = dependent_packages[package][0]
        if min_version != version:
            incorrect_mins.append((package, version, min_version))

if incorrect_mins:
    print("Found incorrect minimum versions:")
    for package, version, expected_min in incorrect_mins:
        print(f"{package}={version} (expected {expected_min})")
    sys.exit(1)

# Use conda-lock to find pinned versions
platform = args.platform
with tempfile.TemporaryDirectory() as tmp_file:
    tmp_path = os.path.join(tmp_file, "lock_file.yml")

    subprocess.check_output(
        ["conda-lock", "-f", str(env_path), "-p", platform, "--lockfile", tmp_path]
    )
    with open(tmp_path) as f:
        lock_file_contents = f.read()
        lock_file_yaml = yaml.load(lock_file_contents)


# Create a environment.yml with pinned versions
name_to_package = {package["name"]: package for package in lock_file_yaml["package"]}
dependencies = env_yaml["dependencies"]
pip_dependencies = None
for i, dependency in enumerate(dependencies):
    if isinstance(dependency, dict) and "pip" in dependency:
        pip_dependencies = dependency["pip"]
        continue

    if dependency in name_to_package:
        version = name_to_package[dependency]["version"]
        dependencies[i] = f"{dependency}={version}"

if pip_dependencies is not None:
    for i, dependency in enumerate(pip_dependencies):
        if dependency in name_to_package:
            version = name_to_package[dependency]["version"]
            pip_dependencies[i] = f"{dependency}=={version}"


output_env_file = Path("build_tools") / "azure" / "lock_files" / env_file
with output_env_file.open("w") as f:
    yaml.dump(env_yaml, f)
