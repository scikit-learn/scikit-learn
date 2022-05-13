"""Script to update CI environment files and associated lock files.

To run it you need to be in the root folder of the scikit-learn repo:
python build_tools/azure/update_environments_and_lock_files.py

Two scenarios where this script can be useful:
- make sure that the latest versions of all the dependencies are used in the CI.
  We can run this script regularly and open a PR with the changes to the lock
  files. This workflow will eventually be automated with a bot in the future.
- bump minimum dependencies in sklearn/_min_dependencies.py. Running this
  script will update both the CI environment files and associated lock files.
  You can then open a PR with the changes.
- pin some packages to an older version by adding them to the
  default_package_constraints variable. This is useful when regressions are
  introduced in our dependencies, this has happened for example with pytest 7
  and coverage 6.3.

Environments are conda environment.yml or pip requirements.txt. Lock files are
conda-lock lock files or pip-compile requirements.txt.

pip requirements.txt are used when we install some dependencies (e.g. numpy and
scipy) with apt-get and the rest of the dependencies (e.g. pytest and joblib)
with pip.

To run this script you need:
- conda-lock. The version should match the one used in the CI in
  build_tools/azure/install.sh
- pip-tools
- jinja2

"""

import re
import subprocess
import sys
from pathlib import Path
import shlex
import json

from jinja2 import Environment

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
logger.addHandler(handler)


common_dependencies_without_coverage = [
    "python",
    "numpy",
    "blas",
    "scipy",
    "cython",
    "joblib",
    "threadpoolctl",
    "matplotlib",
    "pandas",
    "pyamg",
    "pytest",
    "pytest-xdist",
    "pillow",
]

common_dependencies = common_dependencies_without_coverage + [
    "codecov",
    "pytest-cov",
    "coverage",
]

docstring_test_dependencies = ["sphinx", "numpydoc"]

default_package_constraints = {
    # XXX: pytest is temporary pinned to 6.2.5 because pytest 7 causes CI
    # issues https://github.com/scikit-learn/scikit-learn/pull/22381
    "pytest": "6.2.5",
    # XXX: coverage is temporary pinned to 6.2 because 6.3 is not fork-safe
    # cf. https://github.com/nedbat/coveragepy/issues/1310
    "coverage": "6.2",
}


def remove_from(alist, to_remove):
    return [each for each in alist if each not in to_remove]


conda_build_metadata_list = [
    {
        "build_name": "pylatest_conda_forge_mkl_linux-64",
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    {
        "build_name": "pylatest_conda_forge_mkl_osx-64",
        "platform": "osx-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies
        + ["ccache", "compilers", "llvm-openmp"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    {
        "build_name": "pylatest_conda_mkl_no_openmp",
        "platform": "osx-64",
        "channel": "defaults",
        "conda_dependencies": common_dependencies + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    {
        "build_name": "pylatest_conda_forge_mkl_no_coverage",
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies_without_coverage + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    {
        "build_name": "py38_conda_defaults_openblas",
        "platform": "linux-64",
        "channel": "defaults",
        "conda_dependencies": common_dependencies + ["ccache"],
        "package_constraints": {
            "python": "3.8",
            "blas": "[build=openblas]",
            "numpy": "min",
            "scipy": "min",
            "matplotlib": "min",
            "threadpoolctl": "2.2.0",
        },
    },
    {
        "build_name": "py38_conda_forge_openblas_ubuntu_1804",
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies_without_coverage + ["ccache"],
        "package_constraints": {"python": "3.8", "blas": "[build=openblas]"},
    },
    {
        "build_name": "pylatest_pip_openblas_pandas",
        "platform": "linux-64",
        "channel": "defaults",
        "conda_dependencies": ["python", "ccache"],
        "pip_dependencies": remove_from(common_dependencies, ["python", "blas"])
        + docstring_test_dependencies
        + ["lightgbm", "scikit-image"],
        "package_constraints": {"python": "3.9"},
    },
    {
        "build_name": "pylatest_pip_scipy_dev",
        "platform": "linux-64",
        "channel": "defaults",
        "conda_dependencies": ["python", "ccache"],
        "pip_dependencies": remove_from(
            common_dependencies,
            [
                "python",
                "blas",
                "matplotlib",
                "pyamg",
                # all the dependencies below have a development version
                # installed in the CI, so they can be removed from the
                # environment.yml
                "numpy",
                "scipy",
                "pandas",
                "cython",
                "joblib",
                "pillow",
            ],
        )
        + docstring_test_dependencies
        # python-dateutil is a dependency of pandas and pandas is removed from
        # the environment.yml. Adding python-dateutil so it is pinned
        + ["python-dateutil"],
    },
    {
        "build_name": "pypy3",
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": ["pypy"]
        + remove_from(
            common_dependencies_without_coverage, ["python", "pandas", "pillow"]
        )
        + ["ccache"],
        "package_constraints": {"blas": "[build=openblas]"},
    },
]


pip_build_metadata_list = [
    {
        "build_name": "debian_atlas_32bit",
        "pip_dependencies": ["cython", "joblib", "threadpoolctl", "pytest"],
        "package_constraints": {
            "joblib": "min",
            "threadpoolctl": "2.2.0",
            "pytest": "min",
            # no pytest-xdist because it causes issue on 32bit
        },
        # same Python version as in debian-32 build
        "python_version": "3.9.2",
    },
    {
        "build_name": "ubuntu_atlas",
        "pip_dependencies": [
            "cython",
            "joblib",
            "threadpoolctl",
            "pytest",
            "pytest-xdist",
        ],
        "package_constraints": {"joblib": "min", "threadpoolctl": "min"},
        # Ubuntu 20.04 has 3.8.2 but only 3.8.5 is available for osx-arm64 on
        # conda-forge. Chosing 3.8.5 so that this script can be run locally on
        # osx-arm64 machines. This should not matter for pining versions with
        # pip-compile
        "python_version": "3.8.5",
    },
]


def execute_command(command_list):
    proc = subprocess.Popen(
        command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    out, err = proc.communicate()
    out, err = out.decode(), err.decode()

    if proc.returncode != 0:
        command_str = " ".join(command_list)
        raise RuntimeError(
            "Command exited with non-zero exit code.\n"
            "Exit code: {}\n"
            "Command:\n{}\n"
            "stdout:\n{}\n"
            "stderr:\n{}\n".format(proc.returncode, command_str, out, err)
        )
    return out


def get_package_with_constraint(package_name, build_metadata, uses_pip=False):
    build_package_constraints = build_metadata.get("package_constraints")
    if build_package_constraints is None:
        constraint = None
    else:
        constraint = build_package_constraints.get(package_name)

    constraint = constraint or default_package_constraints.get(package_name)

    if constraint is None:
        return package_name

    comment = ""
    if constraint == "min":
        constraint = execute_command(
            [sys.executable, "sklearn/_min_dependencies.py", package_name]
        ).strip()
        comment = "  # min"

    if re.match(r"\d[.\d]*", constraint):
        equality = "==" if uses_pip else "="
        constraint = equality + constraint

    return f"{package_name}{constraint}{comment}"


environment = Environment(trim_blocks=True, lstrip_blocks=True)
environment.filters["get_package_with_constraint"] = get_package_with_constraint


def get_conda_environment_content(build_metadata):
    template = environment.from_string(
        """
# DO NOT EDIT: this file is generated from the specification found in the
# following script to centralize the configuration for all Azure CI builds:
# build_tools/azure/update_environments_and_lock_files.py
channels:
  - {{ build_metadata['channel'] }}
dependencies:
  {% for conda_dep in build_metadata['conda_dependencies'] %}
  - {{ conda_dep | get_package_with_constraint(build_metadata) }}
  {% endfor %}
  {% if build_metadata['pip_dependencies'] %}
  - pip
  - pip:
  {% for pip_dep in build_metadata.get('pip_dependencies', []) %}
    - {{ pip_dep | get_package_with_constraint(build_metadata, uses_pip=True) }}
  {% endfor %}
  {% endif %}""".strip()
    )
    return template.render(build_metadata=build_metadata)


def write_conda_environment(build_metadata, folder_path):
    content = get_conda_environment_content(build_metadata)
    build_name = build_metadata["build_name"]
    output_path = folder_path / f"{build_name}_environment.yml"
    output_path.write_text(content)


def write_all_conda_environments(build_metadata_list, folder_path):
    for build_metadata in build_metadata_list:
        write_conda_environment(build_metadata, folder_path)


def conda_lock(environment_path, lock_file_path, platform):
    command = (
        f"conda-lock lock --mamba --kind explicit --platform {platform} "
        f"--file {environment_path} --filename-template {lock_file_path}"
    )

    logger.debug("conda-lock command: %s", command)
    execute_command(shlex.split(command))


def create_conda_lock_file(build_metadata, folder_path):
    build_name = build_metadata["build_name"]
    environment_path = folder_path / f"{build_name}_environment.yml"
    platform = build_metadata["platform"]
    lock_file_basename = build_name
    if not lock_file_basename.endswith(platform):
        lock_file_basename = f"{lock_file_basename}_{platform}"

    lock_file_path = folder_path / f"{lock_file_basename}_conda.lock"
    conda_lock(environment_path, lock_file_path, platform)


def write_all_conda_lock_files(build_metadata_list, folder_path):
    for build_metadata in build_metadata_list:
        logger.info(build_metadata["build_name"])
        create_conda_lock_file(build_metadata, folder_path)


def get_pip_requirements_content(build_metadata):
    template = environment.from_string(
        """
# DO NOT EDIT: this file is generated from the specification found in the
# following script to centralize the configuration for all Azure CI builds:
# build_tools/azure/update_environments_and_lock_files.py
{% for pip_dep in build_metadata['pip_dependencies'] %}
{{ pip_dep | get_package_with_constraint(build_metadata, uses_pip=True) }}
{% endfor %}""".strip()
    )
    return template.render(build_metadata=build_metadata)


def write_pip_requirements(build_metadata, folder_path):
    build_name = build_metadata["build_name"]
    content = get_pip_requirements_content(build_metadata)
    output_path = folder_path / f"{build_name}_requirements.txt"
    output_path.write_text(content)


def write_all_pip_requirements(build_metadata_list, folder_path):
    for build_metadata in build_metadata_list:
        logger.info(build_metadata["build_name"])
        write_pip_requirements(build_metadata, folder_path)


def pip_compile(pip_compile_path, requirements_path, lock_file_path):
    command = f"{pip_compile_path} --upgrade {requirements_path} -o {lock_file_path}"

    logger.debug("pip-compile command: %s", command)
    execute_command(shlex.split(command))


def write_pip_lock_file(build_metadata, folder_path):
    build_name = build_metadata["build_name"]
    python_version = build_metadata["python_version"]
    environment_name = f"pip-tools-python{python_version}"
    # To make sure that the Python used to create the pip lock file is the same
    # as the one used during the CI build where the lock file is used, we first
    # create a conda environment with the correct Python version and
    # pip-compile and run pip-compile in this environment
    command = (
        "conda create -c conda-forge -n"
        f" pip-tools-python{python_version} python={python_version} pip-tools -y"
    )
    execute_command(shlex.split(command))

    json_output = execute_command(shlex.split("conda info --json"))
    conda_info = json.loads(json_output)
    environment_folder = [
        each for each in conda_info["envs"] if each.endswith(environment_name)
    ][0]
    environment_path = Path(environment_folder)
    pip_compile_path = environment_path / "bin" / "pip-compile"

    requirement_path = folder_path / f"{build_name}_requirements.txt"
    lock_file_path = folder_path / f"{build_name}_lock.txt"
    pip_compile(pip_compile_path, requirement_path, lock_file_path)


def write_all_pip_lock_files(build_metadata_list, folder_path):
    for build_metadata in build_metadata_list:
        write_pip_lock_file(build_metadata, folder_path)


if __name__ == "__main__":
    output_path = Path("build_tools/azure/")
    logger.info("Writing conda environments")
    write_all_conda_environments(conda_build_metadata_list, output_path)
    logger.info("Writing conda lock files")
    write_all_conda_lock_files(conda_build_metadata_list, output_path)

    logger.info("Writing pip requirements")
    write_all_pip_requirements(pip_build_metadata_list, output_path)
    logger.info("Writing pip lock files")
    write_all_pip_lock_files(pip_build_metadata_list, output_path)
