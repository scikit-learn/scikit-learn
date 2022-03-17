import re
import subprocess
import sys
from pathlib import Path
import shlex
import json

from jinja2 import Environment


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


def remove(alist, to_remove):
    return [each for each in alist if each not in (to_remove)]


conda_build_metadata_dict = {
    "pylatest_conda_forge_mkl_linux-64": {
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    "pylatest_conda_forge_mkl_osx-64": {
        "platform": "osx-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies
        + ["ccache", "compilers", "llvm-openmp"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    "pylatest_conda_mkl_no_openmp": {
        "platform": "osx-64",
        "channel": "defaults",
        "conda_dependencies": common_dependencies + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    "pylatest_conda_forge_mkl_no_coverage": {
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies_without_coverage + ["ccache"],
        "package_constraints": {
            "blas": "[build=mkl]",
        },
    },
    "py38_conda_defaults_openblas": {
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
    "py38_conda_forge_openblas_ubuntu_1804": {
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": common_dependencies_without_coverage + ["ccache"],
        "package_constraints": {"python": "3.8", "blas": "[build=openblas]"},
    },
    "pylatest_pip_openblas_pandas": {
        "platform": "linux-64",
        "channel": "defaults",
        "conda_dependencies": ["python", "ccache"],
        "pip_dependencies": remove(common_dependencies, ["python", "blas"])
        + docstring_test_dependencies
        + ["lightgbm", "scikit-image"],
        "package_constraints": {"python": "3.9"},
    },
    "pylatest_pip_scipy_dev": {
        "platform": "linux-64",
        "channel": "defaults",
        "conda_dependencies": ["python", "ccache"],
        "pip_dependencies": remove(
            common_dependencies, ["python", "blas", "matplotlib", "pyamg"]
        )
        + docstring_test_dependencies,
    },
    "pypy3": {
        "platform": "linux-64",
        "channel": "conda-forge",
        "conda_dependencies": ["pypy"]
        + remove(common_dependencies_without_coverage, ["python", "pandas", "pillow"])
        + ["ccache"],
        "package_constraints": {"blas": "[build=openblas]"},
    },
}


def get_conda_environment_content(build_metadata):
    template = environment.from_string(
        """channels:
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
  {% endif %}"""
    )
    return template.render(build_metadata=build_metadata)


def write_conda_environment(build_name, build_metadata, folder_path):
    content = get_conda_environment_content(build_metadata)
    output_path = folder_path / f"{build_name}_environment.yml"
    output_path.write_text(content)


def write_all_conda_environments(build_metadata_dict, folder_path):
    for build_name, build_metadata in build_metadata_dict.items():
        write_conda_environment(build_name, build_metadata, folder_path)


def conda_lock(environment_path, lock_file_path, platform):
    command = (
        f"conda-lock lock --mamba --platform {platform} "
        f"--file {environment_path} --lockfile {lock_file_path}"
    )

    execute_command(shlex.split(command))


def create_conda_lock_file(build_name, build_metadata, folder_path):
    environment_path = folder_path / f"{build_name}_environment.yml"
    platform = build_metadata["platform"]
    lock_file_basename = build_name
    if not lock_file_basename.endswith(platform):
        lock_file_basename = f"{lock_file_basename}_{platform}"

    lock_file_path = folder_path / f"{lock_file_basename}_conda-lock.yml"
    conda_lock(environment_path, lock_file_path, platform)


def write_all_conda_lock_files(build_metadata_dict, folder_path):
    for build_name, build_metadata in build_metadata_dict.items():
        create_conda_lock_file(build_name, build_metadata, output_path)


def get_pip_requirements_content(build_metadata):
    template = environment.from_string(
        """{% for pip_dep in build_metadata['pip_dependencies'] %}
{{ pip_dep | get_package_with_constraint(build_metadata, uses_pip=True) }}
{% endfor %}"""
    )
    return template.render(build_metadata=build_metadata)


def write_pip_requirements(build_name, build_metadata, folder_path):
    content = get_pip_requirements_content(build_metadata)
    output_path = folder_path / f"{build_name}_requirements.txt"
    output_path.write_text(content)


def write_all_pip_requirements(build_metadata_dict, folder_path):
    for build_name, build_metadata in build_metadata_dict.items():
        write_pip_requirements(build_name, build_metadata, folder_path)


def pip_compile(pip_compile_path, requirements_path, lock_file_path):
    command = f"{pip_compile_path} {requirements_path} -o {lock_file_path}"
    execute_command(shlex.split(command))


def write_pip_lock_file(build_name, build_metadata, folder_path):
    python_version = build_metadata["python_version"]
    environment_name = f"pip-tools-python{python_version}"
    # create conda env with python with right version + pip-compile
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


def write_all_pip_lock_files(build_metadata_dict, folder_path):
    for build_name, build_metadata in build_metadata_dict.items():
        write_pip_lock_file(build_name, build_metadata, folder_path)


pip_build_metadata_dict = {
    "debian_atlas_32bit": {
        "pip_dependencies": ["cython", "joblib", "threadpoolctl", "pytest"],
        "package_constraints": {
            "joblib": "min",
            "threadpoolctl": "2.2.0",
            "pytest": "min",
        },
        # same Python version as in debian-32 build
        "python_version": "3.9.2",
    },
    "ubuntu_atlas": {
        "pip_dependencies": [
            "cython",
            "joblib",
            "threadpoolctl",
            "pytest",
            "pytest-xdist",
        ],
        "package_constraints": {"joblib": "min", "threadpoolctl": "min"},
        # same Python version as in Ubuntu 20.04 build
        "python_version": "3.8.2",
    },
}


if __name__ == "__main__":
    output_path = Path("build_tools/azure/")
    print("Writing conda environments")
    write_all_conda_environments(conda_build_metadata_dict, output_path)
    print("Writing conda lock files")
    write_all_conda_lock_files(conda_build_metadata_dict, output_path)

    print("Writing pip requirements")
    write_all_pip_requirements(pip_build_metadata_dict, output_path)
    print("Writing pip lock files")
    write_all_pip_lock_files(pip_build_metadata_dict, output_path)
