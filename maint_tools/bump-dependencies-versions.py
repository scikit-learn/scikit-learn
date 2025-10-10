import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests
from packaging import version

df_list = pd.read_html("https://devguide.python.org/versions/")
df = pd.concat(df_list).astype({"Branch": str})
release_dates = {}
python_version_info = {
    version: release_date
    for version, release_date in zip(df["Branch"], df["First release"])
}
python_version_info = {
    version: pd.to_datetime(release_date)
    for version, release_date in python_version_info.items()
}


def get_min_version_with_wheel(package_name, python_version):
    # For compiled dependencies we want the oldest minor version that has
    # wheels for 'python_version'
    url = f"https://pypi.org/pypi/{package_name}/json"
    response = requests.get(url)
    if response.status_code != 200:
        return None

    data = response.json()
    releases = data["releases"]

    compatible_versions = []
    # We want only minor X.Y.0 and not bugfix X.Y.Z
    minor_releases = [
        (ver, release_info)
        for ver, release_info in releases.items()
        if re.match(r"^\d+\.\d+\.0$", ver)
    ]
    for ver, release_info in minor_releases:
        for file_info in release_info:
            if (
                file_info["packagetype"] == "bdist_wheel"
                and f"cp{python_version.replace('.', '')}" in file_info["filename"]
                and not file_info["yanked"]
            ):
                compatible_versions.append(ver)
                break

    if not compatible_versions:
        return None

    return min(compatible_versions, key=version.parse)


def get_min_python_version(scikit_learn_release_date_str="today"):
    # min Python version is the most recent Python release at least 3 years old
    # at the time of the scikit-learn release
    if scikit_learn_release_date_str == "today":
        scikit_learn_release_date = pd.to_datetime(datetime.now().date())
    else:
        scikit_learn_release_date = datetime.strptime(
            scikit_learn_release_date_str, "%Y-%m-%d"
        )
    version_and_releases = [
        {"python_version": python_version, "python_release_date": python_release_date}
        for python_version, python_release_date in python_version_info.items()
        if (scikit_learn_release_date - python_release_date).days > 365 * 3
    ]
    return max(version_and_releases, key=lambda each: each["python_release_date"])[
        "python_version"
    ]


def get_min_version_pure_python(package_name, scikit_learn_release_date_str="today"):
    # for pure Python dependencies we want the most recent minor release that
    # is at least 2 years old
    if scikit_learn_release_date_str == "today":
        scikit_learn_release_date = pd.to_datetime(datetime.now().date())
    else:
        scikit_learn_release_date = datetime.strptime(
            scikit_learn_release_date_str, "%Y-%m-%d"
        )

    url = f"https://pypi.org/pypi/{package_name}/json"
    response = requests.get(url)
    if response.status_code != 200:
        return None

    data = response.json()
    releases = data["releases"]

    compatible_versions = []
    # We want only minor X.Y.0 and not bugfix X.Y.Z
    releases = [
        (ver, release_info)
        for ver, release_info in releases.items()
        if re.match(r"^\d+\.\d+\.0$", ver)
    ]
    for ver, release_info in releases:
        for file_info in release_info:
            if (
                file_info["packagetype"] == "bdist_wheel"
                and not file_info["yanked"]
                and (
                    scikit_learn_release_date - pd.to_datetime(file_info["upload_time"])
                ).days
                > 365 * 2
            ):
                compatible_versions.append(ver)
                break

    if not compatible_versions:
        return None

    return max(compatible_versions, key=version.parse)


def get_current_dependencies_version(dep):
    return (
        subprocess.check_output([sys.executable, "sklearn/_min_dependencies.py", dep])
        .decode()
        .strip()
    )


def get_current_min_python_version():
    content = Path("pyproject.toml").read_text()
    min_python = re.findall(r'requires-python\s*=\s*">=(\d+\.\d+)"', content)[0]

    return min_python


def show_versions_update(scikit_learn_release_date="today"):
    future_versions = {"python": get_min_python_version(scikit_learn_release_date)}

    compiled_dependencies = ["numpy", "scipy", "pandas", "matplotlib", "pyamg"]
    future_versions.update(
        {
            dep: get_min_version_with_wheel(dep, future_versions["python"])
            for dep in compiled_dependencies
        }
    )

    pure_python_dependencies = ["joblib", "threadpoolctl"]
    future_versions.update(
        {
            dep: get_min_version_pure_python(dep, scikit_learn_release_date)
            for dep in pure_python_dependencies
        }
    )

    current_versions = {"python": get_current_min_python_version()}
    current_versions.update(
        {
            dep: get_current_dependencies_version(dep)
            for dep in compiled_dependencies + pure_python_dependencies
        }
    )

    print(f"For future release at date {scikit_learn_release_date}")
    for k in future_versions:
        if future_versions[k] != current_versions[k]:
            print(f"- {k}: {current_versions[k]} -> {future_versions[k]}")


if __name__ == "__main__":
    scikit_learn_release_date = sys.argv[1] if len(sys.argv) > 1 else "today"
    show_versions_update(scikit_learn_release_date)
