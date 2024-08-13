#!/usr/bin/env python3

# Write the available versions page (--rst) and the version switcher JSON (--json).
# Version switcher see:
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/version-dropdown.html
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/announcements.html#announcement-banners

import argparse
import json
import re
import sys
from urllib.request import urlopen

from sklearn.utils.fixes import parse_version


def json_urlread(url):
    try:
        return json.loads(urlopen(url).read().decode("utf8"))
    except Exception:
        print("Error reading", url, file=sys.stderr)
        raise


def human_readable_data_quantity(quantity, multiple=1024):
    # https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    if quantity == 0:
        quantity = +0
    SUFFIXES = ["B"] + [i + {1000: "B", 1024: "iB"}[multiple] for i in "KMGTPEZY"]
    for suffix in SUFFIXES:
        if quantity < multiple or suffix == SUFFIXES[-1]:
            if suffix == SUFFIXES[0]:
                return "%d %s" % (quantity, suffix)
            else:
                return "%.1f %s" % (quantity, suffix)
        else:
            quantity /= multiple


def get_file_extension(version):
    if "dev" in version:
        # The 'dev' branch should be explicitly handled
        return "zip"

    current_version = parse_version(version)
    min_zip_version = parse_version("0.24")

    return "zip" if current_version >= min_zip_version else "pdf"


def get_file_size(version):
    api_url = ROOT_URL + "%s/_downloads" % version
    for path_details in json_urlread(api_url):
        file_extension = get_file_extension(version)
        file_path = f"scikit-learn-docs.{file_extension}"
        if path_details["name"] == file_path:
            return human_readable_data_quantity(path_details["size"], 1000)


parser = argparse.ArgumentParser()
parser.add_argument("--rst", type=str, required=True)
parser.add_argument("--json", type=str, required=True)
args = parser.parse_args()

heading = "Available documentation for scikit-learn"
json_content = []
rst_content = [
    ":orphan:\n",
    heading,
    "=" * len(heading) + "\n",
    "Web-based documentation is available for versions listed below:\n",
]

ROOT_URL = (
    "https://api.github.com/repos/scikit-learn/scikit-learn.github.io/contents/"  # noqa
)
RAW_FMT = "https://raw.githubusercontent.com/scikit-learn/scikit-learn.github.io/master/%s/index.html"  # noqa
VERSION_RE = re.compile(r"scikit-learn ([\w\.\-]+) documentation</title>")
NAMED_DIRS = ["dev", "stable"]

# Gather data for each version directory, including symlinks
dirs = {}
symlinks = {}
root_listing = json_urlread(ROOT_URL)
for path_details in root_listing:
    name = path_details["name"]
    if not (name[:1].isdigit() or name in NAMED_DIRS):
        continue
    if path_details["type"] == "dir":
        html = urlopen(RAW_FMT % name).read().decode("utf8")
        version_num = VERSION_RE.search(html).group(1)
        file_size = get_file_size(name)
        dirs[name] = (version_num, file_size)

    if path_details["type"] == "symlink":
        symlinks[name] = json_urlread(path_details["_links"]["self"])["target"]


# Symlinks should have same data as target
for src, dst in symlinks.items():
    if dst in dirs:
        dirs[src] = dirs[dst]

# Output in order: dev, stable, decreasing other version
seen = set()
for i, name in enumerate(
    NAMED_DIRS
    + sorted((k for k in dirs if k[:1].isdigit()), key=parse_version, reverse=True)
):
    version_num, file_size = dirs[name]
    if version_num in seen:
        # symlink came first
        continue
    else:
        seen.add(version_num)

    full_name = f"{version_num}" if name[:1].isdigit() else f"{version_num} ({name})"
    path = f"https://scikit-learn.org/{name}/"

    # Update JSON for the version switcher; only keep the 8 latest versions to avoid
    # overloading the version switcher dropdown
    if i < 8:
        info = {"name": full_name, "version": version_num, "url": path}
        if name == "stable":
            info["preferred"] = True
        json_content.append(info)

    # Printout for the historical version page
    out = f"* `scikit-learn {full_name} documentation <{path}>`_"
    if file_size is not None:
        file_extension = get_file_extension(version_num)
        out += (
            f" (`{file_extension.upper()} {file_size} <{path}/"
            f"_downloads/scikit-learn-docs.{file_extension}>`_)"
        )
    rst_content.append(out)

with open(args.rst, "w", encoding="utf-8") as f:
    f.write("\n".join(rst_content) + "\n")
print(f"Written {args.rst}")

with open(args.json, "w", encoding="utf-8") as f:
    json.dump(json_content, f, indent=2)
print(f"Written {args.json}")
