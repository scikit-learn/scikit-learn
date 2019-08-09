#!/usr/bin/env python3

# List all available versions of the documentation
import json
import re
import sys

from distutils.version import LooseVersion
from urllib.request import urlopen

def json_urlread(url):
    try:
        return json.loads(urlopen(url).read().decode('utf8'))
    except Exception:
        print('Error reading', url, file=sys.stderr)
        raise


def human_readable_data_quantity(quantity, multiple=1024):
    # https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    if quantity == 0:
        quantity = +0
    SUFFIXES = ["B"] + [i + {1000: "B", 1024: "iB"}[multiple]
                        for i in "KMGTPEZY"]
    for suffix in SUFFIXES:
        if quantity < multiple or suffix == SUFFIXES[-1]:
            if suffix == SUFFIXES[0]:
                return "%d %s" % (quantity, suffix)
            else:
                return "%.1f %s" % (quantity, suffix)
        else:
            quantity /= multiple


def get_pdf_size(version):
    api_url = ROOT_URL + '%s/_downloads' % version
    for path_details in json_urlread(api_url):
        if path_details['name'] == 'scikit-learn-docs.pdf':
            return human_readable_data_quantity(path_details['size'], 1000)


print(':orphan:')
print()
heading = 'Available documentation for Scikit-learn'
print(heading)
print('=' * len(heading))
print()
print('Web-based documentation is available for versions listed below:')
print()

ROOT_URL = 'https://api.github.com/repos/scikit-learn/scikit-learn.github.io/contents/'  # noqa
RAW_FMT = 'https://raw.githubusercontent.com/scikit-learn/scikit-learn.github.io/master/%s/documentation.html'  # noqa
VERSION_RE = re.compile(r"scikit-learn ([\w\.\-]+) documentation</title>")
NAMED_DIRS = ['dev', 'stable']

# Gather data for each version directory, including symlinks
dirs = {}
symlinks = {}
root_listing = json_urlread(ROOT_URL)
for path_details in root_listing:
    name = path_details['name']
    if not (name[:1].isdigit() or name in NAMED_DIRS):
        continue
    if path_details['type'] == 'dir':
        html = urlopen(RAW_FMT % name).read().decode('utf8')
        version_num = VERSION_RE.search(html).group(1)
        pdf_size = get_pdf_size(name)
        dirs[name] = (version_num, pdf_size)

    if path_details['type'] == 'symlink':
        symlinks[name] = json_urlread(path_details['_links']['self'])['target']


# Symlinks should have same data as target
for src, dst in symlinks.items():
    if dst in dirs:
        dirs[src] = dirs[dst]

# Output in order: dev, stable, decreasing other version
seen = set()
for name in (NAMED_DIRS +
             sorted((k for k in dirs if k[:1].isdigit()),
                    key=LooseVersion, reverse=True)):
    version_num, pdf_size = dirs[name]
    if version_num in seen:
        # symlink came first
        continue
    else:
        seen.add(version_num)
    name_display = '' if name[:1].isdigit() else ' (%s)' % name
    path = 'http://scikit-learn.org/%s' % name
    out = ('* `Scikit-learn %s%s documentation <%s/documentation.html>`_'
           % (version_num, name_display, path))
    if pdf_size is not None:
        out += (' (`PDF %s <%s/_downloads/scikit-learn-docs.pdf>`_)'
                % (pdf_size, path))
    print(out)
