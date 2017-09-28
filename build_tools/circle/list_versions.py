#!/usr/bin/env python3

# List all available versions of the documentation

from urllib.request import urlopen
import sys
import json
import re

from distutils.version import LooseVersion


def json_urlread(url):
    return json.loads(urlopen(url).read().decode('utf8'))

heading = 'Available documentation for Scikit-learn'
print(heading)
print('=' * len(heading))
print()

ROOT_URL = 'https://api.github.com/repos/scikit-learn/scikit-learn.github.io/contents/'
RAW_FMT = 'https://raw.githubusercontent.com/scikit-learn/scikit-learn.github.io/master/%s/documentation.html'
VERSION_RE = re.compile(r"\bVERSION:\s*'([^']+)'")


def human_readable_data_quantity(quantity, multiple=1024):
    # from https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    if quantity == 0:
        quantity = +0
    SUFFIXES = ["B"] + [i + {1000: "B", 1024: "iB"}[multiple] for i in "KMGTPEZY"]
    for suffix in SUFFIXES:
        if quantity < multiple or suffix == SUFFIXES[-1]:
            if suffix == SUFFIXES[0]:
                return "%d%s" % (quantity, suffix)
            else:
                return "%.1f%s" % (quantity, suffix)
        else:
            quantity /= multiple


def get_pdf_size(version):
    api_url = ROOT_URL + '%s/_downloads' % version
    for path_details in json_urlread(api_url):
        if path_details['name'] == 'scikit-learn-docs.pdf':
            return human_readable_data_quantity(path_details['size'], 1000)


dirs = {}
symlinks = {}

root_listing = json_urlread(ROOT_URL)
for path_details in root_listing:
    name = path_details['name']
    if path_details['type'] == 'dir':
        try:
            html = urlopen(RAW_FMT % name).read().decode('utf8')
        except Exception:
            print('Failed to fetch %s' % (RAW_FMT % name), file=sys.stderr)
            continue
        version_num = VERSION_RE.search(html).group(1)
        pdf_size = get_pdf_size(name)
        dirs[name] = (version_num, pdf_size)

    if path_details['type'] == 'symlink':
        symlinks[name] = json_urlread(path_details['_links']['self'])['target']


for src, dst in symlinks.items():
    if dst in dirs:
        dirs[src] = dirs[dst]


digit_names = [k for k in dirs if k[:1].isdigit()]
word_names = [k for k in dirs if not k[:1].isdigit()]


seen = set()
for name in sorted(word_names) + sorted(digit_names, key=LooseVersion, reverse=True):
    version_num, pdf_size = dirs[name]
    if version_num in seen:
        continue
    else:
        seen.add(version_num)
    name_display = '' if name[:1].isdigit() else ' (%s)' % name
    out = '* `Scikit-learn %s%s documentation <http://scikit-learn.org/%s/documentation.html>`' % (version_num, name_display, name)
    if pdf_size is not None:
        out += ' (`PDF %s <http://scikit-learn.org/%s/_downloads/scikit-learn-docs.pdf>`)' % (pdf_size, name)
    print(out)
