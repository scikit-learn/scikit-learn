#!/usr/bin/env python3

# List all available versions of the documentation
import json
import re
import sys

from distutils.version import LooseVersion
from urllib.request import urlopen


def json_urlread(url):
    try:
        return json.loads(urlopen(url).read().decode("utf8"))
    except Exception:
        print("Error reading", url, file=sys.stderr)
        raise
