#!/usr/bin/env python
# Sorts what's new entries with per-module headings.
# Pass what's new entries on stdin.

import sys
import re
from collections import defaultdict

LABEL_ORDER = ['MajorFeature', 'Feature', 'Enhancement', 'Efficiency',
               'Fix', 'API']


def entry_sort_key(s):
    if s.startswith('- |'):
        return LABEL_ORDER.index(s.split('|')[1])
    else:
        return -1


# discard headings and other non-entry lines
text = ''.join(l for l in sys.stdin
               if l.startswith('- ') or l.startswith(' '))

bucketed = defaultdict(list)

for entry in re.split('\n(?=- )', text.strip()):
    modules = re.findall(r':(?:func|meth|mod|class):'
                         r'`(?:[^<`]*<|~)?(?:sklearn.)?([a-z]\w+)',
                         entry)
    modules = set(modules)
    if len(modules) > 1:
        key = 'Multiple modules'
    elif modules:
        key = ':mod:`sklearn.%s`' % next(iter(modules))
    else:
        key = 'Miscellaneous'
    bucketed[key].append(entry)
    entry = entry.strip() + '\n'

everything = []
for key, bucket in sorted(bucketed.items()):
    everything.append(key + '\n' + '.' * len(key))
    bucket.sort(key=entry_sort_key)
    everything.extend(bucket)
print('\n\n'.join(everything))
