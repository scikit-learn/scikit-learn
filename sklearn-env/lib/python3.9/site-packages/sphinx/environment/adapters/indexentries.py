"""
    sphinx.environment.adapters.indexentries
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Index entries adapters for sphinx.environment.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
import unicodedata
from itertools import groupby
from typing import Any, Dict, List, Pattern, Tuple, cast

from sphinx.builders import Builder
from sphinx.domains.index import IndexDomain
from sphinx.environment import BuildEnvironment
from sphinx.errors import NoUri
from sphinx.locale import _, __
from sphinx.util import logging, split_into

logger = logging.getLogger(__name__)


class IndexEntries:
    def __init__(self, env: BuildEnvironment) -> None:
        self.env = env

    def create_index(self, builder: Builder, group_entries: bool = True,
                     _fixre: Pattern = re.compile(r'(.*) ([(][^()]*[)])')
                     ) -> List[Tuple[str, List[Tuple[str, Any]]]]:
        """Create the real index from the collected index entries."""
        new: Dict[str, List] = {}

        def add_entry(word: str, subword: str, main: str, link: bool = True,
                      dic: Dict = new, key: str = None) -> None:
            # Force the word to be unicode if it's a ASCII bytestring.
            # This will solve problems with unicode normalization later.
            # For instance the RFC role will add bytestrings at the moment
            word = str(word)
            entry = dic.get(word)
            if not entry:
                dic[word] = entry = [[], {}, key]
            if subword:
                add_entry(subword, '', main, link=link, dic=entry[1], key=key)
            elif link:
                try:
                    uri = builder.get_relative_uri('genindex', fn) + '#' + tid
                except NoUri:
                    pass
                else:
                    entry[0].append((main, uri))

        domain = cast(IndexDomain, self.env.get_domain('index'))
        for fn, entries in domain.entries.items():
            # new entry types must be listed in directives/other.py!
            for type, value, tid, main, index_key in entries:  # noqa: B007
                try:
                    if type == 'single':
                        try:
                            entry, subentry = split_into(2, 'single', value)
                        except ValueError:
                            entry, = split_into(1, 'single', value)
                            subentry = ''
                        add_entry(entry, subentry, main, key=index_key)
                    elif type == 'pair':
                        first, second = split_into(2, 'pair', value)
                        add_entry(first, second, main, key=index_key)
                        add_entry(second, first, main, key=index_key)
                    elif type == 'triple':
                        first, second, third = split_into(3, 'triple', value)
                        add_entry(first, second + ' ' + third, main, key=index_key)
                        add_entry(second, third + ', ' + first, main, key=index_key)
                        add_entry(third, first + ' ' + second, main, key=index_key)
                    elif type == 'see':
                        first, second = split_into(2, 'see', value)
                        add_entry(first, _('see %s') % second, None,
                                  link=False, key=index_key)
                    elif type == 'seealso':
                        first, second = split_into(2, 'see', value)
                        add_entry(first, _('see also %s') % second, None,
                                  link=False, key=index_key)
                    else:
                        logger.warning(__('unknown index entry type %r'), type, location=fn)
                except ValueError as err:
                    logger.warning(str(err), location=fn)

        # sort the index entries for same keyword.
        def keyfunc0(entry: Tuple[str, str]) -> Tuple[bool, str]:
            main, uri = entry
            return (not main, uri)  # show main entries at first

        for indexentry in new.values():
            indexentry[0].sort(key=keyfunc0)
            for subentry in indexentry[1].values():
                subentry[0].sort(key=keyfunc0)  # type: ignore

        # sort the index entries
        def keyfunc(entry: Tuple[str, List]) -> Tuple[Tuple[int, str], str]:
            key, (void, void, category_key) = entry
            if category_key:
                # using specified category key to sort
                key = category_key
            lckey = unicodedata.normalize('NFD', key.lower())
            if lckey.startswith('\N{RIGHT-TO-LEFT MARK}'):
                lckey = lckey[1:]

            if lckey[0:1].isalpha() or lckey.startswith('_'):
                # put non-symbol characters at the following group (1)
                sortkey = (1, lckey)
            else:
                # put symbols at the front of the index (0)
                sortkey = (0, lckey)
            # ensure a deterministic order *within* letters by also sorting on
            # the entry itself
            return (sortkey, entry[0])
        newlist = sorted(new.items(), key=keyfunc)

        if group_entries:
            # fixup entries: transform
            #   func() (in module foo)
            #   func() (in module bar)
            # into
            #   func()
            #     (in module foo)
            #     (in module bar)
            oldkey = ''
            oldsubitems: Dict[str, List] = None
            i = 0
            while i < len(newlist):
                key, (targets, subitems, _key) = newlist[i]
                # cannot move if it has subitems; structure gets too complex
                if not subitems:
                    m = _fixre.match(key)
                    if m:
                        if oldkey == m.group(1):
                            # prefixes match: add entry as subitem of the
                            # previous entry
                            oldsubitems.setdefault(m.group(2), [[], {}, _key])[0].\
                                extend(targets)
                            del newlist[i]
                            continue
                        oldkey = m.group(1)
                    else:
                        oldkey = key
                oldsubitems = subitems
                i += 1

        # sort the sub-index entries
        def keyfunc2(entry: Tuple[str, List]) -> str:
            key = unicodedata.normalize('NFD', entry[0].lower())
            if key.startswith('\N{RIGHT-TO-LEFT MARK}'):
                key = key[1:]
            if key[0:1].isalpha() or key.startswith('_'):
                key = chr(127) + key
            return key

        # group the entries by letter
        def keyfunc3(item: Tuple[str, List]) -> str:
            # hack: mutating the subitems dicts to a list in the keyfunc
            k, v = item
            v[1] = sorted(((si, se) for (si, (se, void, void)) in v[1].items()),
                          key=keyfunc2)
            if v[2] is None:
                # now calculate the key
                if k.startswith('\N{RIGHT-TO-LEFT MARK}'):
                    k = k[1:]
                letter = unicodedata.normalize('NFD', k[0])[0].upper()
                if letter.isalpha() or letter == '_':
                    return letter
                else:
                    # get all other symbols under one heading
                    return _('Symbols')
            else:
                return v[2]
        return [(key_, list(group))
                for (key_, group) in groupby(newlist, keyfunc3)]
