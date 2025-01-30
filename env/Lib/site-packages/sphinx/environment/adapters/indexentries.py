"""Index entries adapters for sphinx.environment."""

from __future__ import annotations

import re
import unicodedata
from itertools import groupby
from typing import TYPE_CHECKING

from sphinx.errors import NoUri
from sphinx.locale import _, __
from sphinx.util import logging
from sphinx.util.index_entries import _split_into

if TYPE_CHECKING:
    from typing import Literal, TypeAlias

    from sphinx.builders import Builder
    from sphinx.environment import BuildEnvironment

    _IndexEntryTarget: TypeAlias = tuple[str | None, str | Literal[False]]
    _IndexEntryTargets: TypeAlias = list[_IndexEntryTarget]
    _IndexEntryCategoryKey: TypeAlias = str | None
    _IndexEntrySubItems: TypeAlias = dict[
        str,
        tuple[_IndexEntryTargets, _IndexEntryCategoryKey],
    ]
    _IndexEntry: TypeAlias = tuple[
        _IndexEntryTargets,
        _IndexEntrySubItems,
        _IndexEntryCategoryKey,
    ]
    _IndexEntryMap: TypeAlias = dict[str, _IndexEntry]
    _Index: TypeAlias = list[
        tuple[
            str,
            list[
                tuple[
                    str,
                    tuple[
                        _IndexEntryTargets,
                        list[tuple[str, _IndexEntryTargets]],
                        _IndexEntryCategoryKey,
                    ],
                ]
            ],
        ]
    ]

logger = logging.getLogger(__name__)


class IndexEntries:
    def __init__(self, env: BuildEnvironment) -> None:
        self.env = env
        self.builder: Builder

    def create_index(
        self,
        builder: Builder,
        group_entries: bool = True,
        _fixre: re.Pattern[str] = re.compile(r'(.*) ([(][^()]*[)])'),
    ) -> _Index:
        """Create the real index from the collected index entries."""
        new: _IndexEntryMap = {}

        rel_uri: str | Literal[False]
        index_domain = self.env.domains.index_domain
        for docname, entries in index_domain.entries.items():
            try:
                rel_uri = builder.get_relative_uri('genindex', docname)
            except NoUri:
                rel_uri = False

            # new entry types must be listed in directives/other.py!
            for entry_type, value, target_id, main, category_key in entries:
                uri = rel_uri is not False and f'{rel_uri}#{target_id}'
                try:
                    if entry_type == 'single':
                        try:
                            entry, sub_entry = _split_into(2, 'single', value)
                        except ValueError:
                            (entry,) = _split_into(1, 'single', value)
                            sub_entry = ''
                        _add_entry(
                            entry, sub_entry, main, dic=new, link=uri, key=category_key
                        )
                    elif entry_type == 'pair':
                        first, second = _split_into(2, 'pair', value)
                        _add_entry(
                            first, second, main, dic=new, link=uri, key=category_key
                        )
                        _add_entry(
                            second, first, main, dic=new, link=uri, key=category_key
                        )
                    elif entry_type == 'triple':
                        first, second, third = _split_into(3, 'triple', value)
                        _add_entry(
                            first,
                            second + ' ' + third,
                            main,
                            dic=new,
                            link=uri,
                            key=category_key,
                        )
                        _add_entry(
                            second,
                            third + ', ' + first,
                            main,
                            dic=new,
                            link=uri,
                            key=category_key,
                        )
                        _add_entry(
                            third,
                            first + ' ' + second,
                            main,
                            dic=new,
                            link=uri,
                            key=category_key,
                        )
                    elif entry_type == 'see':
                        first, second = _split_into(2, 'see', value)
                        _add_entry(
                            first,
                            _('see %s') % second,
                            None,
                            dic=new,
                            link=False,
                            key=category_key,
                        )
                    elif entry_type == 'seealso':
                        first, second = _split_into(2, 'see', value)
                        _add_entry(
                            first,
                            _('see also %s') % second,
                            None,
                            dic=new,
                            link=False,
                            key=category_key,
                        )
                    else:
                        logger.warning(
                            __('unknown index entry type %r'),
                            entry_type,
                            location=docname,
                        )
                except ValueError as err:
                    logger.warning(str(err), location=docname)

        for targets, sub_items, _category_key in new.values():
            targets.sort(key=_key_func_0)
            for sub_targets, _sub_category_key in sub_items.values():
                sub_targets.sort(key=_key_func_0)

        new_list: list[tuple[str, _IndexEntry]] = sorted(new.items(), key=_key_func_1)

        if group_entries:
            # fixup entries: transform
            #   func() (in module foo)
            #   func() (in module bar)
            # into
            #   func()
            #     (in module foo)
            #     (in module bar)
            old_key = ''
            old_sub_items: _IndexEntrySubItems = {}
            i = 0
            while i < len(new_list):
                key, (targets, sub_items, category_key) = new_list[i]
                # cannot move if it has sub_items; structure gets too complex
                if not sub_items:
                    m = _fixre.match(key)
                    if m:
                        if old_key == m.group(1):
                            # prefixes match: add entry as subitem of the
                            # previous entry
                            prev = old_sub_items.setdefault(m[2], ([], category_key))
                            prev[0].extend(targets)
                            del new_list[i]
                            continue
                        old_key = m.group(1)
                    else:
                        old_key = key
                old_sub_items = sub_items
                i += 1

        grouped = []
        for group_key, group in groupby(new_list, _group_by_func):
            group_list = []
            for group_entry in group:
                entry_key, (targets, sub_items, category_key) = group_entry
                pairs = [
                    (sub_key, sub_targets)
                    for (sub_key, (sub_targets, _sub_category_key)) in sub_items.items()
                ]
                pairs.sort(key=_key_func_2)
                group_list.append((entry_key, (targets, pairs, category_key)))
            grouped.append((group_key, group_list))
        return grouped


def _add_entry(
    word: str,
    subword: str,
    main: str | None,
    *,
    dic: _IndexEntryMap,
    link: str | Literal[False],
    key: _IndexEntryCategoryKey,
) -> None:
    entry = dic.setdefault(word, ([], {}, key))
    if subword:
        targets = entry[1].setdefault(subword, ([], key))[0]
    else:
        targets = entry[0]
    if link:
        targets.append((main, link))


def _key_func_0(entry: _IndexEntryTarget) -> tuple[bool, str | Literal[False]]:
    """Sort the index entries for same keyword."""
    main, uri = entry
    return not main, uri  # show main entries at first


def _key_func_1(entry: tuple[str, _IndexEntry]) -> tuple[tuple[int, str], str]:
    """Sort the index entries"""
    key, (_targets, _sub_items, category_key) = entry
    if category_key:
        # using the specified category key to sort
        key = category_key
    lc_key = unicodedata.normalize('NFD', key.lower())
    if lc_key.startswith('\N{RIGHT-TO-LEFT MARK}'):
        lc_key = lc_key[1:]

    if not lc_key[0:1].isalpha() and not lc_key.startswith('_'):
        # put symbols at the front of the index (0)
        group = 0
    else:
        # put non-symbol characters at the following group (1)
        group = 1
    # ensure a deterministic order *within* letters by also sorting on
    # the entry itself
    return (group, lc_key), entry[0]


def _key_func_2(entry: tuple[str, _IndexEntryTargets]) -> str:
    """Sort the sub-index entries"""
    key = unicodedata.normalize('NFD', entry[0].lower())
    if key.startswith('\N{RIGHT-TO-LEFT MARK}'):
        key = key[1:]
    if key[0:1].isalpha() or key.startswith('_'):
        key = chr(127) + key
    return key


def _group_by_func(entry: tuple[str, _IndexEntry]) -> str:
    """Group the entries by letter or category key."""
    key, (targets, sub_items, category_key) = entry

    if category_key is not None:
        return category_key

    # now calculate the key
    if key.startswith('\N{RIGHT-TO-LEFT MARK}'):
        key = key[1:]
    letter = unicodedata.normalize('NFD', key[0])[0].upper()
    if letter.isalpha() or letter == '_':
        return letter

    # get all other symbols under one heading
    return _('Symbols')
