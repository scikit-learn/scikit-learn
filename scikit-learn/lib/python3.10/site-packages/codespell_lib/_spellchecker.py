#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see
# https://www.gnu.org/licenses/old-licenses/gpl-2.0.html.
"""
Copyright (C) 2010-2011  Lucas De Marchi <lucas.de.marchi@gmail.com>
Copyright (C) 2011  ProFUSION embedded systems
"""

from typing import (
    Dict,
    Set,
)

# Pass all misspellings through this translation table to generate
# alternative misspellings and fixes.
alt_chars = (("'", "’"),)  # noqa: RUF001


class Misspelling:
    def __init__(self, data: str, fix: bool, reason: str) -> None:
        self.data = data
        self.fix = fix
        self.reason = reason


def add_misspelling(
    key: str,
    data: str,
    misspellings: Dict[str, Misspelling],
) -> None:
    data = data.strip()

    if "," in data:
        fix = False
        data, reason = data.rsplit(",", 1)
        reason = reason.lstrip()
    else:
        fix = True
        reason = ""

    misspellings[key] = Misspelling(data, fix, reason)


def build_dict(
    filename: str,
    misspellings: Dict[str, Misspelling],
    ignore_words: Set[str],
) -> None:
    with open(filename, encoding="utf-8") as f:
        translate_tables = [(x, str.maketrans(x, y)) for x, y in alt_chars]
        for line in f:
            [key, data] = line.split("->")
            # TODO: For now, convert both to lower.
            #       Someday we can maybe add support for fixing caps.
            key = key.lower()
            data = data.lower()
            if key not in ignore_words:
                add_misspelling(key, data, misspellings)
            # generate alternative misspellings/fixes
            for x, table in translate_tables:
                if x in key:
                    alt_key = key.translate(table)
                    alt_data = data.translate(table)
                    if alt_key not in ignore_words:
                        add_misspelling(alt_key, alt_data, misspellings)
