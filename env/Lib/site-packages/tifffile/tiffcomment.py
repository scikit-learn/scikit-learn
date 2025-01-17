#!/usr/bin/env python3
# tifffile/tiffcomment.py

"""Print or replace ImageDescription in first page of TIFF file.

Usage: ``tiffcomment [--set comment] file``

"""

from __future__ import annotations

import os
import sys

try:
    from .tifffile import tiffcomment
except ImportError:
    try:
        from tifffile.tifffile import tiffcomment
    except ImportError:
        from tifffile import tiffcomment


def main(argv: list[str] | None = None) -> int:
    """Tiffcomment command line usage main function."""
    comment: str | bytes | None

    if argv is None:
        argv = sys.argv

    if len(argv) > 2 and argv[1] in '--set':
        comment = argv[2]
        files = argv[3:]
    else:
        comment = None
        files = argv[1:]

    if len(files) == 0 or any(f.startswith('-') for f in files):
        print()
        print(__doc__.strip())
        return 0

    if comment is None:
        pass
    elif os.path.exists(comment):
        with open(comment, 'rb') as fh:
            comment = fh.read()
    else:
        try:
            comment = comment.encode('ascii')
        except UnicodeEncodeError as exc:
            print(f'{exc}')
            assert isinstance(comment, str)
            comment = comment.encode()

    for file in files:
        try:
            result = tiffcomment(file, comment)
        except Exception as exc:
            print(f'{file}: {exc}')
        else:
            if result:
                print(result)
    return 0


if __name__ == '__main__':
    sys.exit(main())
