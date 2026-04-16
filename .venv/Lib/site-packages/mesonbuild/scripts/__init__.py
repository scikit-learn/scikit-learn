# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

from pathlib import PurePath

def destdir_join(d1: str, d2: str) -> str:
    if not d1:
        return d2
    # c:\destdir + c:\prefix must produce c:\destdir\prefix
    return str(PurePath(d1, *PurePath(d2).parts[1:]))
