# SPDX-FileCopyrightText: 2021 Filipe Laíns <lains@riseup.net>
# SPDX-FileCopyrightText: 2021 Quansight, LLC
# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import contextlib
import gzip
import os
import tarfile
import typing

from typing import IO


if typing.TYPE_CHECKING:  # pragma: no cover
    from mesonpy._compat import Iterator, Path


@contextlib.contextmanager
def chdir(path: Path) -> Iterator[Path]:
    """Context manager helper to change the current working directory -- cd."""
    old_cwd = os.getcwd()
    os.chdir(os.fspath(path))
    try:
        yield path
    finally:
        os.chdir(old_cwd)


@contextlib.contextmanager
def create_targz(path: Path) -> Iterator[tarfile.TarFile]:
    """Opens a .tar.gz file in the file system for edition.."""

    os.makedirs(os.path.dirname(path), exist_ok=True)
    file = typing.cast(IO[bytes], gzip.GzipFile(
        path,
        mode='w',
        # Set the stream last modification time to 0.  This mimics
        # what 'git archive' does and makes the archives byte-for-byte
        # reproducible.
        mtime=0,
    ))
    tar = tarfile.TarFile(
        mode='w',
        fileobj=file,
        format=tarfile.PAX_FORMAT,  # changed in 3.8 to GNU
    )

    with contextlib.closing(file), tar:
        yield tar


def setup_windows_console() -> bool:
    from ctypes import byref, windll  # type: ignore
    from ctypes.wintypes import DWORD

    STD_OUTPUT_HANDLE = -11
    ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x04

    kernel = windll.kernel32
    stdout = kernel.GetStdHandle(STD_OUTPUT_HANDLE)
    mode = DWORD()

    if not kernel.GetConsoleMode(stdout, byref(mode)):
        return False

    if not kernel.SetConsoleMode(stdout, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING):
        return False

    return True
