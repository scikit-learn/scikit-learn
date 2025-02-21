# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

"""Windows specific implementations of mesonlib functionality."""

import msvcrt
import typing as T

from .core import MesonException
from .platform import BuildDirLock as BuildDirLockBase

__all__ = ['BuildDirLock']

class BuildDirLock(BuildDirLockBase):

    def __enter__(self) -> None:
        self.lockfile = open(self.lockfilename, 'w', encoding='utf-8')
        try:
            msvcrt.locking(self.lockfile.fileno(), msvcrt.LK_NBLCK, 1)
        except (BlockingIOError, PermissionError):
            self.lockfile.close()
            raise MesonException('Some other Meson process is already using this build directory. Exiting.')

    def __exit__(self, *args: T.Any) -> None:
        msvcrt.locking(self.lockfile.fileno(), msvcrt.LK_UNLCK, 1)
        self.lockfile.close()
