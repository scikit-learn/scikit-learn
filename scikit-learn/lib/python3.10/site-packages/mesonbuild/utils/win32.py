# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

"""Windows specific implementations of mesonlib functionality."""

import msvcrt
import typing as T

from .core import MesonException
from .platform import DirectoryLockBase, DirectoryLockAction

__all__ = ['DirectoryLock', 'DirectoryLockAction']

class DirectoryLock(DirectoryLockBase):

    def __enter__(self) -> None:
        try:
            self.lockfile = open(self.lockpath, 'w+', encoding='utf-8')
        except (FileNotFoundError, IsADirectoryError):
            # For FileNotFoundError, there is nothing to lock.
            # For IsADirectoryError, something is seriously wrong.
            raise
        except OSError:
            if self.action == DirectoryLockAction.IGNORE or self.optional:
                return

        try:
            mode = msvcrt.LK_LOCK
            if self.action != DirectoryLockAction.WAIT:
                mode = msvcrt.LK_NBLCK
            msvcrt.locking(self.lockfile.fileno(), mode, 1)
        except BlockingIOError:
            self.lockfile.close()
            if self.action == DirectoryLockAction.IGNORE:
                return
            raise MesonException(self.err)
        except PermissionError:
            self.lockfile.close()
            raise MesonException(self.err)

    def __exit__(self, *args: T.Any) -> None:
        if self.lockfile is None or self.lockfile.closed:
            return
        msvcrt.locking(self.lockfile.fileno(), msvcrt.LK_UNLCK, 1)
        self.lockfile.close()
