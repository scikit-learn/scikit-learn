# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

"""Posix specific implementations of mesonlib functionality."""

import fcntl
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
            flags = fcntl.LOCK_EX
            if self.action != DirectoryLockAction.WAIT:
                flags = flags | fcntl.LOCK_NB
            fcntl.flock(self.lockfile, flags)
        except BlockingIOError:
            self.lockfile.close()
            if self.action == DirectoryLockAction.IGNORE:
                return
            raise MesonException(self.err)
        except PermissionError:
            self.lockfile.close()
            raise MesonException(self.err)
        except OSError as e:
            self.lockfile.close()
            raise MesonException(f'Failed to lock directory {self.lockpath}: {e.strerror}')

    def __exit__(self, *args: T.Any) -> None:
        if self.lockfile is None or self.lockfile.closed:
            return
        fcntl.flock(self.lockfile, fcntl.LOCK_UN)
        self.lockfile.close()
