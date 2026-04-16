# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

"""Utility functions with platform specific implementations."""

from __future__ import annotations

import enum
import os
import sys
import typing as T

from .. import mlog
from .core import MesonException

__all__ = ['DirectoryLock', 'DirectoryLockAction', 'path_has_root']

class DirectoryLockAction(enum.Enum):
    IGNORE = 0
    WAIT = 1
    FAIL = 2

class DirectoryLockBase:

    lockfile: T.Optional[T.TextIO] = None

    def __init__(self, directory: str, lockfile: str, action: DirectoryLockAction, err: str,
                 optional: bool = False) -> None:
        self.action = action
        self.err = err
        self.lockpath = os.path.join(directory, lockfile)
        self.optional = optional

    def __enter__(self) -> None:
        mlog.debug('Calling the no-op version of DirectoryLock')

    def __exit__(self, *args: T.Any) -> None:
        pass


if sys.platform == 'win32':
    import msvcrt

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
                raise

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

    # os.path.isabs changed in Python 3.13 (wtf) and now excludes Windows
    # path with a root component but no drive.  However, checking
    # whether something is from outside the source tree, or is installed
    # outside the prefix fails this new test.
    def path_has_root(path: str) -> bool:
        # Check for root-relative paths (and also UNC paths, which is okay);
        # unnecessary on Python pre-3.13, but not a problem to do it anyway
        if path and path[0] in '/\\':
            return True
        # Check for paths including a drive
        return os.path.isabs(path)

else:
    import fcntl

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
                raise

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

    # on POSIX systems don't go through pathlib as it's slower
    path_has_root = os.path.isabs
