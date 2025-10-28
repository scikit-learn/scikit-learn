# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2021 The Meson development team
# Copyright © 2021-2023 Intel Corporation

from __future__ import annotations

"""base classes providing no-op functionality.."""

import enum
import os
import typing as T

from .. import mlog

__all__ = ['DirectoryLock', 'DirectoryLockAction', 'DirectoryLockBase']

class DirectoryLockAction(enum.Enum):
    IGNORE = 0
    WAIT = 1
    FAIL = 2

class DirectoryLockBase:
    def __init__(self, directory: str, lockfile: str, action: DirectoryLockAction, err: str,
                 optional: bool = False) -> None:
        self.action = action
        self.err = err
        self.lockpath = os.path.join(directory, lockfile)
        self.optional = optional
        self.lockfile: T.Optional[T.TextIO] = None

    def __enter__(self) -> None:
        mlog.debug('Calling the no-op version of DirectoryLock')

    def __exit__(self, *args: T.Any) -> None:
        pass

class DirectoryLock(DirectoryLockBase):
    pass
