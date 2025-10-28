# SPDX-License-Identifier: Apache-2.0
# Copyright 2015 The Meson development team

from __future__ import annotations
import typing as T

from ._qt import QtBaseModule
from . import ModuleInfo

if T.TYPE_CHECKING:
    from ..interpreter import Interpreter


class Qt5Module(QtBaseModule):

    INFO = ModuleInfo('qt5')

    def __init__(self, interpreter: Interpreter):
        QtBaseModule.__init__(self, interpreter, qt_version=5)


def initialize(interp: Interpreter) -> Qt5Module:
    return Qt5Module(interp)
