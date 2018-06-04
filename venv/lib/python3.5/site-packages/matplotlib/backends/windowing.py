"""
MS Windows-specific helper for the TkAgg backend.

With rcParams['tk.window_focus'] default of False, it is
effectively disabled.

It uses a tiny C++ extension module to access MS Win functions.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib import rcParams

try:
    if not rcParams['tk.window_focus']:
        raise ImportError
    from matplotlib._windowing import GetForegroundWindow, SetForegroundWindow
except ImportError:
    def GetForegroundWindow():
        return 0
    def SetForegroundWindow(hwnd):
        pass

class FocusManager(object):
    def __init__(self):
        self._shellWindow = GetForegroundWindow()

    def __del__(self):
        SetForegroundWindow(self._shellWindow)
