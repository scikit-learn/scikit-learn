# encoding: utf-8
"""
Enable Gtk3 to be used interacive by IPython.

Authors: Thomi Richards
"""
#-----------------------------------------------------------------------------
# Copyright (c) 2012, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import sys
from gi.repository import Gtk, GLib

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def _main_quit(*args, **kwargs):
    Gtk.main_quit()
    return False


def inputhook_gtk3():
    GLib.io_add_watch(sys.stdin, GLib.PRIORITY_DEFAULT, GLib.IO_IN, _main_quit)
    Gtk.main()
    return 0
