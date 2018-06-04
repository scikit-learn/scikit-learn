"""
GObject compatibility loader; supports ``gi`` and ``pgi``.

The binding selection rules are as follows:
- if ``gi`` has already been imported, use it; else
- if ``pgi`` has already been imported, use it; else
- if ``gi`` can be imported, use it; else
- if ``pgi`` can be imported, use it; else
- error out.

Thus, to force usage of PGI when both bindings are installed, import it first.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import importlib
import sys


if "gi" in sys.modules:
    import gi
elif "pgi" in sys.modules:
    import pgi as gi
else:
    try:
        import gi
    except ImportError:
        try:
            import pgi as gi
        except ImportError:
            raise ImportError("The Gtk3 backend requires PyGObject or pgi")


gi.require_version("Gtk", "3.0")
globals().update(
    {name:
     importlib.import_module("{}.repository.{}".format(gi.__name__, name))
     for name in ["GLib", "GObject", "Gtk", "Gdk"]})
