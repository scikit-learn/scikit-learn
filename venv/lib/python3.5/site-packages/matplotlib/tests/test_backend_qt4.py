from __future__ import absolute_import, division, print_function

from matplotlib import pyplot as plt
from matplotlib._pylab_helpers import Gcf
import matplotlib
import copy

import pytest
try:
    # mock in python 3.3+
    from unittest import mock
except ImportError:
    import mock

with matplotlib.rc_context(rc={'backend': 'Qt4Agg'}):
    qt_compat = pytest.importorskip('matplotlib.backends.qt_compat')
from matplotlib.backends.backend_qt4 import (
    MODIFIER_KEYS, SUPER, ALT, CTRL, SHIFT)  # noqa

QtCore = qt_compat.QtCore
_, ControlModifier, ControlKey = MODIFIER_KEYS[CTRL]
_, AltModifier, AltKey = MODIFIER_KEYS[ALT]
_, SuperModifier, SuperKey = MODIFIER_KEYS[SUPER]
_, ShiftModifier, ShiftKey = MODIFIER_KEYS[SHIFT]

try:
    py_qt_ver = int(QtCore.PYQT_VERSION_STR.split('.')[0])
except AttributeError:
    py_qt_ver = QtCore.__version_info__[0]

if py_qt_ver != 4:
    pytestmark = pytest.mark.xfail(reason='Qt4 is not available')


@pytest.mark.backend('Qt4Agg')
def test_fig_close():
    # save the state of Gcf.figs
    init_figs = copy.copy(Gcf.figs)

    # make a figure using pyplot interface
    fig = plt.figure()

    # simulate user clicking the close button by reaching in
    # and calling close on the underlying Qt object
    fig.canvas.manager.window.close()

    # assert that we have removed the reference to the FigureManager
    # that got added by plt.figure()
    assert init_figs == Gcf.figs


@pytest.mark.parametrize(
    'qt_key, qt_mods, answer',
    [
        (QtCore.Qt.Key_A, ShiftModifier, 'A'),
        (QtCore.Qt.Key_A, QtCore.Qt.NoModifier, 'a'),
        (QtCore.Qt.Key_A, ControlModifier, 'ctrl+a'),
        (QtCore.Qt.Key_Aacute, ShiftModifier,
         '\N{LATIN CAPITAL LETTER A WITH ACUTE}'),
        (QtCore.Qt.Key_Aacute, QtCore.Qt.NoModifier,
         '\N{LATIN SMALL LETTER A WITH ACUTE}'),
        (ControlKey, AltModifier, 'alt+control'),
        (AltKey, ControlModifier, 'ctrl+alt'),
        (QtCore.Qt.Key_Aacute, (ControlModifier | AltModifier | SuperModifier),
         'ctrl+alt+super+\N{LATIN SMALL LETTER A WITH ACUTE}'),
        (QtCore.Qt.Key_Backspace, QtCore.Qt.NoModifier, 'backspace'),
        (QtCore.Qt.Key_Backspace, ControlModifier, 'ctrl+backspace'),
        (QtCore.Qt.Key_Play, QtCore.Qt.NoModifier, None),
    ],
    ids=[
        'shift',
        'lower',
        'control',
        'unicode_upper',
        'unicode_lower',
        'alt_control',
        'control_alt',
        'modifier_order',
        'backspace',
        'backspace_mod',
        'non_unicode_key',
    ]
)
@pytest.mark.backend('Qt4Agg')
def test_correct_key(qt_key, qt_mods, answer):
    """
    Make a figure
    Send a key_press_event event (using non-public, qt4 backend specific api)
    Catch the event
    Assert sent and caught keys are the same
    """
    qt_canvas = plt.figure().canvas

    event = mock.Mock()
    event.isAutoRepeat.return_value = False
    event.key.return_value = qt_key
    event.modifiers.return_value = qt_mods

    def receive(event):
        assert event.key == answer

    qt_canvas.mpl_connect('key_press_event', receive)
    qt_canvas.keyPressEvent(event)
