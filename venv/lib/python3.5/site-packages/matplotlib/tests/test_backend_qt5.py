from __future__ import absolute_import, division, print_function

import copy

import matplotlib
from matplotlib import pyplot as plt
from matplotlib._pylab_helpers import Gcf

from numpy.testing import assert_equal

import pytest
try:
    # mock in python 3.3+
    from unittest import mock
except ImportError:
    import mock

with matplotlib.rc_context(rc={'backend': 'Qt5Agg'}):
    qt_compat = pytest.importorskip('matplotlib.backends.qt_compat',
                                    minversion='5')
from matplotlib.backends.backend_qt5 import (
    MODIFIER_KEYS, SUPER, ALT, CTRL, SHIFT)  # noqa

QtCore = qt_compat.QtCore
_, ControlModifier, ControlKey = MODIFIER_KEYS[CTRL]
_, AltModifier, AltKey = MODIFIER_KEYS[ALT]
_, SuperModifier, SuperKey = MODIFIER_KEYS[SUPER]
_, ShiftModifier, ShiftKey = MODIFIER_KEYS[SHIFT]


@pytest.mark.backend('Qt5Agg')
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
@pytest.mark.backend('Qt5Agg')
def test_correct_key(qt_key, qt_mods, answer):
    """
    Make a figure
    Send a key_press_event event (using non-public, qt5 backend specific api)
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


@pytest.mark.backend('Qt5Agg')
def test_dpi_ratio_change():
    """
    Make sure that if _dpi_ratio changes, the figure dpi changes but the
    widget remains the same physical size.
    """

    prop = 'matplotlib.backends.backend_qt5.FigureCanvasQT._dpi_ratio'

    with mock.patch(prop, new_callable=mock.PropertyMock) as p:

        p.return_value = 3

        fig = plt.figure(figsize=(5, 2), dpi=120)
        qt_canvas = fig.canvas
        qt_canvas.show()

        from matplotlib.backends.backend_qt5 import qApp

        # Make sure the mocking worked
        assert qt_canvas._dpi_ratio == 3

        size = qt_canvas.size()

        qt_canvas.manager.show()
        qt_canvas.draw()
        qApp.processEvents()

        # The DPI and the renderer width/height change
        assert fig.dpi == 360
        assert qt_canvas.renderer.width == 1800
        assert qt_canvas.renderer.height == 720

        # The actual widget size and figure physical size don't change
        assert size.width() == 600
        assert size.height() == 240
        assert_equal(qt_canvas.get_width_height(), (600, 240))
        assert_equal(fig.get_size_inches(), (5, 2))

        p.return_value = 2

        assert qt_canvas._dpi_ratio == 2

        qt_canvas.draw()
        qApp.processEvents()
        # this second processEvents is required to fully run the draw.
        # On `update` we notice the DPI has changed and trigger a
        # resize event to refresh, the second processEvents is
        # required to process that and fully update the window sizes.
        qApp.processEvents()

        # The DPI and the renderer width/height change
        assert fig.dpi == 240
        assert qt_canvas.renderer.width == 1200
        assert qt_canvas.renderer.height == 480

        # The actual widget size and figure physical size don't change
        assert size.width() == 600
        assert size.height() == 240
        assert_equal(qt_canvas.get_width_height(), (600, 240))
        assert_equal(fig.get_size_inches(), (5, 2))


@pytest.mark.backend('Qt5Agg')
def test_subplottool():
    fig, ax = plt.subplots()
    with mock.patch(
            "matplotlib.backends.backend_qt5.SubplotToolQt.exec_",
            lambda self: None):
        fig.canvas.manager.toolbar.configure_subplots()


@pytest.mark.backend('Qt5Agg')
def test_figureoptions():
    fig, ax = plt.subplots()
    ax.plot([1, 2])
    ax.imshow([[1]])
    with mock.patch(
            "matplotlib.backends.qt_editor.formlayout.FormDialog.exec_",
            lambda self: None):
        fig.canvas.manager.toolbar.edit_parameters()
