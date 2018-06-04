# -*- coding: utf-8 -*-
from skimage.viewer import utils
from skimage.viewer.utils import dialogs
from skimage.viewer.qt import QtCore, QtWidgets, has_qt
from skimage._shared import testing


@testing.skipif(not has_qt, reason="Qt not installed")
def test_event_loop():
    utils.init_qtapp()
    timer = QtCore.QTimer()
    timer.singleShot(10, QtWidgets.QApplication.quit)
    utils.start_qtapp()


@testing.skipif(not has_qt, reason="Qt not installed")
def test_format_filename():
    fname = dialogs._format_filename(('apple', 2))
    assert fname == 'apple'
    fname = dialogs._format_filename('')
    assert fname is None


@testing.skipif(True, reason="Can't automatically close window. See #3081.")
@testing.skipif(not has_qt, reason="Qt not installed")
def test_open_file_dialog():
    QApp = utils.init_qtapp()
    timer = QtCore.QTimer()
    timer.singleShot(100, lambda: QApp.quit())
    filename = dialogs.open_file_dialog()
    assert filename is None


@testing.skipif(True, reason="Can't automatically close window. See #3081.")
@testing.skipif(not has_qt, reason="Qt not installed")
def test_save_file_dialog():
    QApp = utils.init_qtapp()
    timer = QtCore.QTimer()
    timer.singleShot(100, lambda: QApp.quit())
    filename = dialogs.save_file_dialog()
    assert filename is None
