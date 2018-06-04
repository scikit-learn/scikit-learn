""" A Qt API selector that can be used to switch between PyQt and PySide.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import os
import logging
import sys
from matplotlib import rcParams

_log = logging.getLogger(__name__)

# Available APIs.
QT_API_PYQT = 'PyQt4'       # API is not set here; Python 2.x default is V 1
QT_API_PYQTv2 = 'PyQt4v2'   # forced to Version 2 API
QT_API_PYSIDE = 'PySide'    # only supports Version 2 API
QT_API_PYQT5 = 'PyQt5'      # use PyQt5 API; Version 2 with module shim
QT_API_PYSIDE2 = 'PySide2'  # Version 2 API with module shim

ETS = dict(pyqt=(QT_API_PYQTv2, 4), pyside=(QT_API_PYSIDE, 4),
           pyqt5=(QT_API_PYQT5, 5), pyside2=(QT_API_PYSIDE2, 5))
# ETS is a dict of env variable to (QT_API, QT_MAJOR_VERSION)
# If the ETS QT_API environment variable is set, use it, but only
# if the varible if of the same major QT version.  Note that
# ETS requires the version 2 of PyQt4, which is not the platform
# default for Python 2.x.

QT_API_ENV = os.environ.get('QT_API')

if rcParams['backend'] == 'Qt5Agg':
    QT_RC_MAJOR_VERSION = 5
elif rcParams['backend'] == 'Qt4Agg':
    QT_RC_MAJOR_VERSION = 4
else:
    # A different backend was specified, but we still got here because a Qt
    # related file was imported. This is allowed, so lets try and guess
    # what we should be using.
    if "PyQt4" in sys.modules or "PySide" in sys.modules:
        # PyQt4 or PySide is actually used.
        QT_RC_MAJOR_VERSION = 4
    else:
        # This is a fallback: PyQt5
        QT_RC_MAJOR_VERSION = 5

QT_API = None

# check if any binding is already imported, if so silently ignore the
# rcparams/ENV settings and use what ever is already imported.
if 'PySide' in sys.modules:
    # user has imported PySide before importing mpl
    QT_API = QT_API_PYSIDE

if 'PySide2' in sys.modules:
    # user has imported PySide before importing mpl
    QT_API = QT_API_PYSIDE2

if 'PyQt4' in sys.modules:
    # user has imported PyQt4 before importing mpl
    # this case also handles the PyQt4v2 case as once sip is imported
    # the API versions can not be changed so do not try
    QT_API = QT_API_PYQT

if 'PyQt5' in sys.modules:
    # the user has imported PyQt5 before importing mpl
    QT_API = QT_API_PYQT5

if (QT_API_ENV is not None) and QT_API is None:
    try:
        QT_ENV_MAJOR_VERSION = ETS[QT_API_ENV][1]
    except KeyError:
        raise RuntimeError(
            ('Unrecognized environment variable %r, valid values are:'
             ' %r, %r, %r or %r'
             % (QT_API_ENV, 'pyqt', 'pyside', 'pyqt5', 'pyside2')))
    if QT_ENV_MAJOR_VERSION == QT_RC_MAJOR_VERSION:
        # Only if backend and env qt major version are
        # compatible use the env variable.
        QT_API = ETS[QT_API_ENV][0]

_fallback_to_qt4 = False
if QT_API is None:
    # No ETS environment or incompatible so use rcParams.
    if rcParams['backend'] == 'Qt5Agg':
        QT_API = QT_API_PYQT5
    elif rcParams['backend'] == 'Qt4Agg':
        QT_API = QT_API_PYQT
    else:
        # A non-Qt backend was specified, no version of the Qt
        # bindings is imported, but we still got here because a Qt
        # related file was imported. This is allowed, fall back to Qt5
        # using which ever binding the rparams ask for.
        _fallback_to_qt4 = True
        QT_API = QT_API_PYQT5

# We will define an appropriate wrapper for the differing versions
# of file dialog.
_getSaveFileName = None

# Flag to check if sip could be imported
_sip_imported = False

# Now perform the imports.
if QT_API in (QT_API_PYQT, QT_API_PYQTv2, QT_API_PYQT5):
    try:
        import sip
        _sip_imported = True
    except ImportError:
        # Try using PySide
        if QT_RC_MAJOR_VERSION == 5:
            QT_API = QT_API_PYSIDE2
        else:
            QT_API = QT_API_PYSIDE
        cond = ("Could not import sip; falling back on PySide\n"
                "in place of PyQt4 or PyQt5.\n")
        _log.info(cond)

if _sip_imported:
    if QT_API == QT_API_PYQTv2:
        if QT_API_ENV == 'pyqt':
            cond = ("Found 'QT_API=pyqt' environment variable. "
                    "Setting PyQt4 API accordingly.\n")
        else:
            cond = "PyQt API v2 specified."
        try:
            sip.setapi('QString', 2)
        except:
            res = 'QString API v2 specification failed. Defaulting to v1.'
            _log.info(cond + res)
            # condition has now been reported, no need to repeat it:
            cond = ""
        try:
            sip.setapi('QVariant', 2)
        except:
            res = 'QVariant API v2 specification failed. Defaulting to v1.'
            _log.info(cond + res)
    if QT_API == QT_API_PYQT5:
        try:
            from PyQt5 import QtCore, QtGui, QtWidgets
            _getSaveFileName = QtWidgets.QFileDialog.getSaveFileName
        except ImportError:
            if _fallback_to_qt4:
                # fell through, tried PyQt5, failed fall back to PyQt4
                QT_API = QT_API_PYQT
                QT_RC_MAJOR_VERSION = 4
            else:
                raise

    # needs to be if so we can re-test the value of QT_API which may
    # have been changed in the above if block
    if QT_API in [QT_API_PYQT, QT_API_PYQTv2]:  # PyQt4 API
        from PyQt4 import QtCore, QtGui

        try:
            if sip.getapi("QString") > 1:
                # Use new getSaveFileNameAndFilter()
                _getSaveFileName = QtGui.QFileDialog.getSaveFileNameAndFilter
            else:

                # Use old getSaveFileName()
                def _getSaveFileName(*args, **kwargs):
                    return (QtGui.QFileDialog.getSaveFileName(*args, **kwargs),
                            None)

        except (AttributeError, KeyError):

            # call to getapi() can fail in older versions of sip
            def _getSaveFileName(*args, **kwargs):
                return QtGui.QFileDialog.getSaveFileName(*args, **kwargs), None
    try:
        # Alias PyQt-specific functions for PySide compatibility.
        QtCore.Signal = QtCore.pyqtSignal
        try:
            QtCore.Slot = QtCore.pyqtSlot
        except AttributeError:
            # Not a perfect match but works in simple cases
            QtCore.Slot = QtCore.pyqtSignature

        QtCore.Property = QtCore.pyqtProperty
        __version__ = QtCore.PYQT_VERSION_STR
    except NameError:
        # QtCore did not get imported, fall back to pyside
        if QT_RC_MAJOR_VERSION == 5:
            QT_API = QT_API_PYSIDE2
        else:
            QT_API = QT_API_PYSIDE


if QT_API == QT_API_PYSIDE2:
    try:
        from PySide2 import QtCore, QtGui, QtWidgets, __version__
        _getSaveFileName = QtWidgets.QFileDialog.getSaveFileName
    except ImportError:
        # tried PySide2, failed, fall back to PySide
        QT_RC_MAJOR_VERSION = 4
        QT_API = QT_API_PYSIDE

if QT_API == QT_API_PYSIDE:  # try importing pyside
    try:
        from PySide import QtCore, QtGui, __version__, __version_info__
    except ImportError:
        raise ImportError(
            "Matplotlib qt-based backends require an external PyQt4, PyQt5,\n"
            "PySide or PySide2 package to be installed, but it was not found.")

    if __version_info__ < (1, 0, 3):
        raise ImportError(
            "Matplotlib backend_qt4 and backend_qt4agg require PySide >=1.0.3")

    _getSaveFileName = QtGui.QFileDialog.getSaveFileName


# Apply shim to Qt4 APIs to make them look like Qt5
if QT_API in (QT_API_PYQT, QT_API_PYQTv2, QT_API_PYSIDE):
    '''Import all used QtGui objects into QtWidgets

    Here I've opted to simple copy QtGui into QtWidgets as that
    achieves the same result as copying over the objects, and will
    continue to work if other objects are used.

    '''
    QtWidgets = QtGui


def is_pyqt5():
    return QT_API == QT_API_PYQT5
