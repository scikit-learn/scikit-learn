"""
This module contains factory functions that attempt
to return Qt submodules from the various python Qt bindings.

It also protects against double-importing Qt with different
bindings, which is unstable and likely to crash

This is used primarily by qt and qt_for_kernel, and shouldn't
be accessed directly from the outside
"""
import sys
import types
from functools import partial
from importlib import import_module

from IPython.utils.version import check_version

# Available APIs.
QT_API_PYQT = 'pyqt' # Force version 2
QT_API_PYQT5 = 'pyqt5'
QT_API_PYQTv1 = 'pyqtv1' # Force version 2
QT_API_PYQT_DEFAULT = 'pyqtdefault' # use system default for version 1 vs. 2
QT_API_PYSIDE = 'pyside'
QT_API_PYSIDE2 = 'pyside2'

api_to_module = {QT_API_PYSIDE2: 'PySide2',
                 QT_API_PYSIDE: 'PySide',
                 QT_API_PYQT: 'PyQt4',
                 QT_API_PYQTv1: 'PyQt4',
                 QT_API_PYQT5: 'PyQt5',
                 QT_API_PYQT_DEFAULT: 'PyQt4',
                }


class ImportDenier(object):
    """Import Hook that will guard against bad Qt imports
    once IPython commits to a specific binding
    """

    def __init__(self):
        self.__forbidden = set()

    def forbid(self, module_name):
        sys.modules.pop(module_name, None)
        self.__forbidden.add(module_name)

    def find_module(self, fullname, path=None):
        if path:
            return
        if fullname in self.__forbidden:
            return self

    def load_module(self, fullname):
        raise ImportError("""
    Importing %s disabled by IPython, which has
    already imported an Incompatible QT Binding: %s
    """ % (fullname, loaded_api()))

ID = ImportDenier()
sys.meta_path.insert(0, ID)


def commit_api(api):
    """Commit to a particular API, and trigger ImportErrors on subsequent
       dangerous imports"""

    if api == QT_API_PYSIDE2:
        ID.forbid('PySide')
        ID.forbid('PyQt4')
        ID.forbid('PyQt5')
    if api == QT_API_PYSIDE:
        ID.forbid('PySide2')
        ID.forbid('PyQt4')
        ID.forbid('PyQt5')
    elif api == QT_API_PYQT5:
        ID.forbid('PySide2')
        ID.forbid('PySide')
        ID.forbid('PyQt4')
    else:   # There are three other possibilities, all representing PyQt4
        ID.forbid('PyQt5')
        ID.forbid('PySide2')
        ID.forbid('PySide')


def loaded_api():
    """Return which API is loaded, if any

    If this returns anything besides None,
    importing any other Qt binding is unsafe.

    Returns
    -------
    None, 'pyside2', 'pyside', 'pyqt', 'pyqt5', or 'pyqtv1'
    """
    if 'PyQt4.QtCore' in sys.modules:
        if qtapi_version() == 2:
            return QT_API_PYQT
        else:
            return QT_API_PYQTv1
    elif 'PySide.QtCore' in sys.modules:
        return QT_API_PYSIDE
    elif 'PySide2.QtCore' in sys.modules:
        return QT_API_PYSIDE2
    elif 'PyQt5.QtCore' in sys.modules:
        return QT_API_PYQT5
    return None


def has_binding(api):
    """Safely check for PyQt4/5, PySide or PySide2, without importing submodules

    Supports Python <= 3.3

       Parameters
       ----------
       api : str [ 'pyqtv1' | 'pyqt' | 'pyqt5' | 'pyside' | 'pyside2' | 'pyqtdefault']
            Which module to check for

       Returns
       -------
       True if the relevant module appears to be importable
    """
    # we can't import an incomplete pyside and pyqt4
    # this will cause a crash in sip (#1431)
    # check for complete presence before importing
    module_name = api_to_module[api]

    import imp
    try:
        #importing top level PyQt4/PySide module is ok...
        mod = import_module(module_name)
        #...importing submodules is not
        imp.find_module('QtCore', mod.__path__)
        imp.find_module('QtGui', mod.__path__)
        imp.find_module('QtSvg', mod.__path__)
        if api in (QT_API_PYQT5, QT_API_PYSIDE2):
            # QT5 requires QtWidgets too
            imp.find_module('QtWidgets', mod.__path__)

        #we can also safely check PySide version
        if api == QT_API_PYSIDE:
            return check_version(mod.__version__, '1.0.3')
        else:
            return True
    except ImportError:
        return False

def has_binding_new(api):
    """Safely check for PyQt4/5, PySide or PySide2, without importing submodules

    Supports Python >= 3.4

        Parameters
        ----------
        api : str [ 'pyqtv1' | 'pyqt' | 'pyqt5' | 'pyside' | 'pyside2' | 'pyqtdefault']
             Which module to check for

        Returns
        -------
        True if the relevant module appears to be importable
     """
    module_name = api_to_module[api]
    from importlib.util import find_spec

    required = ['QtCore', 'QtGui', 'QtSvg']
    if api in (QT_API_PYQT5, QT_API_PYSIDE2):
        # QT5 requires QtWidgets too
        required.append('QtWidgets')

    for submod in required:
        try:
            spec = find_spec('%s.%s' % (module_name, submod))
        except ImportError:
            # Package (e.g. PyQt5) not found
            return False
        else:
            if spec is None:
                # Submodule (e.g. PyQt5.QtCore) not found
                return False

    if api == QT_API_PYSIDE:
        # We can also safely check PySide version
        import PySide
        return check_version(PySide.__version__, '1.0.3')

    return True

if sys.version_info >= (3, 4):
    has_binding = has_binding_new

def qtapi_version():
    """Return which QString API has been set, if any

    Returns
    -------
    The QString API version (1 or 2), or None if not set
    """
    try:
        import sip
    except ImportError:
        return
    try:
        return sip.getapi('QString')
    except ValueError:
        return


def can_import(api):
    """Safely query whether an API is importable, without importing it"""
    if not has_binding(api):
        return False

    current = loaded_api()
    if api == QT_API_PYQT_DEFAULT:
        return current in [QT_API_PYQT, QT_API_PYQTv1, None]
    else:
        return current in [api, None]


def import_pyqt4(version=2):
    """
    Import PyQt4

    Parameters
    ----------
    version : 1, 2, or None
      Which QString/QVariant API to use. Set to None to use the system
      default

    ImportErrors rasied within this function are non-recoverable
    """
    # The new-style string API (version=2) automatically
    # converts QStrings to Unicode Python strings. Also, automatically unpacks
    # QVariants to their underlying objects.
    import sip

    if version is not None:
        sip.setapi('QString', version)
        sip.setapi('QVariant', version)

    from PyQt4 import QtGui, QtCore, QtSvg

    if not check_version(QtCore.PYQT_VERSION_STR, '4.7'):
        raise ImportError("IPython requires PyQt4 >= 4.7, found %s" %
                          QtCore.PYQT_VERSION_STR)

    # Alias PyQt-specific functions for PySide compatibility.
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.Slot = QtCore.pyqtSlot

    # query for the API version (in case version == None)
    version = sip.getapi('QString')
    api = QT_API_PYQTv1 if version == 1 else QT_API_PYQT
    return QtCore, QtGui, QtSvg, api


def import_pyqt5():
    """
    Import PyQt5

    ImportErrors rasied within this function are non-recoverable
    """
    import sip

    from PyQt5 import QtCore, QtSvg, QtWidgets, QtGui

    # Alias PyQt-specific functions for PySide compatibility.
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.Slot = QtCore.pyqtSlot

    # Join QtGui and QtWidgets for Qt4 compatibility.
    QtGuiCompat = types.ModuleType('QtGuiCompat')
    QtGuiCompat.__dict__.update(QtGui.__dict__)
    QtGuiCompat.__dict__.update(QtWidgets.__dict__)

    api = QT_API_PYQT5
    return QtCore, QtGuiCompat, QtSvg, api


def import_pyside():
    """
    Import PySide

    ImportErrors raised within this function are non-recoverable
    """
    from PySide import QtGui, QtCore, QtSvg
    return QtCore, QtGui, QtSvg, QT_API_PYSIDE

def import_pyside2():
    """
    Import PySide2

    ImportErrors raised within this function are non-recoverable
    """
    from PySide2 import QtGui, QtCore, QtSvg, QtWidgets, QtPrintSupport

    # Join QtGui and QtWidgets for Qt4 compatibility.
    QtGuiCompat = types.ModuleType('QtGuiCompat')
    QtGuiCompat.__dict__.update(QtGui.__dict__)
    QtGuiCompat.__dict__.update(QtWidgets.__dict__)
    QtGuiCompat.__dict__.update(QtPrintSupport.__dict__)

    return QtCore, QtGuiCompat, QtSvg, QT_API_PYSIDE2


def load_qt(api_options):
    """
    Attempt to import Qt, given a preference list
    of permissible bindings

    It is safe to call this function multiple times.

    Parameters
    ----------
    api_options: List of strings
        The order of APIs to try. Valid items are 'pyside', 'pyside2',
        'pyqt', 'pyqt5', 'pyqtv1' and 'pyqtdefault'

    Returns
    -------

    A tuple of QtCore, QtGui, QtSvg, QT_API
    The first three are the Qt modules. The last is the
    string indicating which module was loaded.

    Raises
    ------
    ImportError, if it isn't possible to import any requested
    bindings (either becaues they aren't installed, or because
    an incompatible library has already been installed)
    """
    loaders = {
               QT_API_PYSIDE2: import_pyside2,
               QT_API_PYSIDE: import_pyside,
               QT_API_PYQT: import_pyqt4,
               QT_API_PYQT5: import_pyqt5,
               QT_API_PYQTv1: partial(import_pyqt4, version=1),
               QT_API_PYQT_DEFAULT: partial(import_pyqt4, version=None)
              }

    for api in api_options:

        if api not in loaders:
            raise RuntimeError(
                "Invalid Qt API %r, valid values are: %s" %
                (api, ", ".join(["%r" % k for k in loaders.keys()])))

        if not can_import(api):
            continue

        #cannot safely recover from an ImportError during this
        result = loaders[api]()
        api = result[-1]  # changed if api = QT_API_PYQT_DEFAULT
        commit_api(api)
        return result
    else:
        raise ImportError("""
    Could not load requested Qt binding. Please ensure that
    PyQt4 >= 4.7, PyQt5, PySide >= 1.0.3 or PySide2 is available,
    and only one is imported per session.

    Currently-imported Qt library:                              %r
    PyQt4 available (requires QtCore, QtGui, QtSvg):            %s
    PyQt5 available (requires QtCore, QtGui, QtSvg, QtWidgets): %s
    PySide >= 1.0.3 installed:                                  %s
    PySide2 installed:                                          %s
    Tried to load:                                              %r
    """ % (loaded_api(),
           has_binding(QT_API_PYQT),
           has_binding(QT_API_PYQT5),
           has_binding(QT_API_PYSIDE),
           has_binding(QT_API_PYSIDE2),
           api_options))
