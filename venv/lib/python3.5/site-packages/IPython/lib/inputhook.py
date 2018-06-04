# coding: utf-8
"""
Deprecated since IPython 5.0

Inputhook management for GUI event loop integration.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

try:
    import ctypes
except ImportError:
    ctypes = None
except SystemError: # IronPython issue, 2/8/2014
    ctypes = None
import os
import platform
import sys
from distutils.version import LooseVersion as V

from warnings import warn


warn("`IPython.lib.inputhook` is deprecated since IPython 5.0 and will be removed in future versions.",
        DeprecationWarning, stacklevel=2)


#-----------------------------------------------------------------------------
# Constants
#-----------------------------------------------------------------------------

# Constants for identifying the GUI toolkits.
GUI_WX = 'wx'
GUI_QT = 'qt'
GUI_QT4 = 'qt4'
GUI_GTK = 'gtk'
GUI_TK = 'tk'
GUI_OSX = 'osx'
GUI_GLUT = 'glut'
GUI_PYGLET = 'pyglet'
GUI_GTK3 = 'gtk3'
GUI_NONE = 'none' # i.e. disable

#-----------------------------------------------------------------------------
# Utilities
#-----------------------------------------------------------------------------

def _stdin_ready_posix():
    """Return True if there's something to read on stdin (posix version)."""
    infds, outfds, erfds = select.select([sys.stdin],[],[],0)
    return bool(infds)

def _stdin_ready_nt():
    """Return True if there's something to read on stdin (nt version)."""
    return msvcrt.kbhit()

def _stdin_ready_other():
    """Return True, assuming there's something to read on stdin."""
    return True

def _use_appnope():
    """Should we use appnope for dealing with OS X app nap?

    Checks if we are on OS X 10.9 or greater.
    """
    return sys.platform == 'darwin' and V(platform.mac_ver()[0]) >= V('10.9')

def _ignore_CTRL_C_posix():
    """Ignore CTRL+C (SIGINT)."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _allow_CTRL_C_posix():
    """Take CTRL+C into account (SIGINT)."""
    signal.signal(signal.SIGINT, signal.default_int_handler)

def _ignore_CTRL_C_other():
    """Ignore CTRL+C (not implemented)."""
    pass

def _allow_CTRL_C_other():
    """Take CTRL+C into account (not implemented)."""
    pass

if os.name == 'posix':
    import select
    import signal
    stdin_ready = _stdin_ready_posix
    ignore_CTRL_C = _ignore_CTRL_C_posix
    allow_CTRL_C = _allow_CTRL_C_posix
elif os.name == 'nt':
    import msvcrt
    stdin_ready = _stdin_ready_nt
    ignore_CTRL_C = _ignore_CTRL_C_other
    allow_CTRL_C = _allow_CTRL_C_other
else:
    stdin_ready = _stdin_ready_other
    ignore_CTRL_C = _ignore_CTRL_C_other
    allow_CTRL_C = _allow_CTRL_C_other


#-----------------------------------------------------------------------------
# Main InputHookManager class
#-----------------------------------------------------------------------------


class InputHookManager(object):
    """DEPRECATED since IPython 5.0

    Manage PyOS_InputHook for different GUI toolkits.

    This class installs various hooks under ``PyOSInputHook`` to handle
    GUI event loop integration.
    """
    
    def __init__(self):
        if ctypes is None:
            warn("IPython GUI event loop requires ctypes, %gui will not be available")
        else:
            self.PYFUNC = ctypes.PYFUNCTYPE(ctypes.c_int)
        self.guihooks = {}
        self.aliases = {}
        self.apps = {}
        self._reset()

    def _reset(self):
        self._callback_pyfunctype = None
        self._callback = None
        self._installed = False
        self._current_gui = None

    def get_pyos_inputhook(self):
        """DEPRECATED since IPython 5.0

        Return the current PyOS_InputHook as a ctypes.c_void_p."""
        warn("`get_pyos_inputhook` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        return ctypes.c_void_p.in_dll(ctypes.pythonapi,"PyOS_InputHook")

    def get_pyos_inputhook_as_func(self):
        """DEPRECATED since IPython 5.0

        Return the current PyOS_InputHook as a ctypes.PYFUNCYPE."""
        warn("`get_pyos_inputhook_as_func` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        return self.PYFUNC.in_dll(ctypes.pythonapi,"PyOS_InputHook")

    def set_inputhook(self, callback):
        """DEPRECATED since IPython 5.0

        Set PyOS_InputHook to callback and return the previous one."""
        # On platforms with 'readline' support, it's all too likely to
        # have a KeyboardInterrupt signal delivered *even before* an
        # initial ``try:`` clause in the callback can be executed, so
        # we need to disable CTRL+C in this situation.
        ignore_CTRL_C()
        self._callback = callback
        self._callback_pyfunctype = self.PYFUNC(callback)
        pyos_inputhook_ptr = self.get_pyos_inputhook()
        original = self.get_pyos_inputhook_as_func()
        pyos_inputhook_ptr.value = \
            ctypes.cast(self._callback_pyfunctype, ctypes.c_void_p).value
        self._installed = True
        return original

    def clear_inputhook(self, app=None):
        """DEPRECATED since IPython 5.0

        Set PyOS_InputHook to NULL and return the previous one.

        Parameters
        ----------
        app : optional, ignored
          This parameter is allowed only so that clear_inputhook() can be
          called with a similar interface as all the ``enable_*`` methods.  But
          the actual value of the parameter is ignored.  This uniform interface
          makes it easier to have user-level entry points in the main IPython
          app like :meth:`enable_gui`."""
        warn("`clear_inputhook` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        pyos_inputhook_ptr = self.get_pyos_inputhook()
        original = self.get_pyos_inputhook_as_func()
        pyos_inputhook_ptr.value = ctypes.c_void_p(None).value
        allow_CTRL_C()
        self._reset()
        return original

    def clear_app_refs(self, gui=None):
        """DEPRECATED since IPython 5.0

        Clear IPython's internal reference to an application instance.

        Whenever we create an app for a user on qt4 or wx, we hold a
        reference to the app.  This is needed because in some cases bad things
        can happen if a user doesn't hold a reference themselves.  This
        method is provided to clear the references we are holding.

        Parameters
        ----------
        gui : None or str
            If None, clear all app references.  If ('wx', 'qt4') clear
            the app for that toolkit.  References are not held for gtk or tk
            as those toolkits don't have the notion of an app.
        """
        warn("`clear_app_refs` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        if gui is None:
            self.apps = {}
        elif gui in self.apps:
            del self.apps[gui]

    def register(self, toolkitname, *aliases):
        """DEPRECATED since IPython 5.0

        Register a class to provide the event loop for a given GUI.
        
        This is intended to be used as a class decorator. It should be passed
        the names with which to register this GUI integration. The classes
        themselves should subclass :class:`InputHookBase`.
        
        ::
        
            @inputhook_manager.register('qt')
            class QtInputHook(InputHookBase):
                def enable(self, app=None):
                    ...
        """
        warn("`register` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        def decorator(cls):
            if ctypes is not None:
                inst = cls(self)
                self.guihooks[toolkitname] = inst
                for a in aliases:
                    self.aliases[a] = toolkitname
            return cls
        return decorator

    def current_gui(self):
        """DEPRECATED since IPython 5.0

        Return a string indicating the currently active GUI or None."""
        warn("`current_gui` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        return self._current_gui

    def enable_gui(self, gui=None, app=None):
        """DEPRECATED since IPython 5.0

        Switch amongst GUI input hooks by name.

        This is a higher level method than :meth:`set_inputhook` - it uses the
        GUI name to look up a registered object which enables the input hook
        for that GUI.

        Parameters
        ----------
        gui : optional, string or None
          If None (or 'none'), clears input hook, otherwise it must be one
          of the recognized GUI names (see ``GUI_*`` constants in module).

        app : optional, existing application object.
          For toolkits that have the concept of a global app, you can supply an
          existing one.  If not given, the toolkit will be probed for one, and if
          none is found, a new one will be created.  Note that GTK does not have
          this concept, and passing an app if ``gui=="GTK"`` will raise an error.

        Returns
        -------
        The output of the underlying gui switch routine, typically the actual
        PyOS_InputHook wrapper object or the GUI toolkit app created, if there was
        one.
        """
        warn("`enable_gui` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        if gui in (None, GUI_NONE):
            return self.disable_gui()
        
        if gui in self.aliases:
            return self.enable_gui(self.aliases[gui], app)
        
        try:
            gui_hook = self.guihooks[gui]
        except KeyError:
            e = "Invalid GUI request {!r}, valid ones are: {}"
            raise ValueError(e.format(gui, ', '.join(self.guihooks)))
        self._current_gui = gui

        app = gui_hook.enable(app)
        if app is not None:
            app._in_event_loop = True
            self.apps[gui] = app
        return app

    def disable_gui(self):
        """DEPRECATED since IPython 5.0

        Disable GUI event loop integration.
        
        If an application was registered, this sets its ``_in_event_loop``
        attribute to False. It then calls :meth:`clear_inputhook`.
        """
        warn("`disable_gui` is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        gui = self._current_gui
        if gui in self.apps:
            self.apps[gui]._in_event_loop = False
        return self.clear_inputhook()

class InputHookBase(object):
    """DEPRECATED since IPython 5.0

    Base class for input hooks for specific toolkits.
    
    Subclasses should define an :meth:`enable` method with one argument, ``app``,
    which will either be an instance of the toolkit's application class, or None.
    They may also define a :meth:`disable` method with no arguments.
    """
    def __init__(self, manager):
        self.manager = manager

    def disable(self):
        pass

inputhook_manager = InputHookManager()

@inputhook_manager.register('osx')
class NullInputHook(InputHookBase):
    """DEPRECATED since IPython 5.0

    A null inputhook that doesn't need to do anything"""
    def enable(self, app=None):
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)

@inputhook_manager.register('wx')
class WxInputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with wxPython.

        Parameters
        ----------
        app : WX Application, optional.
            Running application to use.  If not given, we probe WX for an
            existing application object, and create a new one if none is found.

        Notes
        -----
        This methods sets the ``PyOS_InputHook`` for wxPython, which allows
        the wxPython to integrate with terminal based applications like
        IPython.

        If ``app`` is not given we probe for an existing one, and return it if
        found.  If no existing app is found, we create an :class:`wx.App` as
        follows::

            import wx
            app = wx.App(redirect=False, clearSigInt=False)
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        import wx
        
        wx_version = V(wx.__version__).version
        
        if wx_version < [2, 8]:
            raise ValueError("requires wxPython >= 2.8, but you have %s" % wx.__version__)
        
        from IPython.lib.inputhookwx import inputhook_wx
        self.manager.set_inputhook(inputhook_wx)
        if _use_appnope():
            from appnope import nope
            nope()

        import wx
        if app is None:
            app = wx.GetApp()
        if app is None:
            app = wx.App(redirect=False, clearSigInt=False)

        return app

    def disable(self):
        """DEPRECATED since IPython 5.0

        Disable event loop integration with wxPython.

        This restores appnapp on OS X
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        if _use_appnope():
            from appnope import nap
            nap()

@inputhook_manager.register('qt', 'qt4')
class Qt4InputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with PyQt4.
        
        Parameters
        ----------
        app : Qt Application, optional.
            Running application to use.  If not given, we probe Qt for an
            existing application object, and create a new one if none is found.

        Notes
        -----
        This methods sets the PyOS_InputHook for PyQt4, which allows
        the PyQt4 to integrate with terminal based applications like
        IPython.

        If ``app`` is not given we probe for an existing one, and return it if
        found.  If no existing app is found, we create an :class:`QApplication`
        as follows::

            from PyQt4 import QtCore
            app = QtGui.QApplication(sys.argv)
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        from IPython.lib.inputhookqt4 import create_inputhook_qt4
        app, inputhook_qt4 = create_inputhook_qt4(self.manager, app)
        self.manager.set_inputhook(inputhook_qt4)
        if _use_appnope():
            from appnope import nope
            nope()

        return app

    def disable_qt4(self):
        """DEPRECATED since IPython 5.0

        Disable event loop integration with PyQt4.

        This restores appnapp on OS X
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        if _use_appnope():
            from appnope import nap
            nap()


@inputhook_manager.register('qt5')
class Qt5InputHook(Qt4InputHook):
    def enable(self, app=None):
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        os.environ['QT_API'] = 'pyqt5'
        return Qt4InputHook.enable(self, app)


@inputhook_manager.register('gtk')
class GtkInputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with PyGTK.

        Parameters
        ----------
        app : ignored
           Ignored, it's only a placeholder to keep the call signature of all
           gui activation methods consistent, which simplifies the logic of
           supporting magics.

        Notes
        -----
        This methods sets the PyOS_InputHook for PyGTK, which allows
        the PyGTK to integrate with terminal based applications like
        IPython.
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        import gtk
        try:
            gtk.set_interactive(True)
        except AttributeError:
            # For older versions of gtk, use our own ctypes version
            from IPython.lib.inputhookgtk import inputhook_gtk
            self.manager.set_inputhook(inputhook_gtk)


@inputhook_manager.register('tk')
class TkInputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with Tk.

        Parameters
        ----------
        app : toplevel :class:`Tkinter.Tk` widget, optional.
            Running toplevel widget to use.  If not given, we probe Tk for an
            existing one, and create a new one if none is found.

        Notes
        -----
        If you have already created a :class:`Tkinter.Tk` object, the only
        thing done by this method is to register with the
        :class:`InputHookManager`, since creating that object automatically
        sets ``PyOS_InputHook``.
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        if app is None:
            try:
                from tkinter import Tk  # Py 3
            except ImportError:
                from Tkinter import Tk  # Py 2
            app = Tk()
            app.withdraw()
            self.manager.apps[GUI_TK] = app
            return app


@inputhook_manager.register('glut')
class GlutInputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with GLUT.

        Parameters
        ----------

        app : ignored
            Ignored, it's only a placeholder to keep the call signature of all
            gui activation methods consistent, which simplifies the logic of
            supporting magics.

        Notes
        -----

        This methods sets the PyOS_InputHook for GLUT, which allows the GLUT to
        integrate with terminal based applications like IPython. Due to GLUT
        limitations, it is currently not possible to start the event loop
        without first creating a window. You should thus not create another
        window but use instead the created one. See 'gui-glut.py' in the
        docs/examples/lib directory.
        
        The default screen mode is set to:
        glut.GLUT_DOUBLE | glut.GLUT_RGBA | glut.GLUT_DEPTH
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)

        import OpenGL.GLUT as glut
        from IPython.lib.inputhookglut import glut_display_mode, \
                                              glut_close, glut_display, \
                                              glut_idle, inputhook_glut

        if GUI_GLUT not in self.manager.apps:
            glut.glutInit( sys.argv )
            glut.glutInitDisplayMode( glut_display_mode )
            # This is specific to freeglut
            if bool(glut.glutSetOption):
                glut.glutSetOption( glut.GLUT_ACTION_ON_WINDOW_CLOSE,
                                    glut.GLUT_ACTION_GLUTMAINLOOP_RETURNS )
            glut.glutCreateWindow( sys.argv[0] )
            glut.glutReshapeWindow( 1, 1 )
            glut.glutHideWindow( )
            glut.glutWMCloseFunc( glut_close )
            glut.glutDisplayFunc( glut_display )
            glut.glutIdleFunc( glut_idle )
        else:
            glut.glutWMCloseFunc( glut_close )
            glut.glutDisplayFunc( glut_display )
            glut.glutIdleFunc( glut_idle)
        self.manager.set_inputhook( inputhook_glut )


    def disable(self):
        """DEPRECATED since IPython 5.0

        Disable event loop integration with glut.
        
        This sets PyOS_InputHook to NULL and set the display function to a
        dummy one and set the timer to a dummy timer that will be triggered
        very far in the future.
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        import OpenGL.GLUT as glut
        from glut_support import glutMainLoopEvent

        glut.glutHideWindow() # This is an event to be processed below
        glutMainLoopEvent()
        super(GlutInputHook, self).disable()

@inputhook_manager.register('pyglet')
class PygletInputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with pyglet.

        Parameters
        ----------
        app : ignored
           Ignored, it's only a placeholder to keep the call signature of all
           gui activation methods consistent, which simplifies the logic of
           supporting magics.

        Notes
        -----
        This methods sets the ``PyOS_InputHook`` for pyglet, which allows
        pyglet to integrate with terminal based applications like
        IPython.

        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        from IPython.lib.inputhookpyglet import inputhook_pyglet
        self.manager.set_inputhook(inputhook_pyglet)
        return app


@inputhook_manager.register('gtk3')
class Gtk3InputHook(InputHookBase):
    def enable(self, app=None):
        """DEPRECATED since IPython 5.0

        Enable event loop integration with Gtk3 (gir bindings).

        Parameters
        ----------
        app : ignored
           Ignored, it's only a placeholder to keep the call signature of all
           gui activation methods consistent, which simplifies the logic of
           supporting magics.

        Notes
        -----
        This methods sets the PyOS_InputHook for Gtk3, which allows
        the Gtk3 to integrate with terminal based applications like
        IPython.
        """
        warn("This function is deprecated since IPython 5.0 and will be removed in future versions.",
                DeprecationWarning, stacklevel=2)
        from IPython.lib.inputhookgtk3 import inputhook_gtk3
        self.manager.set_inputhook(inputhook_gtk3)


clear_inputhook = inputhook_manager.clear_inputhook
set_inputhook = inputhook_manager.set_inputhook
current_gui = inputhook_manager.current_gui
clear_app_refs = inputhook_manager.clear_app_refs
enable_gui = inputhook_manager.enable_gui
disable_gui = inputhook_manager.disable_gui
register = inputhook_manager.register
guis = inputhook_manager.guihooks


def _deprecated_disable():
    warn("This function is deprecated since IPython 4.0 use disable_gui() instead",
            DeprecationWarning, stacklevel=2)
    inputhook_manager.disable_gui()
    
disable_wx = disable_qt4 = disable_gtk = disable_gtk3 = disable_glut = \
        disable_pyglet = disable_osx = _deprecated_disable
