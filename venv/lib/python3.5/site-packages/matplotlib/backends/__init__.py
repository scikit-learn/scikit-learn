from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib
import inspect
import traceback
import warnings
import logging

_log = logging.getLogger(__name__)

backend = matplotlib.get_backend()
_backend_loading_tb = "".join(
    line for line in traceback.format_stack()
    # Filter out line noise from importlib line.
    if not line.startswith('  File "<frozen importlib._bootstrap'))


def pylab_setup(name=None):
    '''return new_figure_manager, draw_if_interactive and show for pyplot

    This provides the backend-specific functions that are used by
    pyplot to abstract away the difference between interactive backends.

    Parameters
    ----------
    name : str, optional
        The name of the backend to use.  If `None`, falls back to
        ``matplotlib.get_backend()`` (which return :rc:`backend`).

    Returns
    -------
    backend_mod : module
        The module which contains the backend of choice

    new_figure_manager : function
        Create a new figure manager (roughly maps to GUI window)

    draw_if_interactive : function
        Redraw the current figure if pyplot is interactive

    show : function
        Show (and possibly block) any unshown figures.

    '''
    # Import the requested backend into a generic module object
    if name is None:
        # validates, to match all_backends
        name = matplotlib.get_backend()
    if name.startswith('module://'):
        backend_name = name[9:]
    else:
        backend_name = 'backend_' + name
        backend_name = backend_name.lower()  # until we banish mixed case
        backend_name = 'matplotlib.backends.%s' % backend_name.lower()

    # the last argument is specifies whether to use absolute or relative
    # imports. 0 means only perform absolute imports.
    backend_mod = __import__(backend_name, globals(), locals(),
                             [backend_name], 0)

    # Things we pull in from all backends
    new_figure_manager = backend_mod.new_figure_manager

    # image backends like pdf, agg or svg do not need to do anything
    # for "show" or "draw_if_interactive", so if they are not defined
    # by the backend, just do nothing
    def do_nothing_show(*args, **kwargs):
        frame = inspect.currentframe()
        fname = frame.f_back.f_code.co_filename
        if fname in ('<stdin>', '<ipython console>'):
            warnings.warn("""
Your currently selected backend, '%s' does not support show().
Please select a GUI backend in your matplotlibrc file ('%s')
or with matplotlib.use()""" %
                          (name, matplotlib.matplotlib_fname()))

    def do_nothing(*args, **kwargs):
        pass

    backend_version = getattr(backend_mod, 'backend_version', 'unknown')

    show = getattr(backend_mod, 'show', do_nothing_show)

    draw_if_interactive = getattr(backend_mod, 'draw_if_interactive',
                                  do_nothing)

    _log.debug('backend %s version %s', name, backend_version)

    # need to keep a global reference to the backend for compatibility
    # reasons. See https://github.com/matplotlib/matplotlib/issues/6092
    global backend
    backend = name
    return backend_mod, new_figure_manager, draw_if_interactive, show
