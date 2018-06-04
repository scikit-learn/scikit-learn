from __future__ import absolute_import, division, print_function

import six

"""
Core functions and attributes for the matplotlib style library:

``use``
    Select style sheet to override the current matplotlib settings.
``context``
    Context manager to use a style sheet temporarily.
``available``
    List available style sheets.
``library``
    A dictionary of style names and matplotlib settings.
"""
import os
import re
import contextlib
import warnings

import matplotlib as mpl
from matplotlib import rc_params_from_file, rcParamsDefault


__all__ = ['use', 'context', 'available', 'library', 'reload_library']


BASE_LIBRARY_PATH = os.path.join(mpl.get_data_path(), 'stylelib')
# Users may want multiple library paths, so store a list of paths.
USER_LIBRARY_PATHS = [os.path.join(mpl._get_configdir(), 'stylelib')]
STYLE_EXTENSION = 'mplstyle'
STYLE_FILE_PATTERN = re.compile(r'([\S]+).%s$' % STYLE_EXTENSION)


# A list of rcParams that should not be applied from styles
STYLE_BLACKLIST = {
    'interactive', 'backend', 'backend.qt4', 'webagg.port', 'webagg.address',
    'webagg.port_retries', 'webagg.open_in_browser', 'backend_fallback',
    'toolbar', 'timezone', 'datapath', 'figure.max_open_warning',
    'savefig.directory', 'tk.window_focus', 'docstring.hardcopy'}


def _remove_blacklisted_style_params(d, warn=True):
    o = {}
    for key, val in d.items():
        if key in STYLE_BLACKLIST:
            if warn:
                warnings.warn(
                    "Style includes a parameter, '{0}', that is not related "
                    "to style.  Ignoring".format(key))
        else:
            o[key] = val
    return o


def is_style_file(filename):
    """Return True if the filename looks like a style file."""
    return STYLE_FILE_PATTERN.match(filename) is not None


def _apply_style(d, warn=True):
    mpl.rcParams.update(_remove_blacklisted_style_params(d, warn=warn))


def use(style):
    """Use matplotlib style settings from a style specification.

    The style name of 'default' is reserved for reverting back to
    the default style settings.

    Parameters
    ----------
    style : str, dict, or list
        A style specification. Valid options are:

        +------+-------------------------------------------------------------+
        | str  | The name of a style or a path/URL to a style file. For a    |
        |      | list of available style names, see `style.available`.       |
        +------+-------------------------------------------------------------+
        | dict | Dictionary with valid key/value pairs for                   |
        |      | `matplotlib.rcParams`.                                      |
        +------+-------------------------------------------------------------+
        | list | A list of style specifiers (str or dict) applied from first |
        |      | to last in the list.                                        |
        +------+-------------------------------------------------------------+


    """
    style_alias = {'mpl20': 'default',
                   'mpl15': 'classic'}
    if isinstance(style, six.string_types) or hasattr(style, 'keys'):
        # If name is a single str or dict, make it a single element list.
        styles = [style]
    else:
        styles = style

    styles = (style_alias.get(s, s)
              if isinstance(s, six.string_types)
              else s
              for s in styles)
    for style in styles:
        if not isinstance(style, six.string_types):
            _apply_style(style)
        elif style == 'default':
            _apply_style(rcParamsDefault, warn=False)
        elif style in library:
            _apply_style(library[style])
        else:
            try:
                rc = rc_params_from_file(style, use_default_template=False)
                _apply_style(rc)
            except IOError:
                raise IOError(
                    "{!r} not found in the style library and input is not a "
                    "valid URL or path; see `style.available` for list of "
                    "available styles".format(style))


@contextlib.contextmanager
def context(style, after_reset=False):
    """Context manager for using style settings temporarily.

    Parameters
    ----------
    style : str, dict, or list
        A style specification. Valid options are:

        +------+-------------------------------------------------------------+
        | str  | The name of a style or a path/URL to a style file. For a    |
        |      | list of available style names, see `style.available`.       |
        +------+-------------------------------------------------------------+
        | dict | Dictionary with valid key/value pairs for                   |
        |      | `matplotlib.rcParams`.                                      |
        +------+-------------------------------------------------------------+
        | list | A list of style specifiers (str or dict) applied from first |
        |      | to last in the list.                                        |
        +------+-------------------------------------------------------------+

    after_reset : bool
        If True, apply style after resetting settings to their defaults;
        otherwise, apply style on top of the current settings.
    """
    initial_settings = mpl.rcParams.copy()
    if after_reset:
        mpl.rcdefaults()
    try:
        use(style)
    except:
        # Restore original settings before raising errors during the update.
        mpl.rcParams.update(initial_settings)
        raise
    else:
        yield
    finally:
        mpl.rcParams.update(initial_settings)


def load_base_library():
    """Load style library defined in this package."""
    library = dict()
    library.update(read_style_directory(BASE_LIBRARY_PATH))
    return library


def iter_user_libraries():
    for stylelib_path in USER_LIBRARY_PATHS:
        stylelib_path = os.path.expanduser(stylelib_path)
        if os.path.exists(stylelib_path) and os.path.isdir(stylelib_path):
            yield stylelib_path


def update_user_library(library):
    """Update style library with user-defined rc files"""
    for stylelib_path in iter_user_libraries():
        styles = read_style_directory(stylelib_path)
        update_nested_dict(library, styles)
    return library


def iter_style_files(style_dir):
    """Yield file path and name of styles in the given directory."""
    for path in os.listdir(style_dir):
        filename = os.path.basename(path)
        if is_style_file(filename):
            match = STYLE_FILE_PATTERN.match(filename)
            path = os.path.abspath(os.path.join(style_dir, path))
            yield path, match.groups()[0]


def read_style_directory(style_dir):
    """Return dictionary of styles defined in `style_dir`."""
    styles = dict()
    for path, name in iter_style_files(style_dir):
        with warnings.catch_warnings(record=True) as warns:
            styles[name] = rc_params_from_file(path,
                                               use_default_template=False)

        for w in warns:
            message = 'In %s: %s' % (path, w.message)
            warnings.warn(message)

    return styles


def update_nested_dict(main_dict, new_dict):
    """Update nested dict (only level of nesting) with new values.

    Unlike dict.update, this assumes that the values of the parent dict are
    dicts (or dict-like), so you shouldn't replace the nested dict if it
    already exists. Instead you should update the sub-dict.
    """
    # update named styles specified by user
    for name, rc_dict in six.iteritems(new_dict):
        if name in main_dict:
            main_dict[name].update(rc_dict)
        else:
            main_dict[name] = rc_dict
    return main_dict


# Load style library
# ==================
_base_library = load_base_library()

library = None
available = []


def reload_library():
    """Reload style library."""
    global library
    available[:] = library = update_user_library(_base_library)
reload_library()
