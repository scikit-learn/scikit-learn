"""
The rcsetup module contains the default values and the validation code for
customization using matplotlib's rc settings.

Each rc setting is assigned a default value and a function used to validate
any attempted changes to that setting. The default values and validation
functions are defined in the rcsetup module, and are used to construct the
rcParams global object which stores the settings and is referenced throughout
matplotlib.

These default values should be consistent with the default matplotlibrc file
that actually reflects the values given here. Any additions or deletions to the
parameter set listed here should also be visited to the
:file:`matplotlibrc.template` in matplotlib's root source directory.
"""
from __future__ import absolute_import, division, print_function

import six

from collections import Iterable, Mapping
from functools import reduce
import operator
import os
import warnings
import re

from matplotlib import cbook, testing
from matplotlib.cbook import mplDeprecation, deprecated, ls_mapper
from matplotlib.fontconfig_pattern import parse_fontconfig_pattern
from matplotlib.colors import is_color_like

# Don't let the original cycler collide with our validating cycler
from cycler import Cycler, cycler as ccycler


# The capitalized forms are needed for ipython at present; this may
# change for later versions.
interactive_bk = ['GTK', 'GTKAgg', 'GTKCairo', 'MacOSX',
                  'Qt4Agg', 'Qt5Agg', 'TkAgg', 'WX', 'WXAgg',
                  'GTK3Cairo', 'GTK3Agg', 'WebAgg', 'nbAgg']
interactive_bk = ['GTK', 'GTKAgg', 'GTKCairo', 'GTK3Agg', 'GTK3Cairo',
                  'MacOSX',
                  'nbAgg',
                  'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo',
                  'TkAgg', 'TkCairo',
                  'WebAgg',
                  'WX', 'WXAgg', 'WXCairo']
non_interactive_bk = ['agg', 'cairo', 'gdk',
                      'pdf', 'pgf', 'ps', 'svg', 'template']
all_backends = interactive_bk + non_interactive_bk


class ValidateInStrings(object):
    def __init__(self, key, valid, ignorecase=False):
        'valid is a list of legal strings'
        self.key = key
        self.ignorecase = ignorecase

        def func(s):
            if ignorecase:
                return s.lower()
            else:
                return s
        self.valid = {func(k): k for k in valid}

    def __call__(self, s):
        if self.ignorecase:
            s = s.lower()
        if s in self.valid:
            return self.valid[s]
        raise ValueError('Unrecognized %s string "%s": valid strings are %s'
                         % (self.key, s, list(six.itervalues(self.valid))))


def _listify_validator(scalar_validator, allow_stringlist=False):
    def f(s):
        if isinstance(s, six.string_types):
            try:
                return [scalar_validator(v.strip()) for v in s.split(',')
                        if v.strip()]
            except Exception:
                if allow_stringlist:
                    # Sometimes, a list of colors might be a single string
                    # of single-letter colornames. So give that a shot.
                    return [scalar_validator(v.strip()) for v in s if v.strip()]
                else:
                    raise
        # We should allow any generic sequence type, including generators,
        # Numpy ndarrays, and pandas data structures.  However, unordered
        # sequences, such as sets, should be allowed but discouraged unless the
        # user desires pseudorandom behavior.
        elif isinstance(s, Iterable) and not isinstance(s, Mapping):
            # The condition on this list comprehension will preserve the
            # behavior of filtering out any empty strings (behavior was
            # from the original validate_stringlist()), while allowing
            # any non-string/text scalar values such as numbers and arrays.
            return [scalar_validator(v) for v in s
                    if not isinstance(v, six.string_types) or v]
        else:
            raise ValueError("{!r} must be of type: string or non-dictionary "
                             "iterable".format(s))
    try:
        f.__name__ = "{}list".format(scalar_validator.__name__)
    except AttributeError:  # class instance.
        f.__name__ = "{}List".format(type(scalar_validator).__name__)
    f.__doc__ = scalar_validator.__doc__
    return f


def validate_any(s):
    return s
validate_anylist = _listify_validator(validate_any)


def validate_path_exists(s):
    """If s is a path, return s, else False"""
    if s is None:
        return None
    if os.path.exists(s):
        return s
    else:
        raise RuntimeError('"%s" should be a path but it does not exist' % s)


def validate_bool(b):
    """Convert b to a boolean or raise"""
    if isinstance(b, six.string_types):
        b = b.lower()
    if b in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
        return True
    elif b in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
        return False
    else:
        raise ValueError('Could not convert "%s" to boolean' % b)


def validate_bool_maybe_none(b):
    'Convert b to a boolean or raise'
    if isinstance(b, six.string_types):
        b = b.lower()
    if b is None or b == 'none':
        return None
    if b in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
        return True
    elif b in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
        return False
    else:
        raise ValueError('Could not convert "%s" to boolean' % b)


def deprecate_axes_hold(value):
    if value is None:
        return None  # converted to True where accessed in figure.py,
                     # axes/_base.py
    warnings.warn("axes.hold is deprecated, will be removed in 3.0",
                  mplDeprecation)
    return validate_bool(value)


def validate_float(s):
    """convert s to float or raise"""
    try:
        return float(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to float' % s)
validate_floatlist = _listify_validator(validate_float)


def validate_float_or_None(s):
    """convert s to float, None or raise"""
    # values directly from the rc file can only be strings,
    # so we need to recognize the string "None" and convert
    # it into the object. We will be case-sensitive here to
    # avoid confusion between string values of 'none', which
    # can be a valid string value for some other parameters.
    if s is None or s == 'None':
        return None
    try:
        return float(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to float or None' % s)


def validate_string_or_None(s):
    """convert s to string or raise"""
    if s is None:
        return None
    try:
        return validate_string(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to string' % s)


def validate_axisbelow(s):
    try:
        return validate_bool(s)
    except ValueError:
        if isinstance(s, six.string_types):
            s = s.lower()
            if s.startswith('line'):
                return 'line'
    raise ValueError('%s cannot be interpreted as'
                     ' True, False, or "line"' % s)


def validate_dpi(s):
    """confirm s is string 'figure' or convert s to float or raise"""
    if s == 'figure':
        return s
    try:
        return float(s)
    except ValueError:
        raise ValueError('"%s" is not string "figure" or'
            ' could not convert "%s" to float' % (s, s))


def validate_int(s):
    """convert s to int or raise"""
    try:
        return int(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to int' % s)


def validate_int_or_None(s):
    """if not None, tries to validate as an int"""
    if s=='None':
        s = None
    if s is None:
        return None
    try:
        return int(s)
    except ValueError:
        raise ValueError('Could not convert "%s" to int' % s)


def validate_fonttype(s):
    """
    confirm that this is a Postscript of PDF font type that we know how to
    convert to
    """
    fonttypes = {'type3':    3,
                 'truetype': 42}
    try:
        fonttype = validate_int(s)
    except ValueError:
        try:
            return fonttypes[s.lower()]
        except KeyError:
            raise ValueError(
                'Supported Postscript/PDF font types are %s' % list(fonttypes))
    else:
        if fonttype not in six.itervalues(fonttypes):
            raise ValueError(
                'Supported Postscript/PDF font types are %s' %
                list(six.itervalues(fonttypes)))
        return fonttype


_validate_standard_backends = ValidateInStrings(
    'backend', all_backends, ignorecase=True)


def validate_backend(s):
    if s.startswith('module://'):
        return s
    else:
        return _validate_standard_backends(s)


def validate_qt4(s):
    if s is None:
        return None
    return ValidateInStrings("backend.qt4", ['PyQt4', 'PySide', 'PyQt4v2'])(s)


def validate_qt5(s):
    if s is None:
        return None
    return ValidateInStrings("backend.qt5", ['PyQt5', 'PySide2'])(s)


def validate_toolbar(s):
    validator = ValidateInStrings(
                'toolbar',
                ['None', 'toolbar2', 'toolmanager'],
                ignorecase=True)
    return validator(s)


_seq_err_msg = ('You must supply exactly {n} values, you provided {num} '
                'values: {s}')

_str_err_msg = ('You must supply exactly {n} comma-separated values, you '
                'provided {num} comma-separated values: {s}')


class validate_nseq_float(object):
    def __init__(self, n=None, allow_none=False):
        self.n = n
        self.allow_none = allow_none

    def __call__(self, s):
        """return a seq of n floats or raise"""
        if isinstance(s, six.string_types):
            s = [x.strip() for x in s.split(',')]
            err_msg = _str_err_msg
        else:
            err_msg = _seq_err_msg

        if self.n is not None and len(s) != self.n:
            raise ValueError(err_msg.format(n=self.n, num=len(s), s=s))

        try:
            return [float(val)
                    if not self.allow_none or val is not None
                    else val
                    for val in s]
        except ValueError:
            raise ValueError('Could not convert all entries to floats')


class validate_nseq_int(object):
    def __init__(self, n=None):
        self.n = n

    def __call__(self, s):
        """return a seq of n ints or raise"""
        if isinstance(s, six.string_types):
            s = [x.strip() for x in s.split(',')]
            err_msg = _str_err_msg
        else:
            err_msg = _seq_err_msg

        if self.n is not None and len(s) != self.n:
            raise ValueError(err_msg.format(n=self.n, num=len(s), s=s))

        try:
            return [int(val) for val in s]
        except ValueError:
            raise ValueError('Could not convert all entries to ints')


def validate_color_or_inherit(s):
    'return a valid color arg'
    if s == 'inherit':
        return s
    return validate_color(s)


def validate_color_or_auto(s):
    if s == 'auto':
        return s
    return validate_color(s)


def validate_color_for_prop_cycle(s):
    # Special-case the N-th color cycle syntax, this obviously can not
    # go in the color cycle.
    if isinstance(s, bytes):
        match = re.match(b'^C[0-9]$', s)
        if match is not None:
            raise ValueError('Can not put cycle reference ({cn!r}) in '
                             'prop_cycler'.format(cn=s))
    elif isinstance(s, six.string_types):
        match = re.match('^C[0-9]$', s)
        if match is not None:
            raise ValueError('Can not put cycle reference ({cn!r}) in '
                             'prop_cycler'.format(cn=s))
    return validate_color(s)


def validate_color(s):
    'return a valid color arg'
    try:
        if s.lower() == 'none':
            return 'none'
    except AttributeError:
        pass

    if isinstance(s, six.string_types):
        if len(s) == 6 or len(s) == 8:
            stmp = '#' + s
            if is_color_like(stmp):
                return stmp

    if is_color_like(s):
        return s

    # If it is still valid, it must be a tuple.
    colorarg = s
    msg = ''
    if s.find(',') >= 0:
        # get rid of grouping symbols
        stmp = ''.join([c for c in s if c.isdigit() or c == '.' or c == ','])
        vals = stmp.split(',')
        if len(vals) not in [3, 4]:
            msg = '\nColor tuples must be of length 3 or 4'
        else:
            try:
                colorarg = [float(val) for val in vals]
            except ValueError:
                msg = '\nCould not convert all entries to floats'

    if not msg and is_color_like(colorarg):
        return colorarg

    raise ValueError('%s does not look like a color arg%s' % (s, msg))


validate_colorlist = _listify_validator(validate_color, allow_stringlist=True)
validate_colorlist.__doc__ = 'return a list of colorspecs'

def validate_string(s):
    if isinstance(s, (str, six.text_type)):
        # Always leave str as str and unicode as unicode
        return s
    else:
        return str(s)

validate_stringlist = _listify_validator(str)
validate_stringlist.__doc__ = 'return a list'

validate_orientation = ValidateInStrings(
    'orientation', ['landscape', 'portrait'])


def validate_aspect(s):
    if s in ('auto', 'equal'):
        return s
    try:
        return float(s)
    except ValueError:
        raise ValueError('not a valid aspect specification')


def validate_fontsize(s):
    fontsizes = ['xx-small', 'x-small', 'small', 'medium', 'large',
                 'x-large', 'xx-large', 'smaller', 'larger']
    if isinstance(s, six.string_types):
        s = s.lower()
    if s in fontsizes:
        return s
    try:
        return float(s)
    except ValueError:
        raise ValueError("%s is not a valid font size. Valid font sizes "
                         "are %s." % (s, ", ".join(fontsizes)))


validate_fontsizelist = _listify_validator(validate_fontsize)


def validate_font_properties(s):
    parse_fontconfig_pattern(s)
    return s


validate_fontset = ValidateInStrings(
    'fontset',
    ['dejavusans', 'dejavuserif', 'cm', 'stix', 'stixsans', 'custom'])

validate_mathtext_default = ValidateInStrings(
    'default',
    "rm cal it tt sf bf default bb frak circled scr regular".split())

validate_verbose = ValidateInStrings(
    'verbose',
    ['silent', 'helpful', 'debug', 'debug-annoying'])

_validate_alignment = ValidateInStrings(
    'alignment',
    ['center', 'top', 'bottom', 'baseline',
     'center_baseline'])

def validate_whiskers(s):
    if s == 'range':
        return 'range'
    else:
        try:
            v = validate_nseq_float(2)(s)
            return v
        except (TypeError, ValueError):
            try:
                v = float(s)
                return v
            except ValueError:
                raise ValueError("Not a valid whisker value ['range', float, "
                                 "(float, float)]")


def update_savefig_format(value):
    # The old savefig.extension could also have a value of "auto", but
    # the new savefig.format does not.  We need to fix this here.
    value = validate_string(value)
    if value == 'auto':
        value = 'png'
    return value


validate_ps_papersize = ValidateInStrings(
    'ps_papersize',
    ['auto', 'letter', 'legal', 'ledger',
    'a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10',
    'b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10',
    ], ignorecase=True)


def validate_ps_distiller(s):
    if isinstance(s, six.string_types):
        s = s.lower()
    if s in ('none', None):
        return None
    elif s in ('false', False):
        return False
    elif s in ('ghostscript', 'xpdf'):
        return s
    else:
        raise ValueError('matplotlibrc ps.usedistiller must either be none, '
                         'ghostscript or xpdf')

validate_joinstyle = ValidateInStrings('joinstyle',
                                       ['miter', 'round', 'bevel'],
                                       ignorecase=True)
validate_joinstylelist = _listify_validator(validate_joinstyle)

validate_capstyle = ValidateInStrings('capstyle',
                                      ['butt', 'round', 'projecting'],
                                      ignorecase=True)
validate_capstylelist = _listify_validator(validate_capstyle)

validate_fillstyle = ValidateInStrings('markers.fillstyle',
                                       ['full', 'left', 'right', 'bottom',
                                        'top', 'none'])
validate_fillstylelist = _listify_validator(validate_fillstyle)

_validate_negative_linestyle = ValidateInStrings('negative_linestyle',
                                                 ['solid', 'dashed'],
                                                 ignorecase=True)


@deprecated('2.1',
            addendum=(" See 'validate_negative_linestyle_legacy' " +
                      "deprecation warning for more information."))
def validate_negative_linestyle(s):
    return _validate_negative_linestyle(s)


@deprecated('2.1',
            addendum=(" The 'contour.negative_linestyle' rcParam now " +
                      "follows the same validation as the other rcParams " +
                      "that are related to line style."))
def validate_negative_linestyle_legacy(s):
    try:
        res = validate_negative_linestyle(s)
        return res
    except ValueError:
        dashes = validate_nseq_float(2)(s)
        return (0, dashes)  # (offset, (solid, blank))


validate_legend_loc = ValidateInStrings(
    'legend_loc',
    ['best',
     'upper right',
     'upper left',
     'lower left',
     'lower right',
     'right',
     'center left',
     'center right',
     'lower center',
     'upper center',
     'center'], ignorecase=True)


def validate_svg_fonttype(s):
    if s in ["none", "path"]:
        return s
    if s == "svgfont":
        cbook.warn_deprecated(
            "2.2", "'svgfont' support for svg.fonttype is deprecated.")
        return s
    raise ValueError("Unrecognized svg.fonttype string '{}'; "
                     "valid strings are 'none', 'path'")


def validate_hinting(s):
    if s in (True, False):
        return s
    if s.lower() in ('auto', 'native', 'either', 'none'):
        return s.lower()
    raise ValueError("hinting should be 'auto', 'native', 'either' or 'none'")

validate_pgf_texsystem = ValidateInStrings('pgf.texsystem',
                                           ['xelatex', 'lualatex', 'pdflatex'])

validate_movie_writer = ValidateInStrings('animation.writer',
    ['ffmpeg', 'ffmpeg_file',
     'avconv', 'avconv_file',
     'imagemagick', 'imagemagick_file',
     'html'])

validate_movie_frame_fmt = ValidateInStrings('animation.frame_format',
    ['png', 'jpeg', 'tiff', 'raw', 'rgba'])

validate_axis_locator = ValidateInStrings('major', ['minor', 'both', 'major'])

validate_movie_html_fmt = ValidateInStrings('animation.html',
    ['html5', 'jshtml', 'none'])

def validate_bbox(s):
    if isinstance(s, six.string_types):
        s = s.lower()
        if s == 'tight':
            return s
        if s == 'standard':
            return None
        raise ValueError("bbox should be 'tight' or 'standard'")
    elif s is not None:
        # Backwards compatibility. None is equivalent to 'standard'.
        raise ValueError("bbox should be 'tight' or 'standard'")
    return s

def validate_sketch(s):
    if isinstance(s, six.string_types):
        s = s.lower()
    if s == 'none' or s is None:
        return None
    if isinstance(s, six.string_types):
        result = tuple([float(v.strip()) for v in s.split(',')])
    elif isinstance(s, (list, tuple)):
        result = tuple([float(v) for v in s])
    if len(result) != 3:
        raise ValueError("path.sketch must be a tuple (scale, length, randomness)")
    return result

class ValidateInterval(object):
    """
    Value must be in interval
    """
    def __init__(self, vmin, vmax, closedmin=True, closedmax=True):
        self.vmin = vmin
        self.vmax = vmax
        self.cmin = closedmin
        self.cmax = closedmax

    def __call__(self, s):
        try:
            s = float(s)
        except ValueError:
            raise RuntimeError('Value must be a float; found "%s"' % s)

        if self.cmin and s < self.vmin:
            raise RuntimeError('Value must be >= %f; found "%f"' %
                               (self.vmin, s))
        elif not self.cmin and s <= self.vmin:
            raise RuntimeError('Value must be > %f; found "%f"' %
                               (self.vmin, s))

        if self.cmax and s > self.vmax:
            raise RuntimeError('Value must be <= %f; found "%f"' %
                               (self.vmax, s))
        elif not self.cmax and s >= self.vmax:
            raise RuntimeError('Value must be < %f; found "%f"' %
                               (self.vmax, s))
        return s

validate_grid_axis = ValidateInStrings('axes.grid.axis', ['x', 'y', 'both'])


def validate_hatch(s):
    """
    Validate a hatch pattern.
    A hatch pattern string can have any sequence of the following
    characters: ``\\ / | - + * . x o O``.

    """
    if not isinstance(s, six.string_types):
        raise ValueError("Hatch pattern must be a string")
    unknown = set(s) - {'\\', '/', '|', '-', '+', '*', '.', 'x', 'o', 'O'}
    if unknown:
        raise ValueError("Unknown hatch symbol(s): %s" % list(unknown))
    return s
validate_hatchlist = _listify_validator(validate_hatch)
validate_dashlist = _listify_validator(validate_nseq_float(allow_none=True))

_prop_validators = {
        'color': _listify_validator(validate_color_for_prop_cycle,
                                    allow_stringlist=True),
        'linewidth': validate_floatlist,
        'linestyle': validate_stringlist,
        'facecolor': validate_colorlist,
        'edgecolor': validate_colorlist,
        'joinstyle': validate_joinstylelist,
        'capstyle': validate_capstylelist,
        'fillstyle': validate_fillstylelist,
        'markerfacecolor': validate_colorlist,
        'markersize': validate_floatlist,
        'markeredgewidth': validate_floatlist,
        'markeredgecolor': validate_colorlist,
        'alpha': validate_floatlist,
        'marker': validate_stringlist,
        'hatch': validate_hatchlist,
        'dashes': validate_dashlist,
    }
_prop_aliases = {
        'c': 'color',
        'lw': 'linewidth',
        'ls': 'linestyle',
        'fc': 'facecolor',
        'ec': 'edgecolor',
        'mfc': 'markerfacecolor',
        'mec': 'markeredgecolor',
        'mew': 'markeredgewidth',
        'ms': 'markersize',
    }


def cycler(*args, **kwargs):
    """
    Creates a :class:`cycler.Cycler` object much like :func:`cycler.cycler`,
    but includes input validation.

    cycler(arg)
    cycler(label, itr)
    cycler(label1=itr1[, label2=itr2[, ...]])

    Form 1 simply copies a given `Cycler` object.

    Form 2 creates a `Cycler` from a label and an iterable.

    Form 3 composes a `Cycler` as an inner product of the
    pairs of keyword arguments. In other words, all of the
    iterables are cycled simultaneously, as if through zip().

    Parameters
    ----------
    arg : Cycler
        Copy constructor for Cycler.

    label : name
        The property key. Must be a valid `Artist` property.
        For example, 'color' or 'linestyle'. Aliases are allowed,
        such as 'c' for 'color' and 'lw' for 'linewidth'.

    itr : iterable
        Finite-length iterable of the property values. These values
        are validated and will raise a ValueError if invalid.

    Returns
    -------
    cycler : Cycler
        New :class:`cycler.Cycler` for the given properties

    """
    if args and kwargs:
        raise TypeError("cycler() can only accept positional OR keyword "
                        "arguments -- not both.")
    elif not args and not kwargs:
        raise TypeError("cycler() must have positional OR keyword arguments")

    if len(args) == 1:
        if not isinstance(args[0], Cycler):
            raise TypeError("If only one positional argument given, it must "
                            " be a Cycler instance.")
        return validate_cycler(args[0])
    elif len(args) == 2:
        pairs = [(args[0], args[1])]
    elif len(args) > 2:
        raise TypeError("No more than 2 positional arguments allowed")
    else:
        pairs = six.iteritems(kwargs)

    validated = []
    for prop, vals in pairs:
        norm_prop = _prop_aliases.get(prop, prop)
        validator = _prop_validators.get(norm_prop, None)
        if validator is None:
            raise TypeError("Unknown artist property: %s" % prop)
        vals = validator(vals)
        # We will normalize the property names as well to reduce
        # the amount of alias handling code elsewhere.
        validated.append((norm_prop, vals))

    return reduce(operator.add, (ccycler(k, v) for k, v in validated))


def validate_cycler(s):
    'return a Cycler object from a string repr or the object itself'
    if isinstance(s, six.string_types):
        try:
            # TODO: We might want to rethink this...
            # While I think I have it quite locked down,
            # it is execution of arbitrary code without
            # sanitation.
            # Combine this with the possibility that rcparams
            # might come from the internet (future plans), this
            # could be downright dangerous.
            # I locked it down by only having the 'cycler()' function
            # available.
            # UPDATE: Partly plugging a security hole.
            # I really should have read this:
            # http://nedbatchelder.com/blog/201206/eval_really_is_dangerous.html
            # We should replace this eval with a combo of PyParsing and
            # ast.literal_eval()
            if '.__' in s.replace(' ', ''):
                raise ValueError("'%s' seems to have dunder methods. Raising"
                                 " an exception for your safety")
            s = eval(s, {'cycler': cycler, '__builtins__': {}})
        except BaseException as e:
            raise ValueError("'%s' is not a valid cycler construction: %s" %
                             (s, e))
    # Should make sure what comes from the above eval()
    # is a Cycler object.
    if isinstance(s, Cycler):
        cycler_inst = s
    else:
        raise ValueError("object was not a string or Cycler instance: %s" % s)

    unknowns = cycler_inst.keys - (set(_prop_validators) | set(_prop_aliases))
    if unknowns:
        raise ValueError("Unknown artist properties: %s" % unknowns)

    # Not a full validation, but it'll at least normalize property names
    # A fuller validation would require v0.10 of cycler.
    checker = set()
    for prop in cycler_inst.keys:
        norm_prop = _prop_aliases.get(prop, prop)
        if norm_prop != prop and norm_prop in cycler_inst.keys:
            raise ValueError("Cannot specify both '{0}' and alias '{1}'"
                             " in the same prop_cycle".format(norm_prop, prop))
        if norm_prop in checker:
            raise ValueError("Another property was already aliased to '{0}'."
                             " Collision normalizing '{1}'.".format(norm_prop,
                                                                    prop))
        checker.update([norm_prop])

    # This is just an extra-careful check, just in case there is some
    # edge-case I haven't thought of.
    assert len(checker) == len(cycler_inst.keys)

    # Now, it should be safe to mutate this cycler
    for prop in cycler_inst.keys:
        norm_prop = _prop_aliases.get(prop, prop)
        cycler_inst.change_key(prop, norm_prop)

    for key, vals in cycler_inst.by_key().items():
        _prop_validators[key](vals)

    return cycler_inst


def validate_hist_bins(s):
    if isinstance(s, six.string_types) and s == 'auto':
        return s
    try:
        return int(s)
    except (TypeError, ValueError):
        pass

    try:
        return validate_floatlist(s)
    except ValueError:
        pass

    raise ValueError("'hist.bins' must be 'auto', an int or " +
                     "a sequence of floats")

def validate_animation_writer_path(p):
    # Make sure it's a string and then figure out if the animations
    # are already loaded and reset the writers (which will validate
    # the path on next call)
    if not isinstance(p, six.string_types):
        raise ValueError("path must be a (unicode) string")
    from sys import modules
    # set dirty, so that the next call to the registry will re-evaluate
    # the state.
    # only set dirty if already loaded. If not loaded, the load will
    # trigger the checks.
    if "matplotlib.animation" in modules:
        modules["matplotlib.animation"].writers.set_dirty()
    return p

def validate_webagg_address(s):
    if s is not None:
        import socket
        try:
            socket.inet_aton(s)
        except socket.error as e:
            raise ValueError("'webagg.address' is not a valid IP address")
        return s
    raise ValueError("'webagg.address' is not a valid IP address")

# A validator dedicated to the named line styles, based on the items in
# ls_mapper, and a list of possible strings read from Line2D.set_linestyle
_validate_named_linestyle = ValidateInStrings('linestyle',
                                              list(six.iterkeys(ls_mapper)) +
                                              list(six.itervalues(ls_mapper)) +
                                              ['None', 'none', ' ', ''],
                                              ignorecase=True)


def _validate_linestyle(ls):
    """
    A validator for all possible line styles, the named ones *and*
    the on-off ink sequences.
    """
    # Look first for a valid named line style, like '--' or 'solid'
    if isinstance(ls, six.string_types):
        try:
            return _validate_named_linestyle(ls)
        except (UnicodeDecodeError, KeyError):
            # On Python 2, string-like *ls*, like for example
            # 'solid'.encode('utf-16'), may raise a unicode error.
            raise ValueError("the linestyle string {!r} is not a valid "
                             "string.".format(ls))

    if isinstance(ls, (bytes, bytearray)):
        # On Python 2, a string-like *ls* should already have lead to a
        # successful return or to raising an exception. On Python 3, we have
        # to manually raise an exception in the case of a byte-like *ls*.
        # Otherwise, if *ls* is of even-length, it will be passed to the
        # instance of validate_nseq_float, which will return an absurd on-off
        # ink sequence...
        raise ValueError("linestyle {!r} neither looks like an on-off ink "
                         "sequence nor a valid string.".format(ls))

    # Look for an on-off ink sequence (in points) *of even length*.
    # Offset is set to None.
    try:
        if len(ls) % 2 != 0:
            raise ValueError("the linestyle sequence {!r} is not of even "
                             "length.".format(ls))

        return (None, validate_nseq_float()(ls))

    except (ValueError, TypeError):
        # TypeError can be raised inside the instance of validate_nseq_float,
        # by wrong types passed to float(), like NoneType.
        raise ValueError("linestyle {!r} is not a valid on-off ink "
                         "sequence.".format(ls))


# a map from key -> value, converter
defaultParams = {
    'backend':           ['Agg', validate_backend],  # agg is certainly
                                                      # present
    'backend_fallback':  [True, validate_bool],  # agg is certainly present
    'backend.qt4':       [None, validate_qt4],
    'backend.qt5':       [None, validate_qt5],
    'webagg.port':       [8988, validate_int],
    'webagg.address':    ['127.0.0.1', validate_webagg_address],
    'webagg.open_in_browser': [True, validate_bool],
    'webagg.port_retries': [50, validate_int],
    'nbagg.transparent':       [True, validate_bool],
    'toolbar':           ['toolbar2', validate_toolbar],
    'datapath':          [None, validate_path_exists],  # handled by
                                                        # _get_data_path_cached
    'interactive':       [False, validate_bool],
    'timezone':          ['UTC', validate_string],

    # the verbosity setting
    'verbose.level': ['silent', validate_verbose],
    'verbose.fileo': ['sys.stdout', validate_string],

    # line props
    'lines.linewidth':       [1.5, validate_float],  # line width in points
    'lines.linestyle':       ['-', _validate_linestyle],  # solid line
    'lines.color':           ['C0', validate_color],  # first color in color cycle
    'lines.marker':          ['None', validate_string],  # marker name
    'lines.markeredgewidth': [1.0, validate_float],
    'lines.markersize':      [6, validate_float],    # markersize, in points
    'lines.antialiased':     [True, validate_bool],  # antialiased (no jaggies)
    'lines.dash_joinstyle':  ['round', validate_joinstyle],
    'lines.solid_joinstyle': ['round', validate_joinstyle],
    'lines.dash_capstyle':   ['butt', validate_capstyle],
    'lines.solid_capstyle':  ['projecting', validate_capstyle],
    'lines.dashed_pattern':  [[3.7, 1.6], validate_nseq_float(allow_none=True)],
    'lines.dashdot_pattern': [[6.4, 1.6, 1, 1.6],
                              validate_nseq_float(allow_none=True)],
    'lines.dotted_pattern':  [[1, 1.65], validate_nseq_float(allow_none=True)],
    'lines.scale_dashes':  [True, validate_bool],

    # marker props
    'markers.fillstyle': ['full', validate_fillstyle],

    ## patch props
    'patch.linewidth':   [1.0, validate_float],     # line width in points
    'patch.edgecolor':   ['k', validate_color],
    'patch.force_edgecolor' : [False, validate_bool],
    'patch.facecolor':   ['C0', validate_color],    # first color in cycle
    'patch.antialiased': [True, validate_bool],     # antialiased (no jaggies)

    ## hatch props
    'hatch.color': ['k', validate_color],
    'hatch.linewidth': [1.0, validate_float],

    ## Histogram properties
    'hist.bins': [10, validate_hist_bins],

    ## Boxplot properties
    'boxplot.notch': [False, validate_bool],
    'boxplot.vertical': [True, validate_bool],
    'boxplot.whiskers': [1.5, validate_whiskers],
    'boxplot.bootstrap': [None, validate_int_or_None],
    'boxplot.patchartist': [False, validate_bool],
    'boxplot.showmeans': [False, validate_bool],
    'boxplot.showcaps': [True, validate_bool],
    'boxplot.showbox': [True, validate_bool],
    'boxplot.showfliers': [True, validate_bool],
    'boxplot.meanline': [False, validate_bool],

    'boxplot.flierprops.color': ['k', validate_color],
    'boxplot.flierprops.marker': ['o', validate_string],
    'boxplot.flierprops.markerfacecolor': ['none', validate_color_or_auto],
    'boxplot.flierprops.markeredgecolor': ['k', validate_color],
    'boxplot.flierprops.markersize': [6, validate_float],
    'boxplot.flierprops.linestyle': ['none', _validate_linestyle],
    'boxplot.flierprops.linewidth': [1.0, validate_float],

    'boxplot.boxprops.color': ['k', validate_color],
    'boxplot.boxprops.linewidth': [1.0, validate_float],
    'boxplot.boxprops.linestyle': ['-', _validate_linestyle],

    'boxplot.whiskerprops.color': ['k', validate_color],
    'boxplot.whiskerprops.linewidth': [1.0, validate_float],
    'boxplot.whiskerprops.linestyle': ['-', _validate_linestyle],

    'boxplot.capprops.color': ['k', validate_color],
    'boxplot.capprops.linewidth': [1.0, validate_float],
    'boxplot.capprops.linestyle': ['-', _validate_linestyle],

    'boxplot.medianprops.color': ['C1', validate_color],
    'boxplot.medianprops.linewidth': [1.0, validate_float],
    'boxplot.medianprops.linestyle': ['-', _validate_linestyle],

    'boxplot.meanprops.color': ['C2', validate_color],
    'boxplot.meanprops.marker': ['^', validate_string],
    'boxplot.meanprops.markerfacecolor': ['C2', validate_color],
    'boxplot.meanprops.markeredgecolor': ['C2', validate_color],
    'boxplot.meanprops.markersize': [6, validate_float],
    'boxplot.meanprops.linestyle': ['--', _validate_linestyle],
    'boxplot.meanprops.linewidth': [1.0, validate_float],

    ## font props
    'font.family':     [['sans-serif'], validate_stringlist],  # used by text object
    'font.style':      ['normal', validate_string],
    'font.variant':    ['normal', validate_string],
    'font.stretch':    ['normal', validate_string],
    'font.weight':     ['normal', validate_string],
    'font.size':       [10, validate_float],      # Base font size in points
    'font.serif':      [['DejaVu Serif', 'Bitstream Vera Serif',
                         'Computer Modern Roman',
                         'New Century Schoolbook', 'Century Schoolbook L',
                         'Utopia', 'ITC Bookman', 'Bookman',
                         'Nimbus Roman No9 L', 'Times New Roman',
                         'Times', 'Palatino', 'Charter', 'serif'],
                        validate_stringlist],
    'font.sans-serif': [['DejaVu Sans', 'Bitstream Vera Sans',
                         'Computer Modern Sans Serif',
                         'Lucida Grande', 'Verdana', 'Geneva', 'Lucid',
                         'Arial', 'Helvetica', 'Avant Garde', 'sans-serif'],
                        validate_stringlist],
    'font.cursive':    [['Apple Chancery', 'Textile', 'Zapf Chancery',
                         'Sand', 'Script MT', 'Felipa', 'cursive'],
                        validate_stringlist],
    'font.fantasy':    [['Comic Sans MS', 'Chicago', 'Charcoal', 'Impact'
                         'Western', 'Humor Sans', 'xkcd', 'fantasy'],
                        validate_stringlist],
    'font.monospace':  [['DejaVu Sans Mono', 'Bitstream Vera Sans Mono',
                         'Computer Modern Typewriter',
                         'Andale Mono', 'Nimbus Mono L', 'Courier New',
                         'Courier', 'Fixed', 'Terminal', 'monospace'],
                        validate_stringlist],

    # text props
    'text.color':          ['k', validate_color],     # black
    'text.usetex':         [False, validate_bool],
    'text.latex.unicode':  [False, validate_bool],
    'text.latex.preamble': [[''], validate_stringlist],
    'text.latex.preview':  [False, validate_bool],
    'text.dvipnghack':     [None, validate_bool_maybe_none],
    'text.hinting':        ['auto', validate_hinting],
    'text.hinting_factor': [8, validate_int],
    'text.antialiased':    [True, validate_bool],

    'mathtext.cal':            ['cursive', validate_font_properties],
    'mathtext.rm':             ['sans', validate_font_properties],
    'mathtext.tt':             ['monospace', validate_font_properties],
    'mathtext.it':             ['sans:italic', validate_font_properties],
    'mathtext.bf':             ['sans:bold', validate_font_properties],
    'mathtext.sf':             ['sans', validate_font_properties],
    'mathtext.fontset':        ['dejavusans', validate_fontset],
    'mathtext.default':        ['it', validate_mathtext_default],
    'mathtext.fallback_to_cm': [True, validate_bool],

    'image.aspect':        ['equal', validate_aspect],  # equal, auto, a number
    'image.interpolation': ['nearest', validate_string],
    'image.cmap':          ['viridis', validate_string],        # one of gray, jet, etc
    'image.lut':           [256, validate_int],  # lookup table
    'image.origin':        ['upper', validate_string],  # lookup table
    'image.resample':      [True, validate_bool],
    # Specify whether vector graphics backends will combine all images on a
    # set of axes into a single composite image
    'image.composite_image': [True, validate_bool],

    # contour props
    'contour.negative_linestyle': ['dashed', _validate_linestyle],
    'contour.corner_mask':        [True, validate_bool],

    # errorbar props
    'errorbar.capsize':      [0, validate_float],

    # axes props
    'axes.axisbelow':        ['line', validate_axisbelow],
    'axes.hold':             [None, deprecate_axes_hold],
    'axes.facecolor':        ['w', validate_color],  # background color; white
    'axes.edgecolor':        ['k', validate_color],  # edge color; black
    'axes.linewidth':        [0.8, validate_float],  # edge linewidth

    'axes.spines.left':      [True, validate_bool],  # Set visibility of axes
    'axes.spines.right':     [True, validate_bool],  # 'spines', the lines
    'axes.spines.bottom':    [True, validate_bool],  # around the chart
    'axes.spines.top':       [True, validate_bool],  # denoting data boundary

    'axes.titlesize':        ['large', validate_fontsize],  # fontsize of the
                                                            # axes title
    'axes.titleweight':      ['normal', validate_string],  # font weight of axes title
    'axes.titlepad':         [6.0, validate_float],  # pad from axes top to title in points
    'axes.grid':             [False, validate_bool],   # display grid or not
    'axes.grid.which':       ['major', validate_axis_locator],  # set whether the gid are by
                                                                # default draw on 'major'
                                                                # 'minor' or 'both' kind of
                                                                # axis locator
    'axes.grid.axis':        ['both', validate_grid_axis], # grid type.
                                                      # Can be 'x', 'y', 'both'
    'axes.labelsize':        ['medium', validate_fontsize],  # fontsize of the
                                                             # x any y labels
    'axes.labelpad':         [4.0, validate_float], # space between label and axis
    'axes.labelweight':      ['normal', validate_string],  # fontsize of the x any y labels
    'axes.labelcolor':       ['k', validate_color],    # color of axis label
    'axes.formatter.limits': [[-7, 7], validate_nseq_int(2)],
                               # use scientific notation if log10
                               # of the axis range is smaller than the
                               # first or larger than the second
    'axes.formatter.use_locale': [False, validate_bool],
                               # Use the current locale to format ticks
    'axes.formatter.use_mathtext': [False, validate_bool],
    'axes.formatter.min_exponent': [0, validate_int],  # minimum exponent to format in scientific notation
    'axes.formatter.useoffset': [True, validate_bool],
    'axes.formatter.offset_threshold': [4, validate_int],
    'axes.unicode_minus': [True, validate_bool],
    # This entry can be either a cycler object or a
    # string repr of a cycler-object, which gets eval()'ed
    # to create the object.
    'axes.prop_cycle': [
        ccycler('color',
                ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                 '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                 '#bcbd22', '#17becf']),
        validate_cycler],
    # If 'data', axes limits are set close to the data.
    # If 'round_numbers' axes limits are set to the nearest round numbers.
    'axes.autolimit_mode': [
        'data',
        ValidateInStrings('autolimit_mode', ['data', 'round_numbers'])],
    'axes.xmargin': [0.05, ValidateInterval(0, 1,
                                            closedmin=True,
                                            closedmax=True)],  # margin added to xaxis
    'axes.ymargin': [0.05, ValidateInterval(0, 1,
                                            closedmin=True,
                                            closedmax=True)],# margin added to yaxis

    'polaraxes.grid': [True, validate_bool],  # display polar grid or
                                                     # not
    'axes3d.grid': [True, validate_bool],  # display 3d grid

    # scatter props
    'scatter.marker': ['o', validate_string],

    # TODO validate that these are valid datetime format strings
    'date.autoformatter.year': ['%Y', validate_string],
    'date.autoformatter.month': ['%Y-%m', validate_string],
    'date.autoformatter.day': ['%Y-%m-%d', validate_string],
    'date.autoformatter.hour': ['%m-%d %H', validate_string],
    'date.autoformatter.minute': ['%d %H:%M', validate_string],
    'date.autoformatter.second': ['%H:%M:%S', validate_string],
    'date.autoformatter.microsecond': ['%M:%S.%f', validate_string],

    #legend properties
    'legend.fancybox': [True, validate_bool],
    'legend.loc': ['best', validate_legend_loc],
    # the number of points in the legend line
    'legend.numpoints': [1, validate_int],
    # the number of points in the legend line for scatter
    'legend.scatterpoints': [1, validate_int],
    'legend.fontsize': ['medium', validate_fontsize],
     # the relative size of legend markers vs. original
    'legend.markerscale': [1.0, validate_float],
    'legend.shadow': [False, validate_bool],
     # whether or not to draw a frame around legend
    'legend.frameon': [True, validate_bool],
     # alpha value of the legend frame
    'legend.framealpha': [0.8, validate_float_or_None],

    ## the following dimensions are in fraction of the font size
    'legend.borderpad': [0.4, validate_float],  # units are fontsize
    # the vertical space between the legend entries
    'legend.labelspacing': [0.5, validate_float],
    # the length of the legend lines
    'legend.handlelength': [2., validate_float],
    # the length of the legend lines
    'legend.handleheight': [0.7, validate_float],
    # the space between the legend line and legend text
    'legend.handletextpad': [.8, validate_float],
    # the border between the axes and legend edge
    'legend.borderaxespad': [0.5, validate_float],
    # the border between the axes and legend edge
    'legend.columnspacing': [2., validate_float],
    'legend.facecolor': ['inherit', validate_color_or_inherit],
    'legend.edgecolor': ['0.8', validate_color_or_inherit],

    # tick properties
    'xtick.top':         [False, validate_bool],   # draw ticks on the top side
    'xtick.bottom':      [True, validate_bool],   # draw ticks on the bottom side
    'xtick.labeltop':    [False, validate_bool],  # draw label on the top
    'xtick.labelbottom': [True, validate_bool],  # draw label on the bottom
    'xtick.major.size':  [3.5, validate_float],    # major xtick size in points
    'xtick.minor.size':  [2, validate_float],    # minor xtick size in points
    'xtick.major.width': [0.8, validate_float],  # major xtick width in points
    'xtick.minor.width': [0.6, validate_float],  # minor xtick width in points
    'xtick.major.pad':   [3.5, validate_float],    # distance to label in points
    'xtick.minor.pad':   [3.4, validate_float],    # distance to label in points
    'xtick.color':       ['k', validate_color],  # color of the xtick labels
    'xtick.minor.visible':   [False, validate_bool],    # visibility of the x axis minor ticks
    'xtick.minor.top':   [True, validate_bool],  # draw x axis top minor ticks
    'xtick.minor.bottom':    [True, validate_bool],    # draw x axis bottom minor ticks
    'xtick.major.top':   [True, validate_bool],  # draw x axis top major ticks
    'xtick.major.bottom':    [True, validate_bool],    # draw x axis bottom major ticks

    # fontsize of the xtick labels
    'xtick.labelsize':   ['medium', validate_fontsize],
    'xtick.direction':   ['out', validate_string],            # direction of xticks
    'xtick.alignment': ["center", _validate_alignment],

    'ytick.left':        [True, validate_bool],  # draw ticks on the left side
    'ytick.right':       [False, validate_bool],  # draw ticks on the right side
    'ytick.labelleft':   [True, validate_bool],  # draw tick labels on the left side
    'ytick.labelright':  [False, validate_bool],  # draw tick labels on the right side
    'ytick.major.size':  [3.5, validate_float],     # major ytick size in points
    'ytick.minor.size':  [2, validate_float],     # minor ytick size in points
    'ytick.major.width': [0.8, validate_float],   # major ytick width in points
    'ytick.minor.width': [0.6, validate_float],   # minor ytick width in points
    'ytick.major.pad':   [3.5, validate_float],     # distance to label in points
    'ytick.minor.pad':   [3.4, validate_float],     # distance to label in points
    'ytick.color':       ['k', validate_color],   # color of the ytick labels
    'ytick.minor.visible':   [False, validate_bool],    # visibility of the y axis minor ticks
    'ytick.minor.left':   [True, validate_bool],  # draw y axis left minor ticks
    'ytick.minor.right':    [True, validate_bool],    # draw y axis right minor ticks
    'ytick.major.left':   [True, validate_bool],  # draw y axis left major ticks
    'ytick.major.right':    [True, validate_bool],    # draw y axis right major ticks

    # fontsize of the ytick labels
    'ytick.labelsize':   ['medium', validate_fontsize],
    'ytick.direction':   ['out', validate_string],            # direction of yticks
    'ytick.alignment': ["center_baseline", _validate_alignment],


    'grid.color':        ['#b0b0b0', validate_color],  # grid color
    'grid.linestyle':    ['-', _validate_linestyle],  # solid
    'grid.linewidth':    [0.8, validate_float],     # in points
    'grid.alpha':        [1.0, validate_float],


    ## figure props
    # figure title
    'figure.titlesize':   ['large', validate_fontsize],
    'figure.titleweight': ['normal', validate_string],

    # figure size in inches: width by height
    'figure.figsize':    [[6.4, 4.8], validate_nseq_float(2)],
    'figure.dpi':        [100, validate_float],  # DPI
    'figure.facecolor':  ['w', validate_color],  # facecolor; white
    'figure.edgecolor':  ['w', validate_color],  # edgecolor; white
    'figure.frameon':    [True, validate_bool],
    'figure.autolayout': [False, validate_bool],
    'figure.max_open_warning': [20, validate_int],

    'figure.subplot.left': [0.125, ValidateInterval(0, 1, closedmin=True,
                                                       closedmax=True)],
    'figure.subplot.right': [0.9, ValidateInterval(0, 1, closedmin=True,
                                                     closedmax=True)],
    'figure.subplot.bottom': [0.11, ValidateInterval(0, 1, closedmin=True,
                                                     closedmax=True)],
    'figure.subplot.top': [0.88, ValidateInterval(0, 1, closedmin=True,
                                                     closedmax=True)],
    'figure.subplot.wspace': [0.2, ValidateInterval(0, 1, closedmin=True,
                                                     closedmax=False)],
    'figure.subplot.hspace': [0.2, ValidateInterval(0, 1, closedmin=True,
                                                     closedmax=False)],

    # do constrained_layout.
    'figure.constrained_layout.use': [False, validate_bool],
    # wspace and hspace are fraction of adjacent subplots to use
    # for space.  Much smaller than above because we don't need
    # room for the text.
    'figure.constrained_layout.hspace': [0.02, ValidateInterval(
            0, 1, closedmin=True, closedmax=False)],
    'figure.constrained_layout.wspace': [0.02, ValidateInterval(
            0, 1, closedmin=True, closedmax=False)],
    # This is a buffer around the axes in inches.  This is 3pts.
    'figure.constrained_layout.h_pad': [0.04167, validate_float],
    'figure.constrained_layout.w_pad': [0.04167, validate_float],

    ## Saving figure's properties
    'savefig.dpi':         ['figure', validate_dpi],  # DPI
    'savefig.facecolor':   ['w', validate_color],  # facecolor; white
    'savefig.edgecolor':   ['w', validate_color],  # edgecolor; white
    'savefig.frameon':     [True, validate_bool],
    'savefig.orientation': ['portrait', validate_orientation],  # edgecolor;
                                                                 #white
    'savefig.jpeg_quality': [95, validate_int],
    # value checked by backend at runtime
    'savefig.format':     ['png', update_savefig_format],
    # options are 'tight', or 'standard'. 'standard' validates to None.
    'savefig.bbox':       ['standard', validate_bbox],
    'savefig.pad_inches': [0.1, validate_float],
    # default directory in savefig dialog box
    'savefig.directory': ['~', validate_string],
    'savefig.transparent': [False, validate_bool],

    # Maintain shell focus for TkAgg
    'tk.window_focus':  [False, validate_bool],

    # Set the papersize/type
    'ps.papersize':     ['letter', validate_ps_papersize],
    'ps.useafm':        [False, validate_bool],  # Set PYTHONINSPECT
    # use ghostscript or xpdf to distill ps output
    'ps.usedistiller':  [False, validate_ps_distiller],
    'ps.distiller.res': [6000, validate_int],     # dpi
    'ps.fonttype':      [3, validate_fonttype],  # 3 (Type3) or 42 (Truetype)
    # compression level from 0 to 9; 0 to disable
    'pdf.compression':  [6, validate_int],
    # ignore any color-setting commands from the frontend
    'pdf.inheritcolor': [False, validate_bool],
    # use only the 14 PDF core fonts embedded in every PDF viewing application
    'pdf.use14corefonts': [False, validate_bool],
    'pdf.fonttype':     [3, validate_fonttype],  # 3 (Type3) or 42 (Truetype)

    'pgf.debug':     [False, validate_bool],  # output debug information
    # choose latex application for creating pdf files (xelatex/lualatex)
    'pgf.texsystem': ['xelatex', validate_pgf_texsystem],
    # use matplotlib rc settings for font configuration
    'pgf.rcfonts':   [True, validate_bool],
    # provide a custom preamble for the latex process
    'pgf.preamble':  [[''], validate_stringlist],

    # write raster image data directly into the svg file
    'svg.image_inline':     [True, validate_bool],
    # True to save all characters as paths in the SVG
    'svg.fonttype':         ['path', validate_svg_fonttype],
    'svg.hashsalt':         [None, validate_string_or_None],

    # set this when you want to generate hardcopy docstring
    'docstring.hardcopy': [False, validate_bool],
    # where plugin directory is locate
    'plugins.directory':  ['.matplotlib_plugins', validate_string],

    'path.simplify': [True, validate_bool],
    'path.simplify_threshold': [1.0 / 9.0, ValidateInterval(0.0, 1.0)],
    'path.snap': [True, validate_bool],
    'path.sketch': [None, validate_sketch],
    'path.effects': [[], validate_any],
    'agg.path.chunksize': [0, validate_int],       # 0 to disable chunking;

    # key-mappings (multi-character mappings should be a list/tuple)
    'keymap.fullscreen':   [('f', 'ctrl+f'), validate_stringlist],
    'keymap.home':         [['h', 'r', 'home'], validate_stringlist],
    'keymap.back':         [['left', 'c', 'backspace'], validate_stringlist],
    'keymap.forward':      [['right', 'v'], validate_stringlist],
    'keymap.pan':          [['p'], validate_stringlist],
    'keymap.zoom':         [['o'], validate_stringlist],
    'keymap.save':         [['s', 'ctrl+s'], validate_stringlist],
    'keymap.quit':         [['ctrl+w', 'cmd+w', 'q'], validate_stringlist],
    'keymap.quit_all':     [['W', 'cmd+W', 'Q'], validate_stringlist],
    'keymap.grid':         [['g'], validate_stringlist],
    'keymap.grid_minor':   [['G'], validate_stringlist],
    'keymap.yscale':       [['l'], validate_stringlist],
    'keymap.xscale':       [['k', 'L'], validate_stringlist],
    'keymap.all_axes':     [['a'], validate_stringlist],

    # sample data
    'examples.directory': ['', validate_string],

    # Animation settings
    'animation.html':         ['none', validate_movie_html_fmt],
    # Limit, in MB, of size of base64 encoded animation in HTML
    # (i.e. IPython notebook)
    'animation.embed_limit':  [20, validate_float],
    'animation.writer':       ['ffmpeg', validate_movie_writer],
    'animation.codec':        ['h264', validate_string],
    'animation.bitrate':      [-1, validate_int],
    # Controls image format when frames are written to disk
    'animation.frame_format': ['png', validate_movie_frame_fmt],
    # Additional arguments for HTML writer
    'animation.html_args':    [[], validate_stringlist],
    # Path to ffmpeg binary. If just binary name, subprocess uses $PATH.
    'animation.ffmpeg_path':  ['ffmpeg', validate_animation_writer_path],
    # Additional arguments for ffmpeg movie writer (using pipes)
    'animation.ffmpeg_args':   [[], validate_stringlist],
    # Path to AVConv binary. If just binary name, subprocess uses $PATH.
    'animation.avconv_path':   ['avconv', validate_animation_writer_path],
    # Additional arguments for avconv movie writer (using pipes)
    'animation.avconv_args':   [[], validate_stringlist],
     # Path to convert binary. If just binary name, subprocess uses $PATH.
    'animation.convert_path':  ['convert', validate_animation_writer_path],
     # Additional arguments for convert movie writer (using pipes)
    'animation.convert_args':  [[], validate_stringlist],

    # Classic (pre 2.0) compatibility mode
    # This is used for things that are hard to make backward compatible
    # with a sane rcParam alone.  This does *not* turn on classic mode
    # altogether.  For that use `matplotlib.style.use('classic')`.
    '_internal.classic_mode': [False, validate_bool]
}


if __name__ == '__main__':
    rc = defaultParams
    rc['datapath'][0] = '/'
    for key in rc:
        if not rc[key][1](rc[key][0]) == rc[key][0]:
            print("%s: %s != %s" % (key, rc[key][1](rc[key][0]), rc[key][0]))
