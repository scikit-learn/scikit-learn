from collections import OrderedDict
import re

import numpy as np

from sklearn.base import BaseEstimator


def get_params(func):
    """Return dict (name: default_value) describing the parameters of given
    function."""
    try:  # Python 3
        from inspect import signature
        params = signature(func).parameters
        names = list(sorted(params.keys()))
        defaults = [params[name].default for name in names]
        params = OrderedDict(zip(names, defaults))
    except ImportError:  # Python 2
        from inspect import getargspec
        from itertools import izip_longest
        arg_specs = getargspec(func)
        if arg_specs.args[0] == 'self':
            names = arg_specs.args[1:]
        else:
            names = arg_specs.args
        names, defaults = zip(*sorted(izip_longest(names, arg_specs.defaults,
                                                   fillvalue=None)))
        params = OrderedDict(zip(names, defaults))

    return params


def changed_params(estimator):
    """Return dict (name: value) of parameters that were given to estimator
    with non-default values."""

    params = estimator.get_params(deep=False)
    filtered_params = {}
    init = getattr(estimator.__init__, 'deprecated_original',
                   estimator.__init__)
    init_params = get_params(init)
    for k, v in params.items():
        if v != init_params[k]:
            filtered_params[k] = v
    return filtered_params


def _strip_color(string):
    """Remove terminal color characters from a string"""
    return re.sub('\033\[[0-9;]+m', "", string)


class _Formatter(object):
    """Formatter class for pretty printing scikit-learn objects.

    Parameters
    ----------
    indent_est : str, optional (default='step')
        Indentation strategy for long strings describing estimators. If 'step',
        the indentation is constant and the next line is endented with 4 blank
        spaces. If 'name', the next line is aligned after the estimator name.
    changed_only : bool, optional (default=False)
        If True, only show parameters that have non-default values.
    color_changed : bool, optional (default=False)
        If True, color the parameters with non-default values.
    default_color : str, optional
        Color for the parameters with non-default values. Default is light grey
        (r=100,g=100,b=100).
    """
    def __init__(self, indent_est='step', changed_only=False,
                 color_changed=False, default_color='\033[38;2;100;100;100m'):
        self.indent_est = indent_est
        self.changed_only = changed_only
        self.color_changed = color_changed

        self.default_color = default_color
        self.types = {}
        self.htchar = ' '
        self.lfchar = '\n'
        self.step = 4
        self.width = 79

        self.set_formatter(object, self.__class__._format_object)
        self.set_formatter(dict, self.__class__._format_dict)
        self.set_formatter(list, self.__class__._format_list)
        self.set_formatter(tuple, self.__class__._format_tuple)
        self.set_formatter(BaseEstimator, self.__class__._format_estimator)
        self.set_formatter(np.ndarray, self.__class__._format_ndarray)
        self.set_formatter('callable', self.__class__._format_callable)

    def set_formatter(self, cls, callback):
        """Associate given callback with the formatting of objects with class
        cls"""
        self.types[cls] = callback

    def __call__(self, value):
        return self._format_all(value, indent=0)

    def _format_all(self, value, indent):
        """Return formated string for object value"""

        if type(value) in self.types:
            type_ = type(value)
        elif isinstance(value, BaseEstimator):
            type_ = BaseEstimator
        elif callable(value):
            type_ = 'callable'
        else:
            type_ = object

        formatter = self.types[type_]

        return formatter(self, value, indent)

    def _format_object(self, value, indent):
        return repr(value)

    def _format_dict(self, value, indent):
        items = [repr(key) + ': ' + self._format_all(value[key], indent +
                                                     self.step)
                 for key in sorted(value.keys())]
        return "{%s}" % self._join_items(items, indent + self.step)

    def _format_list(self, value, indent):
        items = [self._format_all(item, indent + self.step) for item in value]
        return "[%s]" % self._join_items(items, indent + self.step)

    def _format_tuple(self, value, indent):
        items = [self._format_all(item, indent + self.step) for item in value]
        return "(%s)" % self._join_items(items, indent + self.step)

    def _format_ndarray(self, value, indent):
        return self._format_all(value.tolist(), indent + self.step)

    def _format_callable(self, value, indent):
        # Used for functions and classes
        return value.__name__

    def _format_estimator(self, value, indent):
        if self.changed_only:
            params = changed_params(value)
        else:
            params = value.get_params(deep=False)

        # Sort params for version consistency
        params = OrderedDict((key, params[key]) for key in
                             sorted(params.keys()))

        if self.indent_est == 'step':
            offset = self.step
        elif self.indent_est == "name":
            offset = len(value.__class__.__name__) + 1
        else:
            raise ValueError("Invalid indent_est parameter.")

        # steps representation (only for Pipeline object)
        steps = params.pop("steps", None)
        if steps is not None:
            steps_str = (self.lfchar + self.htchar * (indent + offset) +
                         "steps=[")
            for i, step in enumerate(steps):
                if i > 0:
                    steps_str += (self.lfchar + self.htchar *
                                  (indent + len("steps=[") + offset))
                items = [
                    self._format_all(item, indent + offset +
                                     len("steps=[") + 1)
                    for item in step
                ]
                steps_str += "(%s)" % self._join_items(
                    items, indent + offset + len("steps=[") + 1)
            steps_str += "]," + self.lfchar + self.htchar * (indent + offset)
            init_indent = indent + offset
        else:
            steps_str = ""
            init_indent = indent + len(value.__class__.__name__) + 1
            # + 1 because of opening '('

        # Param representation
        items = [str(key) + '=' + self._format_all(params[key], indent +
                                                   offset)
                 for key in params]
        # add colors if needed
        if self.color_changed:
            changed = changed_params(value)

            def color(string, key):
                if key in changed:
                    return self.default_color + string + '\033[0m'
                else:
                    return string

            items = [color(item, key) for (item, key) in zip(items, params)]
        param_repr = self._join_items(items, indent + offset, init_indent)

        return '%s(%s)' % (value.__class__.__name__, steps_str + param_repr)

    def _join_items(self, items, indent, init_indent=None):
        # This method is used to separate items (typically parameter
        # lists) with commas and to break long lines appropriately
        this_string = ""
        init_indent = init_indent or indent
        pos = len(self.htchar) * (init_indent + 1)
        for i, item in enumerate(items):
            if i > 0:
                this_string += ','
                pos += 1
            if pos + len(_strip_color(item)) + 1 > self.width:
                # + 1 because we'll need the ',' at the end
                this_string += self.lfchar + self.htchar * indent
                pos = len(self.htchar) * (indent + 1)
            elif i > 0:
                this_string += " "
                pos += 1
            this_string += item
            pos += len(_strip_color(item))
        return this_string
