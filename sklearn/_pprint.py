from sklearn.base import BaseEstimator
from sklearn.pipeline import Pipeline
from inspect import signature
import numpy as np


def changed_params(estimator):
    params = estimator.get_params(deep=False)
    filtered_params = {}
    init = getattr(estimator.__init__, 'deprecated_original',
                   estimator.__init__)
    init_params = signature(init).parameters
    for k, v in params.items():
        if v != init_params[k].default:
            filtered_params[k] = v
    return filtered_params


class _Formatter(object):
    def __init__(self, indent_est='step', changed_only=False,
                 color_changed=False):
        self.indent_est = indent_est
        self.changed_only = changed_only
        self.color_changed = color_changed
        self.types = {}
        self.htchar = ' '
        self.lfchar = '\n'
        self.indent = 0
        self.step = 4
        self.width = 79
        self.set_formater(object, self.__class__.format_object)
        self.set_formater(dict, self.__class__.format_dict)
        self.set_formater(list, self.__class__.format_list)
        self.set_formater(tuple, self.__class__.format_tuple)
        self.set_formater(BaseEstimator, self.__class__.format_estimator)
        self.set_formater(np.ndarray, self.__class__.format_ndarray)
        self.set_formater(Pipeline, self.__class__.format_pipeline)

    def set_formater(self, obj, callback):
        self.types[obj] = callback

    def __call__(self, value, **args):
        for key in args:
            setattr(self, key, args[key])

        return self.format_all(value, self.indent)

    def format_all(self, value, indent):
        formater = self.types[type(value) if type(value) in self.types else
                              object]
        if type(value) in self.types:
            formater = self.types[type(value)]
        elif isinstance(value, BaseEstimator):
            formater = self.types[BaseEstimator]
        return formater(self, value, indent)

    def format_object(self, value, indent):
        return repr(value)

    def format_dict(self, value, indent):
        items = [repr(key) + ': ' + self.format_all(value[key],
                                                    indent + self.step)
                 for key in value]
        return "{%s}" % self.join_items(items, indent + self.step)

    def format_list(self, value, indent):
        items = [self.format_all(item, indent + self.step) for item in value]
        return "[%s]" % self.join_items(items, indent + self.step)

    def format_tuple(self, value, indent):
        items = [self.format_all(item, indent + self.step) for item in value]
        return "(%s)" % self.join_items(items, indent + self.step)

    def format_ndarray(self, value, indent):
        return self.format_all(value.tolist(), indent + self.step)

    def format_estimator(self, value, indent):
        if self.changed_only:
            params = changed_params(value)
        else:
            params = value.get_params(deep=False)

        if self.indent_est == 'step':
            offset = self.step
        elif self.indent_est == "name":
            offset = len(value.__class__.__name__) + 1
        else:
            raise ValueError("Invalid indent_est parameter")

        if self.color_changed:
            changed = changed_params(value)
            c = '\033[94m'

            def color(string, key):
                if key in changed:
                    return c + string + '\033[0m'
                else:
                    return string

            items = [color(str(key), key) + '='
                     + self.format_all(params[key], indent + offset)
                     for key in params]
        else:
            items = [str(key) + '='
                     + self.format_all(params[key], indent + offset)
                     for key in params]
        param_repr = self.join_items(items, indent + offset)
        return '%s(%s)' % (value.__class__.__name__, param_repr)

    def format_pipeline(self, value, indent):
        if self.changed_only:
            params = changed_params(value)
        else:
            params = value.get_params(deep=False)
        if self.indent_est == 'step':
            offset = self.step
        elif self.indent_est == "name":
            offset = len(value.__class__.__name__) + 1
        else:
            raise ValueError("Invalid indent_est parameter")
        steps = params.pop("steps")
        steps_str = self.lfchar + self.htchar * (indent + offset) + "steps=["
        for i, step in enumerate(steps):
            if i > 0:
                steps_str += (
                    self.lfchar
                    + self.htchar * (indent + len("steps=[") + offset))
            items = [self.format_all(item, indent + len("steps=[") + offset)
                     for item in step]
            steps_str += "(%s)" % self.join_items(
                items, indent + len("steps=[") + 1 + offset)
        steps_str += "]" + self.lfchar + self.htchar * (indent + offset)
        items = [str(key) + '=' + self.format_all(params[key], indent + offset)
                 for key in params]
        param_repr = self.join_items(items, indent + offset)
        return '%s(%s)' % (value.__class__.__name__, steps_str + param_repr)

    def join_items(self, items, indent):
        this_string = ""
        pos = len(self.htchar * (indent + 1)) + 2
        for i, item in enumerate(items):
            if i > 0:
                this_string += ','
            if pos + len(item) + 1 > self.width:
                this_string += self.lfchar + self.htchar * indent
                pos = len(self.htchar * (indent + 1)) + 1
            elif i > 0:
                this_string += " "
            this_string += item
            pos += len(item) + 1
        return this_string
