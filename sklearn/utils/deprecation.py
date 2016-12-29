import inspect
import warnings

__all__ = ["deprecated", ]


class deprecated(object):
    """Decorator to mark a function or class as deprecated.

    Issue a warning when the function is called/the class is instantiated and
    adds a warning to the docstring.

    The optional extra argument will be appended to the deprecation message
    and the docstring. Note: to use this with the default value for extra, put
    in an empty of parentheses:

    >>> from sklearn.utils import deprecated
    >>> deprecated() # doctest: +ELLIPSIS
    <sklearn.utils.deprecation.deprecated object at ...>

    >>> @deprecated()
    ... def some_function(): pass
    """

    # Adapted from http://wiki.python.org/moin/PythonDecoratorLibrary,
    # but with many changes.

    def __init__(self, extra=''):
        """
        Parameters
        ----------
        extra : string
          to be added to the deprecation messages

        """
        self.extra = extra

    def __call__(self, obj):
        if isinstance(obj, type):
            return self._decorate_class(obj)
        else:
            return self._decorate_fun(obj)

    def _decorate_class(self, cls):
        msg = "Class %s is deprecated" % cls.__name__
        if self.extra:
            msg += "; %s" % self.extra

        # FIXME: we should probably reset __new__ for full generality
        init = cls.__init__

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return init(*args, **kwargs)

        cls.__init__ = wrapped

        wrapped.__name__ = '__init__'
        wrapped.__doc__ = self._update_doc(init.__doc__)
        wrapped.deprecated_original = init

        return cls

    def _decorate_fun(self, fun):
        """Decorate function fun"""

        def wrapped(*args, **kwargs):

            def ismethod(f):
                spec = inspect.getargspec(f)
                return 'self' in spec.args

            name = "Function %s" % fun.__name__
            # if function is a method check if it is connected to an attribute
            if ismethod(fun):
                m = {}
                props = [x for x in inspect.getmembers(
                    args[0].__class__, lambda x: isinstance(x, property))]

                for prop_name, prop in props:
                    m[prop_name] = []
                    for (name, func) in inspect.getmembers(prop):
                        if name in ['fget', 'fset', 'fdel']:

                            if func is not None:
                                m[prop_name] += [func.__name__]

                for k, v in m.items():
                    if fun.__name__ in v:
                        name = "Attribute %s" % k
                        break

            msg = " %s is deprecated" % name
            if self.extra:
                msg += "; %s" % self.extra

            warnings.warn(msg, category=DeprecationWarning)
            return fun(*args, **kwargs)

        wrapped.__name__ = fun.__name__
        wrapped.__dict__ = fun.__dict__
        wrapped.__doc__ = self._update_doc(fun.__doc__)

        return wrapped

    def _update_doc(self, olddoc):
        newdoc = "DEPRECATED"
        if self.extra:
            newdoc = "%s: %s" % (newdoc, self.extra)
        if olddoc:
            newdoc = "%s\n\n%s" % (newdoc, olddoc)
        return newdoc
