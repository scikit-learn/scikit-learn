""" monkeypatching and mocking functionality.  """
from __future__ import absolute_import, division, print_function

import os
import sys
import re
import six
from _pytest.fixtures import fixture

RE_IMPORT_ERROR_NAME = re.compile("^No module named (.*)$")


@fixture
def monkeypatch():
    """The returned ``monkeypatch`` fixture provides these
    helper methods to modify objects, dictionaries or os.environ::

        monkeypatch.setattr(obj, name, value, raising=True)
        monkeypatch.delattr(obj, name, raising=True)
        monkeypatch.setitem(mapping, name, value)
        monkeypatch.delitem(obj, name, raising=True)
        monkeypatch.setenv(name, value, prepend=False)
        monkeypatch.delenv(name, value, raising=True)
        monkeypatch.syspath_prepend(path)
        monkeypatch.chdir(path)

    All modifications will be undone after the requesting
    test function or fixture has finished. The ``raising``
    parameter determines if a KeyError or AttributeError
    will be raised if the set/deletion operation has no target.
    """
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


def resolve(name):
    # simplified from zope.dottedname
    parts = name.split('.')

    used = parts.pop(0)
    found = __import__(used)
    for part in parts:
        used += '.' + part
        try:
            found = getattr(found, part)
        except AttributeError:
            pass
        else:
            continue
        # we use explicit un-nesting of the handling block in order
        # to avoid nested exceptions on python 3
        try:
            __import__(used)
        except ImportError as ex:
            # str is used for py2 vs py3
            expected = str(ex).split()[-1]
            if expected == used:
                raise
            else:
                raise ImportError(
                    'import error in %s: %s' % (used, ex)
                )
        found = annotated_getattr(found, part, used)
    return found


def annotated_getattr(obj, name, ann):
    try:
        obj = getattr(obj, name)
    except AttributeError:
        raise AttributeError(
            '%r object at %s has no attribute %r' % (
                type(obj).__name__, ann, name
            )
        )
    return obj


def derive_importpath(import_path, raising):
    if not isinstance(import_path, six.string_types) or "." not in import_path:
        raise TypeError("must be absolute import path string, not %r" %
                        (import_path,))
    module, attr = import_path.rsplit('.', 1)
    target = resolve(module)
    if raising:
        annotated_getattr(target, attr, ann=module)
    return attr, target


class Notset(object):
    def __repr__(self):
        return "<notset>"


notset = Notset()


class MonkeyPatch(object):
    """ Object returned by the ``monkeypatch`` fixture keeping a record of setattr/item/env/syspath changes.
    """

    def __init__(self):
        self._setattr = []
        self._setitem = []
        self._cwd = None
        self._savesyspath = None

    def setattr(self, target, name, value=notset, raising=True):
        """ Set attribute value on target, memorizing the old value.
        By default raise AttributeError if the attribute did not exist.

        For convenience you can specify a string as ``target`` which
        will be interpreted as a dotted import path, with the last part
        being the attribute name.  Example:
        ``monkeypatch.setattr("os.getcwd", lambda: "/")``
        would set the ``getcwd`` function of the ``os`` module.

        The ``raising`` value determines if the setattr should fail
        if the attribute is not already present (defaults to True
        which means it will raise).
        """
        __tracebackhide__ = True
        import inspect

        if value is notset:
            if not isinstance(target, six.string_types):
                raise TypeError("use setattr(target, name, value) or "
                                "setattr(target, value) with target being a dotted "
                                "import string")
            value = name
            name, target = derive_importpath(target, raising)

        oldval = getattr(target, name, notset)
        if raising and oldval is notset:
            raise AttributeError("%r has no attribute %r" % (target, name))

        # avoid class descriptors like staticmethod/classmethod
        if inspect.isclass(target):
            oldval = target.__dict__.get(name, notset)
        self._setattr.append((target, name, oldval))
        setattr(target, name, value)

    def delattr(self, target, name=notset, raising=True):
        """ Delete attribute ``name`` from ``target``, by default raise
        AttributeError it the attribute did not previously exist.

        If no ``name`` is specified and ``target`` is a string
        it will be interpreted as a dotted import path with the
        last part being the attribute name.

        If ``raising`` is set to False, no exception will be raised if the
        attribute is missing.
        """
        __tracebackhide__ = True
        if name is notset:
            if not isinstance(target, six.string_types):
                raise TypeError("use delattr(target, name) or "
                                "delattr(target) with target being a dotted "
                                "import string")
            name, target = derive_importpath(target, raising)

        if not hasattr(target, name):
            if raising:
                raise AttributeError(name)
        else:
            self._setattr.append((target, name, getattr(target, name, notset)))
            delattr(target, name)

    def setitem(self, dic, name, value):
        """ Set dictionary entry ``name`` to value. """
        self._setitem.append((dic, name, dic.get(name, notset)))
        dic[name] = value

    def delitem(self, dic, name, raising=True):
        """ Delete ``name`` from dict. Raise KeyError if it doesn't exist.

        If ``raising`` is set to False, no exception will be raised if the
        key is missing.
        """
        if name not in dic:
            if raising:
                raise KeyError(name)
        else:
            self._setitem.append((dic, name, dic.get(name, notset)))
            del dic[name]

    def setenv(self, name, value, prepend=None):
        """ Set environment variable ``name`` to ``value``.  If ``prepend``
        is a character, read the current environment variable value
        and prepend the ``value`` adjoined with the ``prepend`` character."""
        value = str(value)
        if prepend and name in os.environ:
            value = value + prepend + os.environ[name]
        self.setitem(os.environ, name, value)

    def delenv(self, name, raising=True):
        """ Delete ``name`` from the environment. Raise KeyError it does not
        exist.

        If ``raising`` is set to False, no exception will be raised if the
        environment variable is missing.
        """
        self.delitem(os.environ, name, raising=raising)

    def syspath_prepend(self, path):
        """ Prepend ``path`` to ``sys.path`` list of import locations. """
        if self._savesyspath is None:
            self._savesyspath = sys.path[:]
        sys.path.insert(0, str(path))

    def chdir(self, path):
        """ Change the current working directory to the specified path.
        Path can be a string or a py.path.local object.
        """
        if self._cwd is None:
            self._cwd = os.getcwd()
        if hasattr(path, "chdir"):
            path.chdir()
        else:
            os.chdir(path)

    def undo(self):
        """ Undo previous changes.  This call consumes the
        undo stack. Calling it a second time has no effect unless
        you do more monkeypatching after the undo call.

        There is generally no need to call `undo()`, since it is
        called automatically during tear-down.

        Note that the same `monkeypatch` fixture is used across a
        single test function invocation. If `monkeypatch` is used both by
        the test function itself and one of the test fixtures,
        calling `undo()` will undo all of the changes made in
        both functions.
        """
        for obj, name, value in reversed(self._setattr):
            if value is not notset:
                setattr(obj, name, value)
            else:
                delattr(obj, name)
        self._setattr[:] = []
        for dictionary, name, value in reversed(self._setitem):
            if value is notset:
                try:
                    del dictionary[name]
                except KeyError:
                    pass  # was already deleted, so we have the desired state
            else:
                dictionary[name] = value
        self._setitem[:] = []
        if self._savesyspath is not None:
            sys.path[:] = self._savesyspath
            self._savesyspath = None

        if self._cwd is not None:
            os.chdir(self._cwd)
            self._cwd = None
