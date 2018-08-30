"""
This class is defined to override standard pickle functionality

The goals of it follow:
-Serialize lambdas and nested functions to compiled byte code
-Deal with main module correctly
-Deal with other non-serializable objects

It does not include an unpickler, as standard python unpickling suffices.

This module was extracted from the `cloud` package, developed by `PiCloud, Inc.
<https://web.archive.org/web/20140626004012/http://www.picloud.com/>`_.

Copyright (c) 2012, Regents of the University of California.
Copyright (c) 2009 `PiCloud, Inc. <https://web.archive.org/web/20140626004012/http://www.picloud.com/>`_.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the University of California, Berkeley nor the
      names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from __future__ import print_function

import dis
from functools import partial
import imp
import io
import itertools
import logging
import opcode
import operator
import pickle
import struct
import sys
import traceback
import types
import weakref


# cloudpickle is meant for inter process communication: we expect all
# communicating processes to run the same Python version hence we favor
# communication speed over compatibility:
DEFAULT_PROTOCOL = pickle.HIGHEST_PROTOCOL


if sys.version < '3':
    from pickle import Pickler
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
    PY3 = False
else:
    types.ClassType = type
    from pickle import _Pickler as Pickler
    from io import BytesIO as StringIO
    PY3 = True


def _make_cell_set_template_code():
    """Get the Python compiler to emit LOAD_FAST(arg); STORE_DEREF

    Notes
    -----
    In Python 3, we could use an easier function:

    .. code-block:: python

       def f():
           cell = None

           def _stub(value):
               nonlocal cell
               cell = value

           return _stub

        _cell_set_template_code = f()

    This function is _only_ a LOAD_FAST(arg); STORE_DEREF, but that is
    invalid syntax on Python 2. If we use this function we also don't need
    to do the weird freevars/cellvars swap below
    """
    def inner(value):
        lambda: cell  # make ``cell`` a closure so that we get a STORE_DEREF
        cell = value

    co = inner.__code__

    # NOTE: we are marking the cell variable as a free variable intentionally
    # so that we simulate an inner function instead of the outer function. This
    # is what gives us the ``nonlocal`` behavior in a Python 2 compatible way.
    if not PY3:
        return types.CodeType(
            co.co_argcount,
            co.co_nlocals,
            co.co_stacksize,
            co.co_flags,
            co.co_code,
            co.co_consts,
            co.co_names,
            co.co_varnames,
            co.co_filename,
            co.co_name,
            co.co_firstlineno,
            co.co_lnotab,
            co.co_cellvars,  # this is the trickery
            (),
        )
    else:
        return types.CodeType(
            co.co_argcount,
            co.co_kwonlyargcount,
            co.co_nlocals,
            co.co_stacksize,
            co.co_flags,
            co.co_code,
            co.co_consts,
            co.co_names,
            co.co_varnames,
            co.co_filename,
            co.co_name,
            co.co_firstlineno,
            co.co_lnotab,
            co.co_cellvars,  # this is the trickery
            (),
        )


_cell_set_template_code = _make_cell_set_template_code()


def cell_set(cell, value):
    """Set the value of a closure cell.
    """
    return types.FunctionType(
        _cell_set_template_code,
        {},
        '_cell_set_inner',
        (),
        (cell,),
    )(value)


# relevant opcodes
STORE_GLOBAL = opcode.opmap['STORE_GLOBAL']
DELETE_GLOBAL = opcode.opmap['DELETE_GLOBAL']
LOAD_GLOBAL = opcode.opmap['LOAD_GLOBAL']
GLOBAL_OPS = (STORE_GLOBAL, DELETE_GLOBAL, LOAD_GLOBAL)
HAVE_ARGUMENT = dis.HAVE_ARGUMENT
EXTENDED_ARG = dis.EXTENDED_ARG


def islambda(func):
    return getattr(func, '__name__') == '<lambda>'


_BUILTIN_TYPE_NAMES = {}
for k, v in types.__dict__.items():
    if type(v) is type:
        _BUILTIN_TYPE_NAMES[v] = k


def _builtin_type(name):
    return getattr(types, name)


def _make__new__factory(type_):
    def _factory():
        return type_.__new__
    return _factory


# NOTE: These need to be module globals so that they're pickleable as globals.
_get_dict_new = _make__new__factory(dict)
_get_frozenset_new = _make__new__factory(frozenset)
_get_list_new = _make__new__factory(list)
_get_set_new = _make__new__factory(set)
_get_tuple_new = _make__new__factory(tuple)
_get_object_new = _make__new__factory(object)

# Pre-defined set of builtin_function_or_method instances that can be
# serialized.
_BUILTIN_TYPE_CONSTRUCTORS = {
    dict.__new__: _get_dict_new,
    frozenset.__new__: _get_frozenset_new,
    set.__new__: _get_set_new,
    list.__new__: _get_list_new,
    tuple.__new__: _get_tuple_new,
    object.__new__: _get_object_new,
}


if sys.version_info < (3, 4):
    def _walk_global_ops(code):
        """
        Yield (opcode, argument number) tuples for all
        global-referencing instructions in *code*.
        """
        code = getattr(code, 'co_code', b'')
        if not PY3:
            code = map(ord, code)

        n = len(code)
        i = 0
        extended_arg = 0
        while i < n:
            op = code[i]
            i += 1
            if op >= HAVE_ARGUMENT:
                oparg = code[i] + code[i + 1] * 256 + extended_arg
                extended_arg = 0
                i += 2
                if op == EXTENDED_ARG:
                    extended_arg = oparg * 65536
                if op in GLOBAL_OPS:
                    yield op, oparg

else:
    def _walk_global_ops(code):
        """
        Yield (opcode, argument number) tuples for all
        global-referencing instructions in *code*.
        """
        for instr in dis.get_instructions(code):
            op = instr.opcode
            if op in GLOBAL_OPS:
                yield op, instr.arg


class CloudPickler(Pickler):

    dispatch = Pickler.dispatch.copy()

    def __init__(self, file, protocol=None):
        if protocol is None:
            protocol = DEFAULT_PROTOCOL
        Pickler.__init__(self, file, protocol=protocol)
        # set of modules to unpickle
        self.modules = set()
        # map ids to dictionary. used to ensure that functions can share global env
        self.globals_ref = {}

    def dump(self, obj):
        self.inject_addons()
        try:
            return Pickler.dump(self, obj)
        except RuntimeError as e:
            if 'recursion' in e.args[0]:
                msg = """Could not pickle object as excessively deep recursion required."""
                raise pickle.PicklingError(msg)
            else:
                raise

    def save_memoryview(self, obj):
        self.save(obj.tobytes())

    dispatch[memoryview] = save_memoryview

    if not PY3:
        def save_buffer(self, obj):
            self.save(str(obj))

        dispatch[buffer] = save_buffer  # noqa: F821 'buffer' was removed in Python 3

    def save_module(self, obj):
        """
        Save a module as an import
        """
        mod_name = obj.__name__
        # If module is successfully found then it is not a dynamically created module
        if hasattr(obj, '__file__'):
            is_dynamic = False
        else:
            try:
                _find_module(mod_name)
                is_dynamic = False
            except ImportError:
                is_dynamic = True

        self.modules.add(obj)
        if is_dynamic:
            self.save_reduce(dynamic_subimport, (obj.__name__, vars(obj)), obj=obj)
        else:
            self.save_reduce(subimport, (obj.__name__,), obj=obj)

    dispatch[types.ModuleType] = save_module

    def save_codeobject(self, obj):
        """
        Save a code object
        """
        if PY3:
            args = (
                obj.co_argcount, obj.co_kwonlyargcount, obj.co_nlocals, obj.co_stacksize,
                obj.co_flags, obj.co_code, obj.co_consts, obj.co_names, obj.co_varnames,
                obj.co_filename, obj.co_name, obj.co_firstlineno, obj.co_lnotab, obj.co_freevars,
                obj.co_cellvars
            )
        else:
            args = (
                obj.co_argcount, obj.co_nlocals, obj.co_stacksize, obj.co_flags, obj.co_code,
                obj.co_consts, obj.co_names, obj.co_varnames, obj.co_filename, obj.co_name,
                obj.co_firstlineno, obj.co_lnotab, obj.co_freevars, obj.co_cellvars
            )
        self.save_reduce(types.CodeType, args, obj=obj)

    dispatch[types.CodeType] = save_codeobject

    def save_function(self, obj, name=None):
        """ Registered with the dispatch to handle all function types.

        Determines what kind of function obj is (e.g. lambda, defined at
        interactive prompt, etc) and handles the pickling appropriately.
        """
        try:
            should_special_case = obj in _BUILTIN_TYPE_CONSTRUCTORS
        except TypeError:
            # Methods of builtin types aren't hashable in python 2.
            should_special_case = False

        if should_special_case:
            # We keep a special-cased cache of built-in type constructors at
            # global scope, because these functions are structured very
            # differently in different python versions and implementations (for
            # example, they're instances of types.BuiltinFunctionType in
            # CPython, but they're ordinary types.FunctionType instances in
            # PyPy).
            #
            # If the function we've received is in that cache, we just
            # serialize it as a lookup into the cache.
            return self.save_reduce(_BUILTIN_TYPE_CONSTRUCTORS[obj], (), obj=obj)

        write = self.write

        if name is None:
            name = obj.__name__
        try:
            # whichmodule() could fail, see
            # https://bitbucket.org/gutworth/six/issues/63/importing-six-breaks-pickling
            modname = pickle.whichmodule(obj, name)
        except Exception:
            modname = None
        # print('which gives %s %s %s' % (modname, obj, name))
        try:
            themodule = sys.modules[modname]
        except KeyError:
            # eval'd items such as namedtuple give invalid items for their function __module__
            modname = '__main__'

        if modname == '__main__':
            themodule = None

        try:
            lookedup_by_name = getattr(themodule, name, None)
        except Exception:
            lookedup_by_name = None

        if themodule:
            self.modules.add(themodule)
            if lookedup_by_name is obj:
                return self.save_global(obj, name)

        # a builtin_function_or_method which comes in as an attribute of some
        # object (e.g., itertools.chain.from_iterable) will end
        # up with modname "__main__" and so end up here. But these functions
        # have no __code__ attribute in CPython, so the handling for
        # user-defined functions below will fail.
        # So we pickle them here using save_reduce; have to do it differently
        # for different python versions.
        if not hasattr(obj, '__code__'):
            if PY3:
                rv = obj.__reduce_ex__(self.proto)
            else:
                if hasattr(obj, '__self__'):
                    rv = (getattr, (obj.__self__, name))
                else:
                    raise pickle.PicklingError("Can't pickle %r" % obj)
            return self.save_reduce(obj=obj, *rv)

        # if func is lambda, def'ed at prompt, is in main, or is nested, then
        # we'll pickle the actual function object rather than simply saving a
        # reference (as is done in default pickler), via save_function_tuple.
        if (islambda(obj)
                or getattr(obj.__code__, 'co_filename', None) == '<stdin>'
                or themodule is None):
            self.save_function_tuple(obj)
            return
        else:
            # func is nested
            if lookedup_by_name is None or lookedup_by_name is not obj:
                self.save_function_tuple(obj)
                return

        if obj.__dict__:
            # essentially save_reduce, but workaround needed to avoid recursion
            self.save(_restore_attr)
            write(pickle.MARK + pickle.GLOBAL + modname + '\n' + name + '\n')
            self.memoize(obj)
            self.save(obj.__dict__)
            write(pickle.TUPLE + pickle.REDUCE)
        else:
            write(pickle.GLOBAL + modname + '\n' + name + '\n')
            self.memoize(obj)

    dispatch[types.FunctionType] = save_function

    def _save_subimports(self, code, top_level_dependencies):
        """
        Ensure de-pickler imports any package child-modules that
        are needed by the function
        """

        # check if any known dependency is an imported package
        for x in top_level_dependencies:
            if isinstance(x, types.ModuleType) and hasattr(x, '__package__') and x.__package__:
                # check if the package has any currently loaded sub-imports
                prefix = x.__name__ + '.'
                # A concurrent thread could mutate sys.modules,
                # make sure we iterate over a copy to avoid exceptions
                for name in list(sys.modules):
                    # Older versions of pytest will add a "None" module to sys.modules.
                    if name is not None and name.startswith(prefix):
                        # check whether the function can address the sub-module
                        tokens = set(name[len(prefix):].split('.'))
                        if not tokens - set(code.co_names):
                            # ensure unpickler executes this import
                            self.save(sys.modules[name])
                            # then discards the reference to it
                            self.write(pickle.POP)

    def save_dynamic_class(self, obj):
        """
        Save a class that can't be stored as module global.

        This method is used to serialize classes that are defined inside
        functions, or that otherwise can't be serialized as attribute lookups
        from global modules.
        """
        clsdict = dict(obj.__dict__)  # copy dict proxy to a dict
        clsdict.pop('__weakref__', None)

        # For ABCMeta in python3.7+, remove _abc_impl as it is not picklable.
        # This is a fix which breaks the cache but this only makes the first
        # calls to issubclass slower.
        if "_abc_impl" in clsdict:
            import abc
            (registry, _, _, _) = abc._get_dump(obj)
            clsdict["_abc_impl"] = [subclass_weakref()
                                    for subclass_weakref in registry]

        # On PyPy, __doc__ is a readonly attribute, so we need to include it in
        # the initial skeleton class.  This is safe because we know that the
        # doc can't participate in a cycle with the original class.
        type_kwargs = {'__doc__': clsdict.pop('__doc__', None)}

        # If type overrides __dict__ as a property, include it in the type kwargs.
        # In Python 2, we can't set this attribute after construction.
        __dict__ = clsdict.pop('__dict__', None)
        if isinstance(__dict__, property):
            type_kwargs['__dict__'] = __dict__

        save = self.save
        write = self.write

        # We write pickle instructions explicitly here to handle the
        # possibility that the type object participates in a cycle with its own
        # __dict__. We first write an empty "skeleton" version of the class and
        # memoize it before writing the class' __dict__ itself. We then write
        # instructions to "rehydrate" the skeleton class by restoring the
        # attributes from the __dict__.
        #
        # A type can appear in a cycle with its __dict__ if an instance of the
        # type appears in the type's __dict__ (which happens for the stdlib
        # Enum class), or if the type defines methods that close over the name
        # of the type, (which is common for Python 2-style super() calls).

        # Push the rehydration function.
        save(_rehydrate_skeleton_class)

        # Mark the start of the args tuple for the rehydration function.
        write(pickle.MARK)

        # Create and memoize an skeleton class with obj's name and bases.
        tp = type(obj)
        self.save_reduce(tp, (obj.__name__, obj.__bases__, type_kwargs), obj=obj)

        # Now save the rest of obj's __dict__. Any references to obj
        # encountered while saving will point to the skeleton class.
        save(clsdict)

        # Write a tuple of (skeleton_class, clsdict).
        write(pickle.TUPLE)

        # Call _rehydrate_skeleton_class(skeleton_class, clsdict)
        write(pickle.REDUCE)

    def save_function_tuple(self, func):
        """  Pickles an actual func object.

        A func comprises: code, globals, defaults, closure, and dict.  We
        extract and save these, injecting reducing functions at certain points
        to recreate the func object.  Keep in mind that some of these pieces
        can contain a ref to the func itself.  Thus, a naive save on these
        pieces could trigger an infinite loop of save's.  To get around that,
        we first create a skeleton func object using just the code (this is
        safe, since this won't contain a ref to the func), and memoize it as
        soon as it's created.  The other stuff can then be filled in later.
        """
        if is_tornado_coroutine(func):
            self.save_reduce(_rebuild_tornado_coroutine, (func.__wrapped__,),
                             obj=func)
            return

        save = self.save
        write = self.write

        code, f_globals, defaults, closure_values, dct, base_globals = self.extract_func_data(func)

        save(_fill_function)  # skeleton function updater
        write(pickle.MARK)    # beginning of tuple that _fill_function expects

        self._save_subimports(
            code,
            itertools.chain(f_globals.values(), closure_values or ()),
        )

        # create a skeleton function object and memoize it
        save(_make_skel_func)
        save((
            code,
            len(closure_values) if closure_values is not None else -1,
            base_globals,
        ))
        write(pickle.REDUCE)
        self.memoize(func)

        # save the rest of the func data needed by _fill_function
        state = {
            'globals': f_globals,
            'defaults': defaults,
            'dict': dct,
            'closure_values': closure_values,
            'module': func.__module__,
            'name': func.__name__,
            'doc': func.__doc__,
        }
        if hasattr(func, '__annotations__'):
            state['annotations'] = func.__annotations__
        if hasattr(func, '__qualname__'):
            state['qualname'] = func.__qualname__
        save(state)
        write(pickle.TUPLE)
        write(pickle.REDUCE)  # applies _fill_function on the tuple

    _extract_code_globals_cache = (
        weakref.WeakKeyDictionary()
        if not hasattr(sys, "pypy_version_info")
        else {})

    @classmethod
    def extract_code_globals(cls, co):
        """
        Find all globals names read or written to by codeblock co
        """
        out_names = cls._extract_code_globals_cache.get(co)
        if out_names is None:
            try:
                names = co.co_names
            except AttributeError:
                # PyPy "builtin-code" object
                out_names = set()
            else:
                out_names = {names[oparg] for _, oparg in _walk_global_ops(co)}

                # see if nested function have any global refs
                if co.co_consts:
                    for const in co.co_consts:
                        if type(const) is types.CodeType:
                            out_names |= cls.extract_code_globals(const)

            cls._extract_code_globals_cache[co] = out_names

        return out_names

    def extract_func_data(self, func):
        """
        Turn the function into a tuple of data necessary to recreate it:
            code, globals, defaults, closure_values, dict
        """
        code = func.__code__

        # extract all global ref's
        func_global_refs = self.extract_code_globals(code)

        # process all variables referenced by global environment
        f_globals = {}
        for var in func_global_refs:
            if var in func.__globals__:
                f_globals[var] = func.__globals__[var]

        # defaults requires no processing
        defaults = func.__defaults__

        # process closure
        closure = (
            list(map(_get_cell_contents, func.__closure__))
            if func.__closure__ is not None
            else None
        )

        # save the dict
        dct = func.__dict__

        base_globals = self.globals_ref.get(id(func.__globals__), None)
        if base_globals is None:
            # For functions defined in __main__, use vars(__main__) for
            # base_global. This is necessary to share the global variables
            # across multiple functions in this module.
            if func.__module__ == "__main__":
                base_globals = "__main__"
            else:
                base_globals = {}
        self.globals_ref[id(func.__globals__)] = base_globals

        return (code, f_globals, defaults, closure, dct, base_globals)

    def save_builtin_function(self, obj):
        if obj.__module__ == "__builtin__":
            return self.save_global(obj)
        return self.save_function(obj)

    dispatch[types.BuiltinFunctionType] = save_builtin_function

    def save_global(self, obj, name=None, pack=struct.pack):
        """
        Save a "global".

        The name of this method is somewhat misleading: all types get
        dispatched here.
        """
        if obj.__module__ == "__main__":
            return self.save_dynamic_class(obj)

        try:
            return Pickler.save_global(self, obj, name=name)
        except Exception:
            if obj.__module__ == "__builtin__" or obj.__module__ == "builtins":
                if obj in _BUILTIN_TYPE_NAMES:
                    return self.save_reduce(
                        _builtin_type, (_BUILTIN_TYPE_NAMES[obj],), obj=obj)

            typ = type(obj)
            if typ is not obj and isinstance(obj, (type, types.ClassType)):
                return self.save_dynamic_class(obj)

            raise

    dispatch[type] = save_global
    dispatch[types.ClassType] = save_global

    def save_instancemethod(self, obj):
        # Memoization rarely is ever useful due to python bounding
        if obj.__self__ is None:
            self.save_reduce(getattr, (obj.im_class, obj.__name__))
        else:
            if PY3:
                self.save_reduce(types.MethodType, (obj.__func__, obj.__self__), obj=obj)
            else:
                self.save_reduce(types.MethodType, (obj.__func__, obj.__self__, obj.__self__.__class__),
                                 obj=obj)

    dispatch[types.MethodType] = save_instancemethod

    def save_inst(self, obj):
        """Inner logic to save instance. Based off pickle.save_inst"""
        cls = obj.__class__

        # Try the dispatch table (pickle module doesn't do it)
        f = self.dispatch.get(cls)
        if f:
            f(self, obj)  # Call unbound method with explicit self
            return

        memo = self.memo
        write = self.write
        save = self.save

        if hasattr(obj, '__getinitargs__'):
            args = obj.__getinitargs__()
            len(args)  # XXX Assert it's a sequence
            pickle._keep_alive(args, memo)
        else:
            args = ()

        write(pickle.MARK)

        if self.bin:
            save(cls)
            for arg in args:
                save(arg)
            write(pickle.OBJ)
        else:
            for arg in args:
                save(arg)
            write(pickle.INST + cls.__module__ + '\n' + cls.__name__ + '\n')

        self.memoize(obj)

        try:
            getstate = obj.__getstate__
        except AttributeError:
            stuff = obj.__dict__
        else:
            stuff = getstate()
            pickle._keep_alive(stuff, memo)
        save(stuff)
        write(pickle.BUILD)

    if not PY3:
        dispatch[types.InstanceType] = save_inst

    def save_property(self, obj):
        # properties not correctly saved in python
        self.save_reduce(property, (obj.fget, obj.fset, obj.fdel, obj.__doc__), obj=obj)

    dispatch[property] = save_property

    def save_classmethod(self, obj):
        orig_func = obj.__func__
        self.save_reduce(type(obj), (orig_func,), obj=obj)

    dispatch[classmethod] = save_classmethod
    dispatch[staticmethod] = save_classmethod

    def save_itemgetter(self, obj):
        """itemgetter serializer (needed for namedtuple support)"""
        class Dummy:
            def __getitem__(self, item):
                return item
        items = obj(Dummy())
        if not isinstance(items, tuple):
            items = (items,)
        return self.save_reduce(operator.itemgetter, items)

    if type(operator.itemgetter) is type:
        dispatch[operator.itemgetter] = save_itemgetter

    def save_attrgetter(self, obj):
        """attrgetter serializer"""
        class Dummy(object):
            def __init__(self, attrs, index=None):
                self.attrs = attrs
                self.index = index
            def __getattribute__(self, item):
                attrs = object.__getattribute__(self, "attrs")
                index = object.__getattribute__(self, "index")
                if index is None:
                    index = len(attrs)
                    attrs.append(item)
                else:
                    attrs[index] = ".".join([attrs[index], item])
                return type(self)(attrs, index)
        attrs = []
        obj(Dummy(attrs))
        return self.save_reduce(operator.attrgetter, tuple(attrs))

    if type(operator.attrgetter) is type:
        dispatch[operator.attrgetter] = save_attrgetter

    def save_file(self, obj):
        """Save a file"""
        try:
            import StringIO as pystringIO  # we can't use cStringIO as it lacks the name attribute
        except ImportError:
            import io as pystringIO

        if not hasattr(obj, 'name') or not hasattr(obj, 'mode'):
            raise pickle.PicklingError("Cannot pickle files that do not map to an actual file")
        if obj is sys.stdout:
            return self.save_reduce(getattr, (sys, 'stdout'), obj=obj)
        if obj is sys.stderr:
            return self.save_reduce(getattr, (sys, 'stderr'), obj=obj)
        if obj is sys.stdin:
            raise pickle.PicklingError("Cannot pickle standard input")
        if obj.closed:
            raise pickle.PicklingError("Cannot pickle closed files")
        if hasattr(obj, 'isatty') and obj.isatty():
            raise pickle.PicklingError("Cannot pickle files that map to tty objects")
        if 'r' not in obj.mode and '+' not in obj.mode:
            raise pickle.PicklingError("Cannot pickle files that are not opened for reading: %s" % obj.mode)

        name = obj.name

        retval = pystringIO.StringIO()

        try:
            # Read the whole file
            curloc = obj.tell()
            obj.seek(0)
            contents = obj.read()
            obj.seek(curloc)
        except IOError:
            raise pickle.PicklingError("Cannot pickle file %s as it cannot be read" % name)
        retval.write(contents)
        retval.seek(curloc)

        retval.name = name
        self.save(retval)
        self.memoize(obj)

    def save_ellipsis(self, obj):
        self.save_reduce(_gen_ellipsis, ())

    def save_not_implemented(self, obj):
        self.save_reduce(_gen_not_implemented, ())

    try:               # Python 2
        dispatch[file] = save_file
    except NameError:  # Python 3
        dispatch[io.TextIOWrapper] = save_file

    dispatch[type(Ellipsis)] = save_ellipsis
    dispatch[type(NotImplemented)] = save_not_implemented

    def save_weakset(self, obj):
        self.save_reduce(weakref.WeakSet, (list(obj),))

    dispatch[weakref.WeakSet] = save_weakset

    def save_logger(self, obj):
        self.save_reduce(logging.getLogger, (obj.name,), obj=obj)

    dispatch[logging.Logger] = save_logger

    def save_root_logger(self, obj):
        self.save_reduce(logging.getLogger, (), obj=obj)

    dispatch[logging.RootLogger] = save_root_logger

    """Special functions for Add-on libraries"""
    def inject_addons(self):
        """Plug in system. Register additional pickling functions if modules already loaded"""
        pass


# Tornado support

def is_tornado_coroutine(func):
    """
    Return whether *func* is a Tornado coroutine function.
    Running coroutines are not supported.
    """
    if 'tornado.gen' not in sys.modules:
        return False
    gen = sys.modules['tornado.gen']
    if not hasattr(gen, "is_coroutine_function"):
        # Tornado version is too old
        return False
    return gen.is_coroutine_function(func)


def _rebuild_tornado_coroutine(func):
    from tornado import gen
    return gen.coroutine(func)


# Shorthands for legacy support

def dump(obj, file, protocol=None):
    """Serialize obj as bytes streamed into file

    protocol defaults to cloudpickle.DEFAULT_PROTOCOL which is an alias to
    pickle.HIGHEST_PROTOCOL. This setting favors maximum communication speed
    between processes running the same Python version.

    Set protocol=pickle.DEFAULT_PROTOCOL instead if you need to ensure
    compatibility with older versions of Python.
    """
    CloudPickler(file, protocol=protocol).dump(obj)


def dumps(obj, protocol=None):
    """Serialize obj as a string of bytes allocated in memory

    protocol defaults to cloudpickle.DEFAULT_PROTOCOL which is an alias to
    pickle.HIGHEST_PROTOCOL. This setting favors maximum communication speed
    between processes running the same Python version.

    Set protocol=pickle.DEFAULT_PROTOCOL instead if you need to ensure
    compatibility with older versions of Python.
    """
    file = StringIO()
    try:
        cp = CloudPickler(file, protocol=protocol)
        cp.dump(obj)
        return file.getvalue()
    finally:
        file.close()


# including pickles unloading functions in this namespace
load = pickle.load
loads = pickle.loads


# hack for __import__ not working as desired
def subimport(name):
    __import__(name)
    return sys.modules[name]


def dynamic_subimport(name, vars):
    mod = imp.new_module(name)
    mod.__dict__.update(vars)
    sys.modules[name] = mod
    return mod


# restores function attributes
def _restore_attr(obj, attr):
    for key, val in attr.items():
        setattr(obj, key, val)
    return obj


def _get_module_builtins():
    return pickle.__builtins__


def print_exec(stream):
    ei = sys.exc_info()
    traceback.print_exception(ei[0], ei[1], ei[2], None, stream)


def _modules_to_main(modList):
    """Force every module in modList to be placed into main"""
    if not modList:
        return

    main = sys.modules['__main__']
    for modname in modList:
        if type(modname) is str:
            try:
                mod = __import__(modname)
            except Exception:
                sys.stderr.write('warning: could not import %s\n.  '
                                 'Your function may unexpectedly error due to this import failing;'
                                 'A version mismatch is likely.  Specific error was:\n' % modname)
                print_exec(sys.stderr)
            else:
                setattr(main, mod.__name__, mod)


# object generators:
def _genpartial(func, args, kwds):
    if not args:
        args = ()
    if not kwds:
        kwds = {}
    return partial(func, *args, **kwds)


def _gen_ellipsis():
    return Ellipsis


def _gen_not_implemented():
    return NotImplemented


def _get_cell_contents(cell):
    try:
        return cell.cell_contents
    except ValueError:
        # sentinel used by ``_fill_function`` which will leave the cell empty
        return _empty_cell_value


def instance(cls):
    """Create a new instance of a class.

    Parameters
    ----------
    cls : type
        The class to create an instance of.

    Returns
    -------
    instance : cls
        A new instance of ``cls``.
    """
    return cls()


@instance
class _empty_cell_value(object):
    """sentinel for empty closures
    """
    @classmethod
    def __reduce__(cls):
        return cls.__name__


def _fill_function(*args):
    """Fills in the rest of function data into the skeleton function object

    The skeleton itself is create by _make_skel_func().
    """
    if len(args) == 2:
        func = args[0]
        state = args[1]
    elif len(args) == 5:
        # Backwards compat for cloudpickle v0.4.0, after which the `module`
        # argument was introduced
        func = args[0]
        keys = ['globals', 'defaults', 'dict', 'closure_values']
        state = dict(zip(keys, args[1:]))
    elif len(args) == 6:
        # Backwards compat for cloudpickle v0.4.1, after which the function
        # state was passed as a dict to the _fill_function it-self.
        func = args[0]
        keys = ['globals', 'defaults', 'dict', 'module', 'closure_values']
        state = dict(zip(keys, args[1:]))
    else:
        raise ValueError('Unexpected _fill_value arguments: %r' % (args,))

    # Only set global variables that do not exist.
    for k, v in state['globals'].items():
        if k not in func.__globals__:
            func.__globals__[k] = v

    func.__defaults__ = state['defaults']
    func.__dict__ = state['dict']
    if 'annotations' in state:
        func.__annotations__ = state['annotations']
    if 'doc' in state:
        func.__doc__  = state['doc']
    if 'name' in state:
        func.__name__ = state['name']
    if 'module' in state:
        func.__module__ = state['module']
    if 'qualname' in state:
        func.__qualname__ = state['qualname']

    cells = func.__closure__
    if cells is not None:
        for cell, value in zip(cells, state['closure_values']):
            if value is not _empty_cell_value:
                cell_set(cell, value)

    return func


def _make_empty_cell():
    if False:
        # trick the compiler into creating an empty cell in our lambda
        cell = None
        raise AssertionError('this route should not be executed')

    return (lambda: cell).__closure__[0]


def _make_skel_func(code, cell_count, base_globals=None):
    """ Creates a skeleton function object that contains just the provided
        code and the correct number of cells in func_closure.  All other
        func attributes (e.g. func_globals) are empty.
    """
    if base_globals is None:
        base_globals = {}
    elif isinstance(base_globals, str):
        base_globals = vars(sys.modules[base_globals])
    base_globals['__builtins__'] = __builtins__

    closure = (
        tuple(_make_empty_cell() for _ in range(cell_count))
        if cell_count >= 0 else
        None
    )
    return types.FunctionType(code, base_globals, None, None, closure)


def _rehydrate_skeleton_class(skeleton_class, class_dict):
    """Put attributes from `class_dict` back on `skeleton_class`.

    See CloudPickler.save_dynamic_class for more info.
    """
    registry = None
    for attrname, attr in class_dict.items():
        if attrname == "_abc_impl":
            registry = attr
        else:
            setattr(skeleton_class, attrname, attr)
    if registry is not None:
        for subclass in registry:
            skeleton_class.register(subclass)

    return skeleton_class


def _find_module(mod_name):
    """
    Iterate over each part instead of calling imp.find_module directly.
    This function is able to find submodules (e.g. scikit.tree)
    """
    path = None
    for part in mod_name.split('.'):
        if path is not None:
            path = [path]
        file, path, description = imp.find_module(part, path)
        if file is not None:
            file.close()
    return path, description


"""Constructors for 3rd party libraries
Note: These can never be renamed due to client compatibility issues"""


def _getobject(modname, attribute):
    mod = __import__(modname, fromlist=[attribute])
    return mod.__dict__[attribute]


""" Use copy_reg to extend global pickle definitions """

if sys.version_info < (3, 4):
    method_descriptor = type(str.upper)

    def _reduce_method_descriptor(obj):
        return (getattr, (obj.__objclass__, obj.__name__))

    try:
        import copy_reg as copyreg
    except ImportError:
        import copyreg
    copyreg.pickle(method_descriptor, _reduce_method_descriptor)
