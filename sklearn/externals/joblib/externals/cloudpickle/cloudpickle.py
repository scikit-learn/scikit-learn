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
import opcode
import operator
import pickle
import struct
import sys
import traceback
import types
import weakref


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


#relevant opcodes
STORE_GLOBAL = opcode.opmap['STORE_GLOBAL']
DELETE_GLOBAL = opcode.opmap['DELETE_GLOBAL']
LOAD_GLOBAL = opcode.opmap['LOAD_GLOBAL']
GLOBAL_OPS = (STORE_GLOBAL, DELETE_GLOBAL, LOAD_GLOBAL)
HAVE_ARGUMENT = dis.HAVE_ARGUMENT
EXTENDED_ARG = dis.EXTENDED_ARG


def islambda(func):
    return getattr(func,'__name__') == '<lambda>'


_BUILTIN_TYPE_NAMES = {}
for k, v in types.__dict__.items():
    if type(v) is type:
        _BUILTIN_TYPE_NAMES[v] = k


def _builtin_type(name):
    return getattr(types, name)


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
        Pickler.__init__(self, file, protocol)
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

    def save_memoryview(self, obj):
        """Fallback to save_string"""
        Pickler.save_string(self, str(obj))

    def save_buffer(self, obj):
        """Fallback to save_string"""
        Pickler.save_string(self,str(obj))
    if PY3:
        dispatch[memoryview] = save_memoryview
    else:
        dispatch[buffer] = save_buffer

    def save_unsupported(self, obj):
        raise pickle.PicklingError("Cannot pickle objects of type %s" % type(obj))
    dispatch[types.GeneratorType] = save_unsupported

    # itertools objects do not pickle!
    for v in itertools.__dict__.values():
        if type(v) is type:
            dispatch[v] = save_unsupported

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
        write = self.write

        if name is None:
            name = obj.__name__
        modname = pickle.whichmodule(obj, name)
        # print('which gives %s %s %s' % (modname, obj, name))
        try:
            themodule = sys.modules[modname]
        except KeyError:
            # eval'd items such as namedtuple give invalid items for their function __module__
            modname = '__main__'

        if modname == '__main__':
            themodule = None

        if themodule:
            self.modules.add(themodule)
            if getattr(themodule, name, None) is obj:
                return self.save_global(obj, name)

        # a builtin_function_or_method which comes in as an attribute of some
        # object (e.g., object.__new__, itertools.chain.from_iterable) will end
        # up with modname "__main__" and so end up here. But these functions
        # have no __code__ attribute in CPython, so the handling for
        # user-defined functions below will fail.
        # So we pickle them here using save_reduce; have to do it differently
        # for different python versions.
        if not hasattr(obj, '__code__'):
            if PY3:
                if sys.version_info < (3, 4):
                    raise pickle.PicklingError("Can't pickle %r" % obj)
                else:
                    rv = obj.__reduce_ex__(self.proto)
            else:
                if hasattr(obj, '__self__'):
                    rv = (getattr, (obj.__self__, name))
                else:
                    raise pickle.PicklingError("Can't pickle %r" % obj)
            return Pickler.save_reduce(self, obj=obj, *rv)

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
            klass = getattr(themodule, name, None)
            if klass is None or klass is not obj:
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
            if isinstance(x, types.ModuleType) and x.__package__:
                # check if the package has any currently loaded sub-imports
                prefix = x.__name__ + '.'
                for name, module in sys.modules.items():
                    if name.startswith(prefix):
                        # check whether the function can address the sub-module
                        tokens = set(name[len(prefix):].split('.'))
                        if not tokens - set(code.co_names):
                            # ensure unpickler executes this import
                            self.save(module)
                            # then discards the reference to it
                            self.write(pickle.POP)

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
        save(f_globals)
        save(defaults)
        save(dct)
        save(closure_values)
        write(pickle.TUPLE)
        write(pickle.REDUCE)  # applies _fill_function on the tuple

    _extract_code_globals_cache = (
        weakref.WeakKeyDictionary()
        if sys.version_info >= (2, 7) and not hasattr(sys, "pypy_version_info")
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
                out_names = set(names[oparg]
                                for op, oparg in _walk_global_ops(co))

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
            [c.cell_contents for c in func.__closure__]
            if func.__closure__ is not None else None
        )

        # save the dict
        dct = func.__dict__

        base_globals = self.globals_ref.get(id(func.__globals__), {})
        self.globals_ref[id(func.__globals__)] = base_globals

        return (code, f_globals, defaults, closure, dct, base_globals)

    def save_builtin_function(self, obj):
        if obj.__module__ == "__builtin__":
            return self.save_global(obj)
        return self.save_function(obj)
    dispatch[types.BuiltinFunctionType] = save_builtin_function

    def save_global(self, obj, name=None, pack=struct.pack):
        if obj.__module__ == "__builtin__" or obj.__module__ == "builtins":
            if obj in _BUILTIN_TYPE_NAMES:
                return self.save_reduce(_builtin_type, (_BUILTIN_TYPE_NAMES[obj],), obj=obj)

        if name is None:
            name = obj.__name__

        modname = getattr(obj, "__module__", None)
        if modname is None:
            modname = pickle.whichmodule(obj, name)

        if modname == '__main__':
            themodule = None
        else:
            __import__(modname)
            themodule = sys.modules[modname]
            self.modules.add(themodule)

        if hasattr(themodule, name) and getattr(themodule, name) is obj:
            return Pickler.save_global(self, obj, name)

        typ = type(obj)
        if typ is not obj and isinstance(obj, (type, types.ClassType)):
            d = dict(obj.__dict__)  # copy dict proxy to a dict
            if not isinstance(d.get('__dict__', None), property):
                # don't extract dict that are properties
                d.pop('__dict__', None)
            d.pop('__weakref__', None)

            # hack as __new__ is stored differently in the __dict__
            new_override = d.get('__new__', None)
            if new_override:
                d['__new__'] = obj.__new__

            self.save_reduce(typ, (obj.__name__, obj.__bases__, d), obj=obj)
        else:
            raise pickle.PicklingError("Can't pickle %r" % obj)

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
        """Inner logic to save instance. Based off pickle.save_inst
        Supports __transient__"""
        cls = obj.__class__

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
            #remove items if transient
            if hasattr(obj, '__transient__'):
                transient = obj.__transient__
                stuff = stuff.copy()
                for k in list(stuff.keys()):
                    if k in transient:
                        del stuff[k]
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
        try:
            orig_func = obj.__func__
        except AttributeError:  # Python 2.6
            orig_func = obj.__get__(None, object)
            if isinstance(obj, classmethod):
                orig_func = orig_func.__func__  # Unbind
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
            items = (items, )
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

    def save_reduce(self, func, args, state=None,
                    listitems=None, dictitems=None, obj=None):
        """Modified to support __transient__ on new objects
        Change only affects protocol level 2 (which is always used by PiCloud"""
        # Assert that args is a tuple or None
        if not isinstance(args, tuple):
            raise pickle.PicklingError("args from reduce() should be a tuple")

        # Assert that func is callable
        if not hasattr(func, '__call__'):
            raise pickle.PicklingError("func from reduce should be callable")

        save = self.save
        write = self.write

        # Protocol 2 special case: if func's name is __newobj__, use NEWOBJ
        if self.proto >= 2 and getattr(func, "__name__", "") == "__newobj__":
            #Added fix to allow transient
            cls = args[0]
            if not hasattr(cls, "__new__"):
                raise pickle.PicklingError(
                    "args[0] from __newobj__ args has no __new__")
            if obj is not None and cls is not obj.__class__:
                raise pickle.PicklingError(
                    "args[0] from __newobj__ args has the wrong class")
            args = args[1:]
            save(cls)

            #Don't pickle transient entries
            if hasattr(obj, '__transient__'):
                transient = obj.__transient__
                state = state.copy()

                for k in list(state.keys()):
                    if k in transient:
                        del state[k]

            save(args)
            write(pickle.NEWOBJ)
        else:
            save(func)
            save(args)
            write(pickle.REDUCE)

        if obj is not None:
            self.memoize(obj)

        # More new special cases (that work with older protocols as
        # well): when __reduce__ returns a tuple with 4 or 5 items,
        # the 4th and 5th item should be iterators that provide list
        # items and dict items (as (key, value) tuples), or None.

        if listitems is not None:
            self._batch_appends(listitems)

        if dictitems is not None:
            self._batch_setitems(dictitems)

        if state is not None:
            save(state)
            write(pickle.BUILD)

    def save_partial(self, obj):
        """Partial objects do not serialize correctly in python2.x -- this fixes the bugs"""
        self.save_reduce(_genpartial, (obj.func, obj.args, obj.keywords))

    if sys.version_info < (2,7):  # 2.7 supports partial pickling
        dispatch[partial] = save_partial


    def save_file(self, obj):
        """Save a file"""
        try:
            import StringIO as pystringIO #we can't use cStringIO as it lacks the name attribute
        except ImportError:
            import io as pystringIO

        if not hasattr(obj, 'name') or  not hasattr(obj, 'mode'):
            raise pickle.PicklingError("Cannot pickle files that do not map to an actual file")
        if obj is sys.stdout:
            return self.save_reduce(getattr, (sys,'stdout'), obj=obj)
        if obj is sys.stderr:
            return self.save_reduce(getattr, (sys,'stderr'), obj=obj)
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

    if PY3:
        dispatch[io.TextIOWrapper] = save_file
    else:
        dispatch[file] = save_file

    dispatch[type(Ellipsis)] = save_ellipsis
    dispatch[type(NotImplemented)] = save_not_implemented

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

def dump(obj, file, protocol=2):
    CloudPickler(file, protocol).dump(obj)


def dumps(obj, protocol=2):
    file = StringIO()

    cp = CloudPickler(file,protocol)
    cp.dump(obj)

    return file.getvalue()

# including pickles unloading functions in this namespace
load = pickle.load
loads = pickle.loads


#hack for __import__ not working as desired
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
            except Exception as e:
                sys.stderr.write('warning: could not import %s\n.  '
                                 'Your function may unexpectedly error due to this import failing;'
                                 'A version mismatch is likely.  Specific error was:\n' % modname)
                print_exec(sys.stderr)
            else:
                setattr(main, mod.__name__, mod)


#object generators:
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

def _fill_function(func, globals, defaults, dict, closure_values):
    """ Fills in the rest of function data into the skeleton function object
        that were created via _make_skel_func().
    """
    func.__globals__.update(globals)
    func.__defaults__ = defaults
    func.__dict__ = dict

    cells = func.__closure__
    if cells is not None:
        for cell, value in zip(cells, closure_values):
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
    base_globals['__builtins__'] = __builtins__

    closure = (
        tuple(_make_empty_cell() for _ in range(cell_count))
        if cell_count >= 0 else
        None
    )
    return types.FunctionType(code, base_globals, None, None, closure)


def _find_module(mod_name):
    """
    Iterate over each part instead of calling imp.find_module directly.
    This function is able to find submodules (e.g. sickit.tree)
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
