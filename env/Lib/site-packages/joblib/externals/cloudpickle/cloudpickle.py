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

import abc
import builtins
import dis
import io
import itertools
import logging
import opcode
import operator
import pickle
import platform
import struct
import sys
import types
import weakref
import uuid
import threading
import typing
from enum import Enum

from typing import Generic, Union, Tuple, Callable
from pickle import _Pickler as Pickler
from pickle import _getattribute
from io import BytesIO
from importlib._bootstrap import _find_spec

try:  # pragma: no branch
    import typing_extensions as _typing_extensions
    from typing_extensions import Literal, Final
except ImportError:
    _typing_extensions = Literal = Final = None

if sys.version_info >= (3, 5, 3):
    from typing import ClassVar
else:  # pragma: no cover
    ClassVar = None


# cloudpickle is meant for inter process communication: we expect all
# communicating processes to run the same Python version hence we favor
# communication speed over compatibility:
DEFAULT_PROTOCOL = pickle.HIGHEST_PROTOCOL

# Track the provenance of reconstructed dynamic classes to make it possible to
# recontruct instances from the matching singleton class definition when
# appropriate and preserve the usual "isinstance" semantics of Python objects.
_DYNAMIC_CLASS_TRACKER_BY_CLASS = weakref.WeakKeyDictionary()
_DYNAMIC_CLASS_TRACKER_BY_ID = weakref.WeakValueDictionary()
_DYNAMIC_CLASS_TRACKER_LOCK = threading.Lock()

PYPY = platform.python_implementation() == "PyPy"

builtin_code_type = None
if PYPY:
    # builtin-code objects only exist in pypy
    builtin_code_type = type(float.__new__.__code__)

_extract_code_globals_cache = weakref.WeakKeyDictionary()


def _get_or_create_tracker_id(class_def):
    with _DYNAMIC_CLASS_TRACKER_LOCK:
        class_tracker_id = _DYNAMIC_CLASS_TRACKER_BY_CLASS.get(class_def)
        if class_tracker_id is None:
            class_tracker_id = uuid.uuid4().hex
            _DYNAMIC_CLASS_TRACKER_BY_CLASS[class_def] = class_tracker_id
            _DYNAMIC_CLASS_TRACKER_BY_ID[class_tracker_id] = class_def
    return class_tracker_id


def _lookup_class_or_track(class_tracker_id, class_def):
    if class_tracker_id is not None:
        with _DYNAMIC_CLASS_TRACKER_LOCK:
            class_def = _DYNAMIC_CLASS_TRACKER_BY_ID.setdefault(
                class_tracker_id, class_def)
            _DYNAMIC_CLASS_TRACKER_BY_CLASS[class_def] = class_tracker_id
    return class_def


def _whichmodule(obj, name):
    """Find the module an object belongs to.

    This function differs from ``pickle.whichmodule`` in two ways:
    - it does not mangle the cases where obj's module is __main__ and obj was
      not found in any module.
    - Errors arising during module introspection are ignored, as those errors
      are considered unwanted side effects.
    """
    if sys.version_info[:2] < (3, 7) and isinstance(obj, typing.TypeVar):  # pragma: no branch  # noqa
        # Workaround bug in old Python versions: prior to Python 3.7,
        # T.__module__ would always be set to "typing" even when the TypeVar T
        # would be defined in a different module.
        #
        # For such older Python versions, we ignore the __module__ attribute of
        # TypeVar instances and instead exhaustively lookup those instances in
        # all currently imported modules.
        module_name = None
    else:
        module_name = getattr(obj, '__module__', None)

    if module_name is not None:
        return module_name
    # Protect the iteration by using a copy of sys.modules against dynamic
    # modules that trigger imports of other modules upon calls to getattr or
    # other threads importing at the same time.
    for module_name, module in sys.modules.copy().items():
        # Some modules such as coverage can inject non-module objects inside
        # sys.modules
        if (
                module_name == '__main__' or
                module is None or
                not isinstance(module, types.ModuleType)
        ):
            continue
        try:
            if _getattribute(module, name)[0] is obj:
                return module_name
        except Exception:
            pass
    return None


def _is_importable_by_name(obj, name=None):
    """Determine if obj can be pickled as attribute of a file-backed module"""
    return _lookup_module_and_qualname(obj, name=name) is not None


def _lookup_module_and_qualname(obj, name=None):
    if name is None:
        name = getattr(obj, '__qualname__', None)
    if name is None:  # pragma: no cover
        # This used to be needed for Python 2.7 support but is probably not
        # needed anymore. However we keep the __name__ introspection in case
        # users of cloudpickle rely on this old behavior for unknown reasons.
        name = getattr(obj, '__name__', None)

    module_name = _whichmodule(obj, name)

    if module_name is None:
        # In this case, obj.__module__ is None AND obj was not found in any
        # imported module. obj is thus treated as dynamic.
        return None

    if module_name == "__main__":
        return None

    module = sys.modules.get(module_name, None)
    if module is None:
        # The main reason why obj's module would not be imported is that this
        # module has been dynamically created, using for example
        # types.ModuleType. The other possibility is that module was removed
        # from sys.modules after obj was created/imported. But this case is not
        # supported, as the standard pickle does not support it either.
        return None

    # module has been added to sys.modules, but it can still be dynamic.
    if _is_dynamic(module):
        return None

    try:
        obj2, parent = _getattribute(module, name)
    except AttributeError:
        # obj was not found inside the module it points to
        return None
    if obj2 is not obj:
        return None
    return module, name


def _extract_code_globals(co):
    """
    Find all globals names read or written to by codeblock co
    """
    out_names = _extract_code_globals_cache.get(co)
    if out_names is None:
        names = co.co_names
        out_names = {names[oparg] for _, oparg in _walk_global_ops(co)}

        # Declaring a function inside another one using the "def ..."
        # syntax generates a constant code object corresonding to the one
        # of the nested function's As the nested function may itself need
        # global variables, we need to introspect its code, extract its
        # globals, (look for code object in it's co_consts attribute..) and
        # add the result to code_globals
        if co.co_consts:
            for const in co.co_consts:
                if isinstance(const, types.CodeType):
                    out_names |= _extract_code_globals(const)

        _extract_code_globals_cache[co] = out_names

    return out_names


def _find_imported_submodules(code, top_level_dependencies):
    """
    Find currently imported submodules used by a function.

    Submodules used by a function need to be detected and referenced for the
    function to work correctly at depickling time. Because submodules can be
    referenced as attribute of their parent package (``package.submodule``), we
    need a special introspection technique that does not rely on GLOBAL-related
    opcodes to find references of them in a code object.

    Example:
    ```
    import concurrent.futures
    import cloudpickle
    def func():
        x = concurrent.futures.ThreadPoolExecutor
    if __name__ == '__main__':
        cloudpickle.dumps(func)
    ```
    The globals extracted by cloudpickle in the function's state include the
    concurrent package, but not its submodule (here, concurrent.futures), which
    is the module used by func. Find_imported_submodules will detect the usage
    of concurrent.futures. Saving this module alongside with func will ensure
    that calling func once depickled does not fail due to concurrent.futures
    not being imported
    """

    subimports = []
    # check if any known dependency is an imported package
    for x in top_level_dependencies:
        if (isinstance(x, types.ModuleType) and
                hasattr(x, '__package__') and x.__package__):
            # check if the package has any currently loaded sub-imports
            prefix = x.__name__ + '.'
            # A concurrent thread could mutate sys.modules,
            # make sure we iterate over a copy to avoid exceptions
            for name in list(sys.modules):
                # Older versions of pytest will add a "None" module to
                # sys.modules.
                if name is not None and name.startswith(prefix):
                    # check whether the function can address the sub-module
                    tokens = set(name[len(prefix):].split('.'))
                    if not tokens - set(code.co_names):
                        subimports.append(sys.modules[name])
    return subimports


def cell_set(cell, value):
    """Set the value of a closure cell.

    The point of this function is to set the cell_contents attribute of a cell
    after its creation. This operation is necessary in case the cell contains a
    reference to the function the cell belongs to, as when calling the
    function's constructor
    ``f = types.FunctionType(code, globals, name, argdefs, closure)``,
    closure will not be able to contain the yet-to-be-created f.

    In Python3.7, cell_contents is writeable, so setting the contents of a cell
    can be done simply using
    >>> cell.cell_contents = value

    In earlier Python3 versions, the cell_contents attribute of a cell is read
    only, but this limitation can be worked around by leveraging the Python 3
    ``nonlocal`` keyword.

    In Python2 however, this attribute is read only, and there is no
    ``nonlocal`` keyword. For this reason, we need to come up with more
    complicated hacks to set this attribute.

    The chosen approach is to create a function with a STORE_DEREF opcode,
    which sets the content of a closure variable. Typically:

    >>> def inner(value):
    ...     lambda: cell  # the lambda makes cell a closure
    ...     cell = value  # cell is a closure, so this triggers a STORE_DEREF

    (Note that in Python2, A STORE_DEREF can never be triggered from an inner
    function. The function g for example here
    >>> def f(var):
    ...     def g():
    ...         var += 1
    ...     return g

    will not modify the closure variable ``var```inplace, but instead try to
    load a local variable var and increment it. As g does not assign the local
    variable ``var`` any initial value, calling f(1)() will fail at runtime.)

    Our objective is to set the value of a given cell ``cell``. So we need to
    somewhat reference our ``cell`` object into the ``inner`` function so that
    this object (and not the smoke cell of the lambda function) gets affected
    by the STORE_DEREF operation.

    In inner, ``cell`` is referenced as a cell variable (an enclosing variable
    that is referenced by the inner function). If we create a new function
    cell_set with the exact same code as ``inner``, but with ``cell`` marked as
    a free variable instead, the STORE_DEREF will be applied on its closure -
    ``cell``, which we can specify explicitly during construction! The new
    cell_set variable thus actually sets the contents of a specified cell!

    Note: we do not make use of the ``nonlocal`` keyword to set the contents of
    a cell in early python3 versions to limit possible syntax errors in case
    test and checker libraries decide to parse the whole file.
    """

    if sys.version_info[:2] >= (3, 7):  # pragma: no branch
        cell.cell_contents = value
    else:
        _cell_set = types.FunctionType(
            _cell_set_template_code, {}, '_cell_set', (), (cell,),)
        _cell_set(value)


def _make_cell_set_template_code():
    def _cell_set_factory(value):
        lambda: cell
        cell = value

    co = _cell_set_factory.__code__

    _cell_set_template_code = types.CodeType(
        co.co_argcount,
        co.co_kwonlyargcount,   # Python 3 only argument
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
        co.co_cellvars,  # co_freevars is initialized with co_cellvars
        (),  # co_cellvars is made empty
    )
    return _cell_set_template_code


if sys.version_info[:2] < (3, 7):
    _cell_set_template_code = _make_cell_set_template_code()

# relevant opcodes
STORE_GLOBAL = opcode.opmap['STORE_GLOBAL']
DELETE_GLOBAL = opcode.opmap['DELETE_GLOBAL']
LOAD_GLOBAL = opcode.opmap['LOAD_GLOBAL']
GLOBAL_OPS = (STORE_GLOBAL, DELETE_GLOBAL, LOAD_GLOBAL)
HAVE_ARGUMENT = dis.HAVE_ARGUMENT
EXTENDED_ARG = dis.EXTENDED_ARG


_BUILTIN_TYPE_NAMES = {}
for k, v in types.__dict__.items():
    if type(v) is type:
        _BUILTIN_TYPE_NAMES[v] = k


def _builtin_type(name):
    if name == "ClassType":  # pragma: no cover
        # Backward compat to load pickle files generated with cloudpickle
        # < 1.3 even if loading pickle files from older versions is not
        # officially supported.
        return type
    return getattr(types, name)


def _walk_global_ops(code):
    """
    Yield (opcode, argument number) tuples for all
    global-referencing instructions in *code*.
    """
    for instr in dis.get_instructions(code):
        op = instr.opcode
        if op in GLOBAL_OPS:
            yield op, instr.arg


def _extract_class_dict(cls):
    """Retrieve a copy of the dict of a class without the inherited methods"""
    clsdict = dict(cls.__dict__)  # copy dict proxy to a dict
    if len(cls.__bases__) == 1:
        inherited_dict = cls.__bases__[0].__dict__
    else:
        inherited_dict = {}
        for base in reversed(cls.__bases__):
            inherited_dict.update(base.__dict__)
    to_remove = []
    for name, value in clsdict.items():
        try:
            base_value = inherited_dict[name]
            if value is base_value:
                to_remove.append(name)
        except KeyError:
            pass
    for name in to_remove:
        clsdict.pop(name)
    return clsdict


if sys.version_info[:2] < (3, 7):  # pragma: no branch
    def _is_parametrized_type_hint(obj):
        # This is very cheap but might generate false positives.
        # general typing Constructs
        is_typing = getattr(obj, '__origin__', None) is not None

        # typing_extensions.Literal
        is_litteral = getattr(obj, '__values__', None) is not None

        # typing_extensions.Final
        is_final = getattr(obj, '__type__', None) is not None

        # typing.Union/Tuple for old Python 3.5
        is_union = getattr(obj, '__union_params__', None) is not None
        is_tuple = getattr(obj, '__tuple_params__', None) is not None
        is_callable = (
            getattr(obj, '__result__', None) is not None and
            getattr(obj, '__args__', None) is not None
        )
        return any((is_typing, is_litteral, is_final, is_union, is_tuple,
                    is_callable))

    def _create_parametrized_type_hint(origin, args):
        return origin[args]


class CloudPickler(Pickler):

    dispatch = Pickler.dispatch.copy()

    def __init__(self, file, protocol=None):
        if protocol is None:
            protocol = DEFAULT_PROTOCOL
        Pickler.__init__(self, file, protocol=protocol)
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

    def save_typevar(self, obj):
        self.save_reduce(*_typevar_reduce(obj), obj=obj)

    dispatch[typing.TypeVar] = save_typevar

    def save_memoryview(self, obj):
        self.save(obj.tobytes())

    dispatch[memoryview] = save_memoryview

    def save_module(self, obj):
        """
        Save a module as an import
        """
        if _is_dynamic(obj):
            obj.__dict__.pop('__builtins__', None)
            self.save_reduce(dynamic_subimport, (obj.__name__, vars(obj)),
                             obj=obj)
        else:
            self.save_reduce(subimport, (obj.__name__,), obj=obj)

    dispatch[types.ModuleType] = save_module

    def save_codeobject(self, obj):
        """
        Save a code object
        """
        if hasattr(obj, "co_posonlyargcount"):  # pragma: no branch
            args = (
                obj.co_argcount, obj.co_posonlyargcount,
                obj.co_kwonlyargcount, obj.co_nlocals, obj.co_stacksize,
                obj.co_flags, obj.co_code, obj.co_consts, obj.co_names,
                obj.co_varnames, obj.co_filename, obj.co_name,
                obj.co_firstlineno, obj.co_lnotab, obj.co_freevars,
                obj.co_cellvars
            )
        else:
            args = (
                obj.co_argcount, obj.co_kwonlyargcount, obj.co_nlocals,
                obj.co_stacksize, obj.co_flags, obj.co_code, obj.co_consts,
                obj.co_names, obj.co_varnames, obj.co_filename,
                obj.co_name, obj.co_firstlineno, obj.co_lnotab,
                obj.co_freevars, obj.co_cellvars
            )
        self.save_reduce(types.CodeType, args, obj=obj)

    dispatch[types.CodeType] = save_codeobject

    def save_function(self, obj, name=None):
        """ Registered with the dispatch to handle all function types.

        Determines what kind of function obj is (e.g. lambda, defined at
        interactive prompt, etc) and handles the pickling appropriately.
        """
        if _is_importable_by_name(obj, name=name):
            return Pickler.save_global(self, obj, name=name)
        elif PYPY and isinstance(obj.__code__, builtin_code_type):
            return self.save_pypy_builtin_func(obj)
        else:
            return self.save_function_tuple(obj)

    dispatch[types.FunctionType] = save_function

    def save_pypy_builtin_func(self, obj):
        """Save pypy equivalent of builtin functions.

        PyPy does not have the concept of builtin-functions. Instead,
        builtin-functions are simple function instances, but with a
        builtin-code attribute.
        Most of the time, builtin functions should be pickled by attribute. But
        PyPy has flaky support for __qualname__, so some builtin functions such
        as float.__new__ will be classified as dynamic. For this reason only,
        we created this special routine. Because builtin-functions are not
        expected to have closure or globals, there is no additional hack
        (compared the one already implemented in pickle) to protect ourselves
        from reference cycles. A simple (reconstructor, newargs, obj.__dict__)
        tuple is save_reduced.

        Note also that PyPy improved their support for __qualname__ in v3.6, so
        this routing should be removed when cloudpickle supports only PyPy 3.6
        and later.
        """
        rv = (types.FunctionType, (obj.__code__, {}, obj.__name__,
                                   obj.__defaults__, obj.__closure__),
              obj.__dict__)
        self.save_reduce(*rv, obj=obj)

    def _save_dynamic_enum(self, obj, clsdict):
        """Special handling for dynamic Enum subclasses

        Use a dedicated Enum constructor (inspired by EnumMeta.__call__) as the
        EnumMeta metaclass has complex initialization that makes the Enum
        subclasses hold references to their own instances.
        """
        members = dict((e.name, e.value) for e in obj)

        self.save_reduce(
                _make_skeleton_enum,
                (obj.__bases__, obj.__name__, obj.__qualname__,
                 members, obj.__module__, _get_or_create_tracker_id(obj), None),
                obj=obj
         )

        # Cleanup the clsdict that will be passed to _rehydrate_skeleton_class:
        # Those attributes are already handled by the metaclass.
        for attrname in ["_generate_next_value_", "_member_names_",
                         "_member_map_", "_member_type_",
                         "_value2member_map_"]:
            clsdict.pop(attrname, None)
        for member in members:
            clsdict.pop(member)

    def save_dynamic_class(self, obj):
        """Save a class that can't be stored as module global.

        This method is used to serialize classes that are defined inside
        functions, or that otherwise can't be serialized as attribute lookups
        from global modules.
        """
        clsdict = _extract_class_dict(obj)
        clsdict.pop('__weakref__', None)

        if issubclass(type(obj), abc.ABCMeta):
            # If obj is an instance of an ABCMeta subclass, dont pickle the
            # cache/negative caches populated during isinstance/issubclass
            # checks, but pickle the list of registered subclasses of obj.
            clsdict.pop('_abc_cache', None)
            clsdict.pop('_abc_negative_cache', None)
            clsdict.pop('_abc_negative_cache_version', None)
            registry = clsdict.pop('_abc_registry', None)
            if registry is None:
                # in Python3.7+, the abc caches and registered subclasses of a
                # class are bundled into the single _abc_impl attribute
                clsdict.pop('_abc_impl', None)
                (registry, _, _, _) = abc._get_dump(obj)

                clsdict["_abc_impl"] = [subclass_weakref()
                                        for subclass_weakref in registry]
            else:
                # In the above if clause, registry is a set of weakrefs -- in
                # this case, registry is a WeakSet
                clsdict["_abc_impl"] = [type_ for type_ in registry]

        # On PyPy, __doc__ is a readonly attribute, so we need to include it in
        # the initial skeleton class.  This is safe because we know that the
        # doc can't participate in a cycle with the original class.
        type_kwargs = {'__doc__': clsdict.pop('__doc__', None)}

        if "__slots__" in clsdict:
            type_kwargs['__slots__'] = obj.__slots__
            # pickle string length optimization: member descriptors of obj are
            # created automatically from obj's __slots__ attribute, no need to
            # save them in obj's state
            if isinstance(obj.__slots__, str):
                clsdict.pop(obj.__slots__)
            else:
                for k in obj.__slots__:
                    clsdict.pop(k, None)

        # If type overrides __dict__ as a property, include it in the type
        # kwargs. In Python 2, we can't set this attribute after construction.
        # XXX: can this ever happen in Python 3? If so add a test.
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
        if Enum is not None and issubclass(obj, Enum):
            # Special handling of Enum subclasses
            self._save_dynamic_enum(obj, clsdict)
        else:
            # "Regular" class definition:
            tp = type(obj)
            self.save_reduce(_make_skeleton_class,
                             (tp, obj.__name__, _get_bases(obj), type_kwargs,
                              _get_or_create_tracker_id(obj), None),
                             obj=obj)

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

        # Extract currently-imported submodules used by func. Storing these
        # modules in a smoke _cloudpickle_subimports attribute of the object's
        # state will trigger the side effect of importing these modules at
        # unpickling time (which is necessary for func to work correctly once
        # depickled)
        submodules = _find_imported_submodules(
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
            '_cloudpickle_submodules': submodules
        }
        if hasattr(func, '__annotations__'):
            state['annotations'] = func.__annotations__
        if hasattr(func, '__qualname__'):
            state['qualname'] = func.__qualname__
        if hasattr(func, '__kwdefaults__'):
            state['kwdefaults'] = func.__kwdefaults__
        save(state)
        write(pickle.TUPLE)
        write(pickle.REDUCE)  # applies _fill_function on the tuple

    def extract_func_data(self, func):
        """
        Turn the function into a tuple of data necessary to recreate it:
            code, globals, defaults, closure_values, dict
        """
        code = func.__code__

        # extract all global ref's
        func_global_refs = _extract_code_globals(code)

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

        # base_globals represents the future global namespace of func at
        # unpickling time. Looking it up and storing it in globals_ref allow
        # functions sharing the same globals at pickling time to also
        # share them once unpickled, at one condition: since globals_ref is
        # an attribute of a Cloudpickler instance, and that a new CloudPickler is
        # created each time pickle.dump or pickle.dumps is called, functions
        # also need to be saved within the same invokation of
        # cloudpickle.dump/cloudpickle.dumps (for example: cloudpickle.dumps([f1, f2])). There
        # is no such limitation when using Cloudpickler.dump, as long as the
        # multiple invokations are bound to the same Cloudpickler.
        base_globals = self.globals_ref.setdefault(id(func.__globals__), {})

        if base_globals == {}:
            # Add module attributes used to resolve relative imports
            # instructions inside func.
            for k in ["__package__", "__name__", "__path__", "__file__"]:
                # Some built-in functions/methods such as object.__new__  have
                # their __globals__ set to None in PyPy
                if func.__globals__ is not None and k in func.__globals__:
                    base_globals[k] = func.__globals__[k]

        return (code, f_globals, defaults, closure, dct, base_globals)

    def save_getset_descriptor(self, obj):
        return self.save_reduce(getattr, (obj.__objclass__, obj.__name__))

    dispatch[types.GetSetDescriptorType] = save_getset_descriptor

    def save_global(self, obj, name=None, pack=struct.pack):
        """
        Save a "global".

        The name of this method is somewhat misleading: all types get
        dispatched here.
        """
        if obj is type(None):
            return self.save_reduce(type, (None,), obj=obj)
        elif obj is type(Ellipsis):
            return self.save_reduce(type, (Ellipsis,), obj=obj)
        elif obj is type(NotImplemented):
            return self.save_reduce(type, (NotImplemented,), obj=obj)
        elif obj in _BUILTIN_TYPE_NAMES:
            return self.save_reduce(
                _builtin_type, (_BUILTIN_TYPE_NAMES[obj],), obj=obj)

        if sys.version_info[:2] < (3, 7) and _is_parametrized_type_hint(obj):  # noqa  # pragma: no branch
            # Parametrized typing constructs in Python < 3.7 are not compatible
            # with type checks and ``isinstance`` semantics. For this reason,
            # it is easier to detect them using a duck-typing-based check
            # (``_is_parametrized_type_hint``) than to populate the Pickler's
            # dispatch with type-specific savers.
            self._save_parametrized_type_hint(obj)
        elif name is not None:
            Pickler.save_global(self, obj, name=name)
        elif not _is_importable_by_name(obj, name=name):
            self.save_dynamic_class(obj)
        else:
            Pickler.save_global(self, obj, name=name)

    dispatch[type] = save_global

    def save_instancemethod(self, obj):
        # Memoization rarely is ever useful due to python bounding
        if obj.__self__ is None:
            self.save_reduce(getattr, (obj.im_class, obj.__name__))
        else:
            self.save_reduce(types.MethodType, (obj.__func__, obj.__self__), obj=obj)

    dispatch[types.MethodType] = save_instancemethod

    def save_property(self, obj):
        # properties not correctly saved in python
        self.save_reduce(property, (obj.fget, obj.fset, obj.fdel, obj.__doc__),
                         obj=obj)

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

        # TODO: also support binary mode files with io.BytesIO
        retval = io.StringIO()

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

    if hasattr(types, "MappingProxyType"):  # pragma: no branch
        def save_mappingproxy(self, obj):
            self.save_reduce(types.MappingProxyType, (dict(obj),), obj=obj)

        dispatch[types.MappingProxyType] = save_mappingproxy

    """Special functions for Add-on libraries"""
    def inject_addons(self):
        """Plug in system. Register additional pickling functions if modules already loaded"""
        pass

    if sys.version_info < (3, 7):  # pragma: no branch
        def _save_parametrized_type_hint(self, obj):
            # The distorted type check sematic for typing construct becomes:
            # ``type(obj) is type(TypeHint)``, which means "obj is a
            # parametrized TypeHint"
            if type(obj) is type(Literal):  # pragma: no branch
                initargs = (Literal, obj.__values__)
            elif type(obj) is type(Final):  # pragma: no branch
                initargs = (Final, obj.__type__)
            elif type(obj) is type(ClassVar):
                initargs = (ClassVar, obj.__type__)
            elif type(obj) is type(Generic):
                parameters = obj.__parameters__
                if len(obj.__parameters__) > 0:
                    # in early Python 3.5, __parameters__ was sometimes
                    # preferred to __args__
                    initargs = (obj.__origin__, parameters)
                else:
                    initargs = (obj.__origin__, obj.__args__)
            elif type(obj) is type(Union):
                if sys.version_info < (3, 5, 3):  # pragma: no cover
                    initargs = (Union, obj.__union_params__)
                else:
                    initargs = (Union, obj.__args__)
            elif type(obj) is type(Tuple):
                if sys.version_info < (3, 5, 3):  # pragma: no cover
                    initargs = (Tuple, obj.__tuple_params__)
                else:
                    initargs = (Tuple, obj.__args__)
            elif type(obj) is type(Callable):
                if sys.version_info < (3, 5, 3):  # pragma: no cover
                    args = obj.__args__
                    result = obj.__result__
                    if args != Ellipsis:
                        if isinstance(args, tuple):
                            args = list(args)
                        else:
                            args = [args]
                else:
                    (*args, result) = obj.__args__
                    if len(args) == 1 and args[0] is Ellipsis:
                        args = Ellipsis
                    else:
                        args = list(args)
                initargs = (Callable, (args, result))
            else:  # pragma: no cover
                raise pickle.PicklingError(
                    "Cloudpickle Error: Unknown type {}".format(type(obj))
                )
            self.save_reduce(_create_parametrized_type_hint, initargs, obj=obj)


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
    file = BytesIO()
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
    mod = types.ModuleType(name)
    mod.__dict__.update(vars)
    mod.__dict__['__builtins__'] = builtins.__dict__
    return mod


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

    # - At pickling time, any dynamic global variable used by func is
    #   serialized by value (in state['globals']).
    # - At unpickling time, func's __globals__ attribute is initialized by
    #   first retrieving an empty isolated namespace that will be shared
    #   with other functions pickled from the same original module
    #   by the same CloudPickler instance and then updated with the
    #   content of state['globals'] to populate the shared isolated
    #   namespace with all the global variables that are specifically
    #   referenced for this function.
    func.__globals__.update(state['globals'])

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
    if 'kwdefaults' in state:
        func.__kwdefaults__ = state['kwdefaults']
    # _cloudpickle_subimports is a set of submodules that must be loaded for
    # the pickled function to work correctly at unpickling time. Now that these
    # submodules are depickled (hence imported), they can be removed from the
    # object's state (the object state only served as a reference holder to
    # these submodules)
    if '_cloudpickle_submodules' in state:
        state.pop('_cloudpickle_submodules')

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
    # This is backward-compatibility code: for cloudpickle versions between
    # 0.5.4 and 0.7, base_globals could be a string or None. base_globals
    # should now always be a dictionary.
    if base_globals is None or isinstance(base_globals, str):
        base_globals = {}

    base_globals['__builtins__'] = __builtins__

    closure = (
        tuple(_make_empty_cell() for _ in range(cell_count))
        if cell_count >= 0 else
        None
    )
    return types.FunctionType(code, base_globals, None, None, closure)


def _make_skeleton_class(type_constructor, name, bases, type_kwargs,
                         class_tracker_id, extra):
    """Build dynamic class with an empty __dict__ to be filled once memoized

    If class_tracker_id is not None, try to lookup an existing class definition
    matching that id. If none is found, track a newly reconstructed class
    definition under that id so that other instances stemming from the same
    class id will also reuse this class definition.

    The "extra" variable is meant to be a dict (or None) that can be used for
    forward compatibility shall the need arise.
    """
    skeleton_class = types.new_class(
        name, bases, {'metaclass': type_constructor},
        lambda ns: ns.update(type_kwargs)
    )
    return _lookup_class_or_track(class_tracker_id, skeleton_class)


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


def _make_skeleton_enum(bases, name, qualname, members, module,
                        class_tracker_id, extra):
    """Build dynamic enum with an empty __dict__ to be filled once memoized

    The creation of the enum class is inspired by the code of
    EnumMeta._create_.

    If class_tracker_id is not None, try to lookup an existing enum definition
    matching that id. If none is found, track a newly reconstructed enum
    definition under that id so that other instances stemming from the same
    class id will also reuse this enum definition.

    The "extra" variable is meant to be a dict (or None) that can be used for
    forward compatibility shall the need arise.
    """
    # enums always inherit from their base Enum class at the last position in
    # the list of base classes:
    enum_base = bases[-1]
    metacls = enum_base.__class__
    classdict = metacls.__prepare__(name, bases)

    for member_name, member_value in members.items():
        classdict[member_name] = member_value
    enum_class = metacls.__new__(metacls, name, bases, classdict)
    enum_class.__module__ = module
    enum_class.__qualname__ = qualname

    return _lookup_class_or_track(class_tracker_id, enum_class)


def _is_dynamic(module):
    """
    Return True if the module is special module that cannot be imported by its
    name.
    """
    # Quick check: module that have __file__ attribute are not dynamic modules.
    if hasattr(module, '__file__'):
        return False

    if module.__spec__ is not None:
        return False

    # In PyPy, Some built-in modules such as _codecs can have their
    # __spec__ attribute set to None despite being imported.  For such
    # modules, the ``_find_spec`` utility of the standard library is used.
    parent_name = module.__name__.rpartition('.')[0]
    if parent_name:  # pragma: no cover
        # This code handles the case where an imported package (and not
        # module) remains with __spec__ set to None. It is however untested
        # as no package in the PyPy stdlib has __spec__ set to None after
        # it is imported.
        try:
            parent = sys.modules[parent_name]
        except KeyError:
            msg = "parent {!r} not in sys.modules"
            raise ImportError(msg.format(parent_name))
        else:
            pkgpath = parent.__path__
    else:
        pkgpath = None
    return _find_spec(module.__name__, pkgpath, module) is None


def _make_typevar(name, bound, constraints, covariant, contravariant,
                  class_tracker_id):
    tv = typing.TypeVar(
        name, *constraints, bound=bound,
        covariant=covariant, contravariant=contravariant
    )
    if class_tracker_id is not None:
        return _lookup_class_or_track(class_tracker_id, tv)
    else:  # pragma: nocover
        # Only for Python 3.5.3 compat.
        return tv


def _decompose_typevar(obj):
    try:
        class_tracker_id = _get_or_create_tracker_id(obj)
    except TypeError:  # pragma: nocover
        # TypeVar instances are not weakref-able in Python 3.5.3
        class_tracker_id = None
    return (
        obj.__name__, obj.__bound__, obj.__constraints__,
        obj.__covariant__, obj.__contravariant__,
        class_tracker_id,
    )


def _typevar_reduce(obj):
    # TypeVar instances have no __qualname__ hence we pass the name explicitly.
    module_and_name = _lookup_module_and_qualname(obj, name=obj.__name__)
    if module_and_name is None:
        return (_make_typevar, _decompose_typevar(obj))
    return (getattr, module_and_name)


def _get_bases(typ):
    if hasattr(typ, '__orig_bases__'):
        # For generic types (see PEP 560)
        bases_attr = '__orig_bases__'
    else:
        # For regular class objects
        bases_attr = '__bases__'
    return getattr(typ, bases_attr)
