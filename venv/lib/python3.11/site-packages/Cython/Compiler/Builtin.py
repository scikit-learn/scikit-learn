#
#   Builtin Definitions
#

from __future__ import absolute_import

from .StringEncoding import EncodedString
from .Symtab import BuiltinScope, StructOrUnionScope, ModuleScope, Entry
from .Code import UtilityCode, TempitaUtilityCode
from .TypeSlots import Signature
from . import PyrexTypes


# C-level implementations of builtin types, functions and methods

iter_next_utility_code = UtilityCode.load("IterNext", "ObjectHandling.c")
getattr_utility_code = UtilityCode.load("GetAttr", "ObjectHandling.c")
getattr3_utility_code = UtilityCode.load("GetAttr3", "Builtins.c")
pyexec_utility_code = UtilityCode.load("PyExec", "Builtins.c")
pyexec_globals_utility_code = UtilityCode.load("PyExecGlobals", "Builtins.c")
globals_utility_code = UtilityCode.load("Globals", "Builtins.c")

builtin_utility_code = {
    'StopAsyncIteration': UtilityCode.load_cached("StopAsyncIteration", "Coroutine.c"),
}


# mapping from builtins to their C-level equivalents

class _BuiltinOverride(object):
    def __init__(self, py_name, args, ret_type, cname, py_equiv="*",
                 utility_code=None, sig=None, func_type=None,
                 is_strict_signature=False, builtin_return_type=None,
                 nogil=None):
        self.py_name, self.cname, self.py_equiv = py_name, cname, py_equiv
        self.args, self.ret_type = args, ret_type
        self.func_type, self.sig = func_type, sig
        self.builtin_return_type = builtin_return_type
        self.is_strict_signature = is_strict_signature
        self.utility_code = utility_code
        self.nogil = nogil

    def build_func_type(self, sig=None, self_arg=None):
        if sig is None:
            sig = Signature(self.args, self.ret_type, nogil=self.nogil)
            sig.exception_check = False  # not needed for the current builtins
        func_type = sig.function_type(self_arg)
        if self.is_strict_signature:
            func_type.is_strict_signature = True
        if self.builtin_return_type:
            func_type.return_type = builtin_types[self.builtin_return_type]
        return func_type


class BuiltinAttribute(object):
    def __init__(self, py_name, cname=None, field_type=None, field_type_name=None):
        self.py_name = py_name
        self.cname = cname or py_name
        self.field_type_name = field_type_name  # can't do the lookup before the type is declared!
        self.field_type = field_type

    def declare_in_type(self, self_type):
        if self.field_type_name is not None:
            # lazy type lookup
            field_type = builtin_scope.lookup(self.field_type_name).type
        else:
            field_type = self.field_type or PyrexTypes.py_object_type
        entry = self_type.scope.declare(self.py_name, self.cname, field_type, None, 'private')
        entry.is_variable = True


class BuiltinFunction(_BuiltinOverride):
    def declare_in_scope(self, scope):
        func_type, sig = self.func_type, self.sig
        if func_type is None:
            func_type = self.build_func_type(sig)
        scope.declare_builtin_cfunction(self.py_name, func_type, self.cname,
                                        self.py_equiv, self.utility_code)


class BuiltinMethod(_BuiltinOverride):
    def declare_in_type(self, self_type):
        method_type, sig = self.func_type, self.sig
        if method_type is None:
            # override 'self' type (first argument)
            self_arg = PyrexTypes.CFuncTypeArg("", self_type, None)
            self_arg.not_none = True
            self_arg.accept_builtin_subtypes = True
            method_type = self.build_func_type(sig, self_arg)
        self_type.scope.declare_builtin_cfunction(
            self.py_name, method_type, self.cname, utility_code=self.utility_code)


class BuiltinProperty(object):
    # read only for now
    def __init__(self, py_name, property_type, call_cname,
                 exception_value=None, exception_check=None, utility_code=None):
        self.py_name = py_name
        self.property_type = property_type
        self.call_cname = call_cname
        self.utility_code = utility_code
        self.exception_value = exception_value
        self.exception_check = exception_check

    def declare_in_type(self, self_type):
        self_type.scope.declare_cproperty(
            self.py_name,
            self.property_type,
            self.call_cname,
            exception_value=self.exception_value,
            exception_check=self.exception_check,
            utility_code=self.utility_code
        )


builtin_function_table = [
    # name,        args,   return,  C API func,           py equiv = "*"
    BuiltinFunction('abs',        "d",    "d",     "fabs",
                    is_strict_signature=True, nogil=True),
    BuiltinFunction('abs',        "f",    "f",     "fabsf",
                    is_strict_signature=True, nogil=True),
    BuiltinFunction('abs',        "i",    "i",     "abs",
                    is_strict_signature=True, nogil=True),
    BuiltinFunction('abs',        "l",    "l",     "labs",
                    is_strict_signature=True, nogil=True),
    BuiltinFunction('abs',        None,    None,   "__Pyx_abs_longlong",
                utility_code = UtilityCode.load("abs_longlong", "Builtins.c"),
                func_type = PyrexTypes.CFuncType(
                    PyrexTypes.c_longlong_type, [
                        PyrexTypes.CFuncTypeArg("arg", PyrexTypes.c_longlong_type, None)
                        ],
                    is_strict_signature = True, nogil=True)),
    ] + list(
        BuiltinFunction('abs',        None,    None,   "/*abs_{0}*/".format(t.specialization_name()),
                    func_type = PyrexTypes.CFuncType(
                        t,
                        [PyrexTypes.CFuncTypeArg("arg", t, None)],
                        is_strict_signature = True, nogil=True))
                            for t in (PyrexTypes.c_uint_type, PyrexTypes.c_ulong_type, PyrexTypes.c_ulonglong_type)
             ) + list(
        BuiltinFunction('abs',        None,    None,   "__Pyx_c_abs{0}".format(t.funcsuffix),
                    func_type = PyrexTypes.CFuncType(
                        t.real_type, [
                            PyrexTypes.CFuncTypeArg("arg", t, None)
                            ],
                            is_strict_signature = True, nogil=True))
                        for t in (PyrexTypes.c_float_complex_type,
                                  PyrexTypes.c_double_complex_type,
                                  PyrexTypes.c_longdouble_complex_type)
                        ) + [
    BuiltinFunction('abs',        "O",    "O",     "__Pyx_PyNumber_Absolute",
                    utility_code=UtilityCode.load("py_abs", "Builtins.c")),
    #('all',       "",     "",      ""),
    #('any',       "",     "",      ""),
    #('ascii',     "",     "",      ""),
    #('bin',       "",     "",      ""),
    BuiltinFunction('callable',   "O",    "b",     "__Pyx_PyCallable_Check",
                    utility_code = UtilityCode.load("CallableCheck", "ObjectHandling.c")),
    #('chr',       "",     "",      ""),
    #('cmp', "",   "",     "",      ""), # int PyObject_Cmp(PyObject *o1, PyObject *o2, int *result)
    #('compile',   "",     "",      ""), # PyObject* Py_CompileString(    char *str, char *filename, int start)
    BuiltinFunction('delattr',    "OO",   "r",     "PyObject_DelAttr"),
    BuiltinFunction('dir',        "O",    "O",     "PyObject_Dir"),
    BuiltinFunction('divmod',     "OO",   "O",     "PyNumber_Divmod"),
    BuiltinFunction('exec',       "O",    "O",     "__Pyx_PyExecGlobals",
                    utility_code = pyexec_globals_utility_code),
    BuiltinFunction('exec',       "OO",   "O",     "__Pyx_PyExec2",
                    utility_code = pyexec_utility_code),
    BuiltinFunction('exec',       "OOO",  "O",     "__Pyx_PyExec3",
                    utility_code = pyexec_utility_code),
    #('eval',      "",     "",      ""),
    #('execfile',  "",     "",      ""),
    #('filter',    "",     "",      ""),
    BuiltinFunction('getattr3',   "OOO",  "O",     "__Pyx_GetAttr3",     "getattr",
                    utility_code=getattr3_utility_code),  # Pyrex legacy
    BuiltinFunction('getattr',    "OOO",  "O",     "__Pyx_GetAttr3",
                    utility_code=getattr3_utility_code),
    BuiltinFunction('getattr',    "OO",   "O",     "__Pyx_GetAttr",
                    utility_code=getattr_utility_code),
    BuiltinFunction('hasattr',    "OO",   "b",     "__Pyx_HasAttr",
                    utility_code = UtilityCode.load("HasAttr", "Builtins.c")),
    BuiltinFunction('hash',       "O",    "h",     "PyObject_Hash"),
    #('hex',       "",     "",      ""),
    #('id',        "",     "",      ""),
    #('input',     "",     "",      ""),
    BuiltinFunction('intern',     "O",    "O",     "__Pyx_Intern",
                    utility_code = UtilityCode.load("Intern", "Builtins.c")),
    BuiltinFunction('isinstance', "OO",   "b",     "PyObject_IsInstance"),
    BuiltinFunction('issubclass', "OO",   "b",     "PyObject_IsSubclass"),
    BuiltinFunction('iter',       "OO",   "O",     "PyCallIter_New"),
    BuiltinFunction('iter',       "O",    "O",     "PyObject_GetIter"),
    BuiltinFunction('len',        "O",    "z",     "PyObject_Length"),
    BuiltinFunction('locals',     "",     "O",     "__pyx_locals"),
    #('map',       "",     "",      ""),
    #('max',       "",     "",      ""),
    #('min',       "",     "",      ""),
    BuiltinFunction('next',       "O",    "O",     "__Pyx_PyIter_Next",
                    utility_code = iter_next_utility_code),   # not available in Py2 => implemented here
    BuiltinFunction('next',      "OO",    "O",     "__Pyx_PyIter_Next2",
                    utility_code = iter_next_utility_code),  # not available in Py2 => implemented here
    #('oct',       "",     "",      ""),
    #('open',       "ss",   "O",     "PyFile_FromString"),   # not in Py3
] + [
    BuiltinFunction('ord',        None,    None,   "__Pyx_long_cast",
                    func_type=PyrexTypes.CFuncType(
                        PyrexTypes.c_long_type, [PyrexTypes.CFuncTypeArg("c", c_type, None)],
                        is_strict_signature=True))
    for c_type in [PyrexTypes.c_py_ucs4_type, PyrexTypes.c_py_unicode_type]
] + [
    BuiltinFunction('ord',        None,    None,   "__Pyx_uchar_cast",
                    func_type=PyrexTypes.CFuncType(
                        PyrexTypes.c_uchar_type, [PyrexTypes.CFuncTypeArg("c", c_type, None)],
                        is_strict_signature=True))
    for c_type in [PyrexTypes.c_char_type, PyrexTypes.c_schar_type, PyrexTypes.c_uchar_type]
] + [
    BuiltinFunction('ord',        None,    None,   "__Pyx_PyObject_Ord",
                    utility_code=UtilityCode.load_cached("object_ord", "Builtins.c"),
                    func_type=PyrexTypes.CFuncType(
                        PyrexTypes.c_long_type, [
                            PyrexTypes.CFuncTypeArg("c", PyrexTypes.py_object_type, None)
                        ],
                        exception_value="(long)(Py_UCS4)-1")),
    BuiltinFunction('pow',        "OOO",  "O",     "PyNumber_Power"),
    BuiltinFunction('pow',        "OO",   "O",     "__Pyx_PyNumber_Power2",
                    utility_code = UtilityCode.load("pow2", "Builtins.c")),
    #('range',     "",     "",      ""),
    #('raw_input', "",     "",      ""),
    #('reduce',    "",     "",      ""),
    BuiltinFunction('reload',     "O",    "O",     "PyImport_ReloadModule"),
    BuiltinFunction('repr',       "O",    "O",     "PyObject_Repr"),  # , builtin_return_type='str'),  # add in Cython 3.1
    #('round',     "",     "",      ""),
    BuiltinFunction('setattr',    "OOO",  "r",     "PyObject_SetAttr"),
    #('sum',       "",     "",      ""),
    #('sorted',    "",     "",      ""),
    #('type',       "O",    "O",     "PyObject_Type"),
    BuiltinFunction('unichr',     "i",    "O",      "PyUnicode_FromOrdinal", builtin_return_type='unicode'),
    #('unicode',   "",     "",      ""),
    #('vars',      "",     "",      ""),
    #('zip',       "",     "",      ""),
    #  Can't do these easily until we have builtin type entries.
    #('typecheck',  "OO",   "i",     "PyObject_TypeCheck", False),
    #('issubtype',  "OO",   "i",     "PyType_IsSubtype",   False),

    # Put in namespace append optimization.
    BuiltinFunction('__Pyx_PyObject_Append', "OO",  "O",     "__Pyx_PyObject_Append"),

    # This is conditionally looked up based on a compiler directive.
    BuiltinFunction('__Pyx_Globals',    "",     "O",     "__Pyx_Globals",
                    utility_code=globals_utility_code),
]


# Builtin types
#  bool
#  buffer
#  classmethod
#  dict
#  enumerate
#  file
#  float
#  int
#  list
#  long
#  object
#  property
#  slice
#  staticmethod
#  super
#  str
#  tuple
#  type
#  xrange

builtin_types_table = [

    ("type",    "PyType_Type",     []),

# This conflicts with the C++ bool type, and unfortunately
# C++ is too liberal about PyObject* <-> bool conversions,
# resulting in unintuitive runtime behavior and segfaults.
#    ("bool",    "PyBool_Type",     []),

    ("int",     "PyInt_Type",      []),
    ("long",    "PyLong_Type",     []),
    ("float",   "PyFloat_Type",    []),

    ("complex", "PyComplex_Type",  [BuiltinAttribute('cval', field_type_name = 'Py_complex'),
                                    BuiltinAttribute('real', 'cval.real', field_type = PyrexTypes.c_double_type),
                                    BuiltinAttribute('imag', 'cval.imag', field_type = PyrexTypes.c_double_type),
                                    ]),

    ("basestring", "PyBaseString_Type", [
                                    BuiltinMethod("join",  "TO",   "T", "__Pyx_PyBaseString_Join",
                                                  utility_code=UtilityCode.load("StringJoin", "StringTools.c")),
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),
    ("bytearray", "PyByteArray_Type", [
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),
    ("bytes",   "PyBytes_Type",    [BuiltinMethod("join",  "TO",   "O", "__Pyx_PyBytes_Join",
                                                  utility_code=UtilityCode.load("StringJoin", "StringTools.c")),
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),
    ("str",     "PyString_Type",   [BuiltinMethod("join",  "TO",   "O", "__Pyx_PyString_Join",
                                                  builtin_return_type='basestring',
                                                  utility_code=UtilityCode.load("StringJoin", "StringTools.c")),
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),
    ("unicode", "PyUnicode_Type",  [BuiltinMethod("__contains__",  "TO",   "b", "PyUnicode_Contains"),
                                    BuiltinMethod("join",  "TO",   "T", "PyUnicode_Join"),
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),

    ("tuple",   "PyTuple_Type",    [BuiltinMethod("__mul__", "Tz", "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),

    ("list",    "PyList_Type",     [BuiltinMethod("insert",  "TzO",  "r", "PyList_Insert"),
                                    BuiltinMethod("reverse", "T",    "r", "PyList_Reverse"),
                                    BuiltinMethod("append",  "TO",   "r", "__Pyx_PyList_Append",
                                                  utility_code=UtilityCode.load("ListAppend", "Optimize.c")),
                                    BuiltinMethod("extend",  "TO",   "r", "__Pyx_PyList_Extend",
                                                  utility_code=UtilityCode.load("ListExtend", "Optimize.c")),
                                    BuiltinMethod("__mul__",  "Tz",   "T", "__Pyx_PySequence_Multiply",
                                                  utility_code=UtilityCode.load("PySequenceMultiply", "ObjectHandling.c")),
                                    ]),

    ("dict",    "PyDict_Type",     [BuiltinMethod("__contains__",  "TO",   "b", "PyDict_Contains"),
                                    BuiltinMethod("has_key",       "TO",   "b", "PyDict_Contains"),
                                    BuiltinMethod("items",  "T",   "O", "__Pyx_PyDict_Items",
                                                  utility_code=UtilityCode.load("py_dict_items", "Builtins.c")),
                                    BuiltinMethod("keys",   "T",   "O", "__Pyx_PyDict_Keys",
                                                  utility_code=UtilityCode.load("py_dict_keys", "Builtins.c")),
                                    BuiltinMethod("values", "T",   "O", "__Pyx_PyDict_Values",
                                                  utility_code=UtilityCode.load("py_dict_values", "Builtins.c")),
                                    BuiltinMethod("iteritems",  "T",   "O", "__Pyx_PyDict_IterItems",
                                                  utility_code=UtilityCode.load("py_dict_iteritems", "Builtins.c")),
                                    BuiltinMethod("iterkeys",   "T",   "O", "__Pyx_PyDict_IterKeys",
                                                  utility_code=UtilityCode.load("py_dict_iterkeys", "Builtins.c")),
                                    BuiltinMethod("itervalues", "T",   "O", "__Pyx_PyDict_IterValues",
                                                  utility_code=UtilityCode.load("py_dict_itervalues", "Builtins.c")),
                                    BuiltinMethod("viewitems",  "T",   "O", "__Pyx_PyDict_ViewItems",
                                                  utility_code=UtilityCode.load("py_dict_viewitems", "Builtins.c")),
                                    BuiltinMethod("viewkeys",   "T",   "O", "__Pyx_PyDict_ViewKeys",
                                                  utility_code=UtilityCode.load("py_dict_viewkeys", "Builtins.c")),
                                    BuiltinMethod("viewvalues", "T",   "O", "__Pyx_PyDict_ViewValues",
                                                  utility_code=UtilityCode.load("py_dict_viewvalues", "Builtins.c")),
                                    BuiltinMethod("clear",  "T",   "r", "__Pyx_PyDict_Clear",
                                                  utility_code=UtilityCode.load("py_dict_clear", "Optimize.c")),
                                    BuiltinMethod("copy",   "T",   "T", "PyDict_Copy")]),

    ("slice",   "PySlice_Type",    [BuiltinAttribute('start'),
                                    BuiltinAttribute('stop'),
                                    BuiltinAttribute('step'),
                                    ]),
#    ("file",    "PyFile_Type",     []),  # not in Py3

    ("set",       "PySet_Type",    [BuiltinMethod("clear",   "T",  "r", "PySet_Clear"),
                                    # discard() and remove() have a special treatment for unhashable values
                                    BuiltinMethod("discard", "TO", "r", "__Pyx_PySet_Discard",
                                                  utility_code=UtilityCode.load("py_set_discard", "Optimize.c")),
                                    BuiltinMethod("remove",  "TO", "r", "__Pyx_PySet_Remove",
                                                  utility_code=UtilityCode.load("py_set_remove", "Optimize.c")),
                                    # update is actually variadic (see Github issue #1645)
#                                    BuiltinMethod("update",     "TO", "r", "__Pyx_PySet_Update",
#                                                  utility_code=UtilityCode.load_cached("PySet_Update", "Builtins.c")),
                                    BuiltinMethod("add",     "TO", "r", "PySet_Add"),
                                    BuiltinMethod("pop",     "T",  "O", "PySet_Pop")]),
    ("frozenset", "PyFrozenSet_Type", []),
    ("Exception", "((PyTypeObject*)PyExc_Exception)[0]", []),
    ("StopAsyncIteration", "((PyTypeObject*)__Pyx_PyExc_StopAsyncIteration)[0]", []),
    ("memoryview", "PyMemoryView_Type", [
        # TODO - format would be nice, but hard to get
        # __len__ can be accessed through a direct lookup of the buffer (but probably in Optimize.c)
        # error checking would ideally be limited api only
        BuiltinProperty("ndim", PyrexTypes.c_int_type, '__Pyx_PyMemoryView_Get_ndim',
                        exception_value="-1", exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="ndim")
                        )
        ),
        BuiltinProperty("readonly", PyrexTypes.c_bint_type, '__Pyx_PyMemoryView_Get_readonly',
                        exception_value="-1", exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="readonly")
                        )
        ),
        BuiltinProperty("itemsize", PyrexTypes.c_py_ssize_t_type, '__Pyx_PyMemoryView_Get_itemsize',
                        exception_value="-1", exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="itemsize")
                        )
        )]
    )
]


types_that_construct_their_instance = frozenset({
    # some builtin types do not always return an instance of
    # themselves - these do:
    'type', 'bool', 'long', 'float', 'complex',
    'bytes', 'unicode', 'bytearray',
    'tuple', 'list', 'dict', 'set', 'frozenset',
    # 'str',             # only in Py3.x
    # 'file',            # only in Py2.x
    'memoryview'
})


builtin_structs_table = [
    ('Py_buffer', 'Py_buffer',
     [("buf",        PyrexTypes.c_void_ptr_type),
      ("obj",        PyrexTypes.py_object_type),
      ("len",        PyrexTypes.c_py_ssize_t_type),
      ("itemsize",   PyrexTypes.c_py_ssize_t_type),
      ("readonly",   PyrexTypes.c_bint_type),
      ("ndim",       PyrexTypes.c_int_type),
      ("format",     PyrexTypes.c_char_ptr_type),
      ("shape",      PyrexTypes.c_py_ssize_t_ptr_type),
      ("strides",    PyrexTypes.c_py_ssize_t_ptr_type),
      ("suboffsets", PyrexTypes.c_py_ssize_t_ptr_type),
      ("smalltable", PyrexTypes.CArrayType(PyrexTypes.c_py_ssize_t_type, 2)),
      ("internal",   PyrexTypes.c_void_ptr_type),
      ]),
    ('Py_complex', 'Py_complex',
     [('real', PyrexTypes.c_double_type),
      ('imag', PyrexTypes.c_double_type),
      ])
]

# set up builtin scope

builtin_scope = BuiltinScope()

def init_builtin_funcs():
    for bf in builtin_function_table:
        bf.declare_in_scope(builtin_scope)

builtin_types = {}

def init_builtin_types():
    global builtin_types
    for name, cname, methods in builtin_types_table:
        utility = builtin_utility_code.get(name)
        if name == 'frozenset':
            objstruct_cname = 'PySetObject'
        elif name == 'bytearray':
            objstruct_cname = 'PyByteArrayObject'
        elif name == 'bool':
            objstruct_cname = None
        elif name == 'Exception':
            objstruct_cname = "PyBaseExceptionObject"
        elif name == 'StopAsyncIteration':
            objstruct_cname = "PyBaseExceptionObject"
        else:
            objstruct_cname = 'Py%sObject' % name.capitalize()
        type_class = PyrexTypes.BuiltinObjectType
        if name in ['dict', 'list', 'set', 'frozenset']:
            type_class = PyrexTypes.BuiltinTypeConstructorObjectType
        elif name == 'tuple':
            type_class = PyrexTypes.PythonTupleTypeConstructor
        the_type = builtin_scope.declare_builtin_type(name, cname, utility, objstruct_cname,
                                                      type_class=type_class)
        builtin_types[name] = the_type
        for method in methods:
            method.declare_in_type(the_type)

def init_builtin_structs():
    for name, cname, attribute_types in builtin_structs_table:
        scope = StructOrUnionScope(name)
        for attribute_name, attribute_type in attribute_types:
            scope.declare_var(attribute_name, attribute_type, None,
                              attribute_name, allow_pyobject=True)
        builtin_scope.declare_struct_or_union(
            name, "struct", scope, 1, None, cname = cname)


def init_builtins():
    #Errors.init_thread()  # hopefully not needed - we should not emit warnings ourselves
    init_builtin_structs()
    init_builtin_types()
    init_builtin_funcs()

    entry = builtin_scope.declare_var(
        '__debug__', PyrexTypes.c_const_type(PyrexTypes.c_bint_type),
        pos=None, cname='__pyx_assertions_enabled()', is_cdef=True)
    entry.utility_code = UtilityCode.load_cached("AssertionsEnabled", "Exceptions.c")

    global type_type, list_type, tuple_type, dict_type, set_type, frozenset_type, slice_type
    global bytes_type, str_type, unicode_type, basestring_type, bytearray_type
    global float_type, int_type, long_type, bool_type, complex_type
    global memoryview_type, py_buffer_type
    global sequence_types
    type_type  = builtin_scope.lookup('type').type
    list_type  = builtin_scope.lookup('list').type
    tuple_type = builtin_scope.lookup('tuple').type
    dict_type  = builtin_scope.lookup('dict').type
    set_type   = builtin_scope.lookup('set').type
    frozenset_type = builtin_scope.lookup('frozenset').type
    slice_type   = builtin_scope.lookup('slice').type

    bytes_type = builtin_scope.lookup('bytes').type
    str_type   = builtin_scope.lookup('str').type
    unicode_type = builtin_scope.lookup('unicode').type
    basestring_type = builtin_scope.lookup('basestring').type
    bytearray_type = builtin_scope.lookup('bytearray').type
    memoryview_type = builtin_scope.lookup('memoryview').type

    float_type = builtin_scope.lookup('float').type
    int_type = builtin_scope.lookup('int').type
    long_type = builtin_scope.lookup('long').type
    bool_type  = builtin_scope.lookup('bool').type
    complex_type  = builtin_scope.lookup('complex').type

    sequence_types = (
        list_type,
        tuple_type,
        bytes_type,
        str_type,
        unicode_type,
        basestring_type,
        bytearray_type,
        memoryview_type,
    )

    # Set up type inference links between equivalent Python/C types
    bool_type.equivalent_type = PyrexTypes.c_bint_type
    PyrexTypes.c_bint_type.equivalent_type = bool_type

    float_type.equivalent_type = PyrexTypes.c_double_type
    PyrexTypes.c_double_type.equivalent_type = float_type

    complex_type.equivalent_type = PyrexTypes.c_double_complex_type
    PyrexTypes.c_double_complex_type.equivalent_type = complex_type

    py_buffer_type = builtin_scope.lookup('Py_buffer').type


init_builtins()

##############################
# Support for a few standard library modules that Cython understands (currently typing and dataclasses)
##############################
_known_module_scopes = {}

def get_known_standard_library_module_scope(module_name):
    mod = _known_module_scopes.get(module_name)
    if mod:
        return mod

    if module_name == "typing":
        mod = ModuleScope(module_name, None, None)
        for name, tp in [
                ('Dict', dict_type),
                ('List', list_type),
                ('Tuple', tuple_type),
                ('Set', set_type),
                ('FrozenSet', frozenset_type),
                ]:
            name = EncodedString(name)
            entry = mod.declare_type(name, tp, pos = None)
            var_entry = Entry(name, None, PyrexTypes.py_object_type)
            var_entry.is_pyglobal = True
            var_entry.is_variable = True
            var_entry.scope = mod
            entry.as_variable = var_entry
            entry.known_standard_library_import = "%s.%s" % (module_name, name)

        for name in ['ClassVar', 'Optional']:
            name = EncodedString(name)
            indexed_type = PyrexTypes.SpecialPythonTypeConstructor(EncodedString("typing."+name))
            entry = mod.declare_type(name, indexed_type, pos = None)
            var_entry = Entry(name, None, PyrexTypes.py_object_type)
            var_entry.is_pyglobal = True
            var_entry.is_variable = True
            var_entry.scope = mod
            entry.as_variable = var_entry
            entry.known_standard_library_import = "%s.%s" % (module_name, name)
        _known_module_scopes[module_name] = mod
    elif module_name == "dataclasses":
        mod = ModuleScope(module_name, None, None)
        indexed_type = PyrexTypes.SpecialPythonTypeConstructor(EncodedString("dataclasses.InitVar"))
        initvar_string = EncodedString("InitVar")
        entry = mod.declare_type(initvar_string, indexed_type, pos = None)
        var_entry = Entry(initvar_string, None, PyrexTypes.py_object_type)
        var_entry.is_pyglobal = True
        var_entry.scope = mod
        entry.as_variable = var_entry
        entry.known_standard_library_import = "%s.InitVar" % module_name
        for name in ["dataclass", "field"]:
            mod.declare_var(EncodedString(name), PyrexTypes.py_object_type, pos=None)
        _known_module_scopes[module_name] = mod
    elif module_name == "functools":
        mod = ModuleScope(module_name, None, None)
        for name in ["total_ordering"]:
            mod.declare_var(EncodedString(name), PyrexTypes.py_object_type, pos=None)
        _known_module_scopes[module_name] = mod

    return mod


def get_known_standard_library_entry(qualified_name):
    name_parts = qualified_name.split(".")
    module_name = EncodedString(name_parts[0])
    rest = name_parts[1:]

    if len(rest) > 1:  # for now, we don't know how to deal with any nested modules
        return None

    mod = get_known_standard_library_module_scope(module_name)

    # eventually handle more sophisticated multiple lookups if needed
    if mod and rest:
        return mod.lookup_here(rest[0])
    return None


def exprnode_to_known_standard_library_name(node, env):
    qualified_name_parts = []
    known_name = None
    while node.is_attribute:
        qualified_name_parts.append(node.attribute)
        node = node.obj
    if node.is_name:
        entry = env.lookup(node.name)
        if entry and entry.known_standard_library_import:
            if get_known_standard_library_entry(
                    entry.known_standard_library_import):
                known_name = entry.known_standard_library_import
            else:
                standard_env = get_known_standard_library_module_scope(
                    entry.known_standard_library_import)
                if standard_env:
                    qualified_name_parts.append(standard_env.name)
                    known_name = ".".join(reversed(qualified_name_parts))
    return known_name
