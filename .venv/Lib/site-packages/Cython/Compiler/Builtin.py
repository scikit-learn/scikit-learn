#
#   Builtin Definitions
#


from .StringEncoding import EncodedString
from .Symtab import BuiltinScope, CClassScope, StructOrUnionScope, ModuleScope, Entry
from .Code import UtilityCode, TempitaUtilityCode, KNOWN_PYTHON_BUILTINS, uncachable_builtins
from .TypeSlots import Signature
from . import PyrexTypes


# C-level implementations of builtin types, functions and methods

iter_next_utility_code = UtilityCode.load("IterNext", "ObjectHandling.c")
getattr_utility_code = UtilityCode.load("GetAttr", "ObjectHandling.c")
getattr3_utility_code = UtilityCode.load("GetAttr3", "Builtins.c")
pyexec_utility_code = UtilityCode.load("PyExec", "Builtins.c")
pyexec_globals_utility_code = UtilityCode.load("PyExecGlobals", "Builtins.c")
globals_utility_code = UtilityCode.load("Globals", "Builtins.c")
range_utility_code = UtilityCode.load("PyRange_Check", "Builtins.c")
include_std_lib_h_utility_code = UtilityCode.load("IncludeStdlibH", "ModuleSetupCode.c")
slice_accessor_utility_code = UtilityCode.load("PySliceAccessors", "Builtins.c")

def make_sequence_multiply_method(typeobj_cname):
    pysequence_multiply_utility_code = TempitaUtilityCode.load(
        "BuiltinSequenceMultiply", "ObjectHandling.c",
        context={'typeobj': typeobj_cname})
    return BuiltinMethod("__mul__",  "Tz",   "T", f"__Pyx_{typeobj_cname}_Multiply",
                         utility_code=pysequence_multiply_utility_code)


# mapping from builtins to their C-level equivalents

class _BuiltinOverride:
    def __init__(self, py_name, args, ret_type, cname, py_equiv="*",
                 utility_code=None, sig=None, func_type=None,
                 is_strict_signature=False, builtin_return_type=None,
                 nogil=None, specialiser=None):
        self.py_name, self.cname, self.py_equiv = py_name, cname, py_equiv
        self.args, self.ret_type = args, ret_type
        self.func_type, self.sig = func_type, sig
        self.builtin_return_type = builtin_return_type
        self.is_strict_signature = is_strict_signature
        self.utility_code = utility_code
        self.nogil = nogil
        self.specialiser = specialiser

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


class BuiltinAttribute:
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
        scope.declare_builtin_cfunction(
            self.py_name, func_type, self.cname, self.py_equiv, self.utility_code,
            specialiser=self.specialiser,
        )


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


class BuiltinProperty:
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


### Special builtin implementations generated at runtime.

def _generate_divmod_function(scope, argument_types):
    if len(argument_types) != 2:
        return None
    type_op1, type_op2 = argument_types

    # Resolve internal typedefs to avoid useless code duplication.
    if type_op1.is_typedef:
        type_op1 = type_op1.resolve_known_type()
    if type_op2.is_typedef:
        type_op2 = type_op2.resolve_known_type()

    if type_op1.is_float or type_op1 is float_type or type_op2.is_float and (type_op1.is_int or type_op1 is int_type):
        impl = "float"
        # TODO: support 'long double'? Currently fails to handle the error return value.
        number_type = PyrexTypes.c_double_type
    elif type_op1.is_int and type_op2.is_int:
        impl = "int"
        number_type = type_op1 if type_op1.rank >= type_op2.rank else type_op2
    else:
        return None

    nogil = scope.nogil
    cfunc_suffix = f"{'nogil_' if nogil else ''}{impl}_{'td_' if number_type.is_typedef else ''}{number_type.specialization_name()}"
    function_cname = f"__Pyx_divmod_{cfunc_suffix}"

    # Reuse an existing specialisation, if available.
    builtin_scope = scope.builtin_scope()
    existing_entry = builtin_scope.lookup_here("divmod")
    if existing_entry is not None:
        for entry in existing_entry.all_alternatives():
            if entry.cname == function_cname:
                return entry

    # Generate a new specialisation.
    ctuple_entry = scope.declare_tuple_type(None, [number_type]*2)
    ctuple_entry.used = True
    return_type = ctuple_entry.type

    function_type = PyrexTypes.CFuncType(
        return_type, [
            PyrexTypes.CFuncTypeArg("a", number_type, None),
            PyrexTypes.CFuncTypeArg("b", number_type, None),
        ],
        exception_value=f"__Pyx_divmod_ERROR_VALUE_{cfunc_suffix}",
        exception_check=True,
        is_strict_signature=True,
        nogil=nogil,
    )

    utility_code = TempitaUtilityCode.load(
        f"divmod_{impl}", "Builtins.c", context={
            'CFUNC_SUFFIX': cfunc_suffix,
            'MATH_SUFFIX': number_type.math_h_modifier if number_type.is_float else '',
            'TYPE': number_type.empty_declaration_code(),
            'RETURN_TYPE': return_type.empty_declaration_code(),
            'NOGIL': nogil,
    })

    entry = builtin_scope.declare_builtin_cfunction(
        "divmod", function_type, function_cname, utility_code=utility_code)

    return entry


### List of builtin functions and their implementation.

builtin_function_table = [
    # name,        args,   return,  C API func,           py equiv = "*"
    BuiltinFunction('abs',        "d",    "d",     "fabs",
                    is_strict_signature=True, nogil=True,
                    utility_code=include_std_lib_h_utility_code),
    BuiltinFunction('abs',        "f",    "f",     "fabsf",
                    is_strict_signature=True, nogil=True,
                    utility_code=include_std_lib_h_utility_code),
    BuiltinFunction('abs',        "i",    "i",     "abs",
                    is_strict_signature=True, nogil=True,
                    utility_code=include_std_lib_h_utility_code),
    BuiltinFunction('abs',        "l",    "l",     "labs",
                    is_strict_signature=True, nogil=True,
                    utility_code=include_std_lib_h_utility_code),
    BuiltinFunction('abs',        None,    None,   "__Pyx_abs_longlong",
                utility_code = UtilityCode.load("abs_longlong", "Builtins.c"),
                func_type = PyrexTypes.CFuncType(
                    PyrexTypes.c_longlong_type, [
                        PyrexTypes.CFuncTypeArg("arg", PyrexTypes.c_longlong_type, None)
                        ],
                    is_strict_signature = True, nogil=True)),
    ] + list(
        BuiltinFunction('abs',        None,    None,   "/*abs_{}*/".format(t.specialization_name()),
                    func_type = PyrexTypes.CFuncType(
                        t,
                        [PyrexTypes.CFuncTypeArg("arg", t, None)],
                        is_strict_signature = True, nogil=True))
                            for t in (PyrexTypes.c_uint_type, PyrexTypes.c_ulong_type, PyrexTypes.c_ulonglong_type)
             ) + list(
        BuiltinFunction('abs',        None,    None,   "__Pyx_c_abs{}".format(t.funcsuffix),
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
    #('aiter',     "",     "",      ""),
    #('anext',     "",     "",      ""),
    BuiltinFunction('ascii',     "O",     "O",      "PyObject_ASCII", builtin_return_type='str'),
    BuiltinFunction('bin',       "O",     "O",      "__Pyx_PyNumber_Bin", builtin_return_type='str',
                    utility_code=UtilityCode(
                        proto="#define __Pyx_PyNumber_Bin(obj) PyNumber_ToBase((obj), 2)",
                        name="PyNumber_Bin")),
    #('breakpoint', "",     "",      ""),
    BuiltinFunction('callable',   "O",    "b",     "__Pyx_PyCallable_Check",
                    utility_code = UtilityCode.load("CallableCheck", "ObjectHandling.c")),
    BuiltinFunction('chr',        "i",    "O",      "PyUnicode_FromOrdinal", builtin_return_type='str'),
    #('compile',   "",     "",      ""), # PyObject* Py_CompileString(    char *str, char *filename, int start)
    BuiltinFunction('delattr',    "OO",   "r",     "__Pyx_PyObject_DelAttr",
                    utility_code=UtilityCode.load("PyObjectDelAttr", "ObjectHandling.c")),
    BuiltinFunction('dir',        "O",    "O",     "PyObject_Dir"),
    BuiltinFunction('divmod',     "OO",   "O",     "PyNumber_Divmod",
                    specialiser=_generate_divmod_function),
    #('enumerate', "",     "",      ""),
    #('eval',      "",     "",      ""),
    BuiltinFunction('exec',       "O",    "O",     "__Pyx_PyExecGlobals",
                    utility_code = pyexec_globals_utility_code),
    BuiltinFunction('exec',       "OO",   "O",     "__Pyx_PyExec2",
                    utility_code = pyexec_utility_code),
    BuiltinFunction('exec',       "OOO",  "O",     "__Pyx_PyExec3",
                    utility_code = pyexec_utility_code),
    #('filter',    "",     "",      ""),
    BuiltinFunction('format',    "OO",     "O",      "PyObject_Format", builtin_return_type='str'),
    BuiltinFunction('format',    "O",      "O",      "__Pyx_PyObject_Format1", builtin_return_type='str',
                    utility_code=UtilityCode(
                        proto="#define __Pyx_PyObject_Format1(obj) PyObject_Format((obj), NULL)",
                        name="PyObject_Format1")),
    BuiltinFunction('getattr3',   "OOO",  "O",     "__Pyx_GetAttr3",     "getattr",
                    utility_code=getattr3_utility_code),  # Pyrex legacy
    BuiltinFunction('getattr',    "OOO",  "O",     "__Pyx_GetAttr3",
                    utility_code=getattr3_utility_code),
    BuiltinFunction('getattr',    "OO",   "O",     "__Pyx_GetAttr",
                    utility_code=getattr_utility_code),
    BuiltinFunction('hasattr',    "OO",   "b",     "__Pyx_HasAttr",
                    utility_code = UtilityCode.load("HasAttr", "Builtins.c")),
    BuiltinFunction('hash',       "O",    "h",     "PyObject_Hash"),
    #('help',      "",     "",      ""),
    BuiltinFunction('hex',       "O",     "O",      "__Pyx_PyNumber_Hex", builtin_return_type='str',
                    utility_code=UtilityCode(
                        proto="#define __Pyx_PyNumber_Hex(obj) PyNumber_ToBase((obj), 16)",
                        name="PyNumber_Hex")),
    #('id',        "",     "",      ""),
    #('input',     "",     "",      ""),
    BuiltinFunction('intern',     "O",    "O",     "__Pyx_Intern",  # Py2 legacy
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
                    utility_code = iter_next_utility_code),
    BuiltinFunction('next',      "OO",    "O",     "__Pyx_PyIter_Next2",
                    utility_code = iter_next_utility_code),
    BuiltinFunction('oct',       "O",     "O",      "__Pyx_PyNumber_Oct", builtin_return_type='str',
                    utility_code=UtilityCode(
                        proto="#define __Pyx_PyNumber_Oct(obj) PyNumber_ToBase((obj), 8)",
                        name="PyNumber_Oct")),
    #('open',      "ss",   "O",     ""),   # no C-API equivalent in Py3
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
    #('print',     "",     "",      ""),
    #('property',  "",     "",      ""),
    BuiltinFunction('reload',     "O",    "O",     "PyImport_ReloadModule"),  # legacy Py2
    BuiltinFunction('repr',       "O",    "O",     "PyObject_Repr", builtin_return_type='str'),
    #('reversed',  "",     "",      ""),
    #('round',     "",     "",      ""),
    BuiltinFunction('setattr',    "OOO",  "r",     "PyObject_SetAttr"),
    #('sorted',    "",     "",      ""),
    #('sum',       "",     "",      ""),
    #('type',       "O",    "O",     "PyObject_Type"),
    BuiltinFunction('unichr',     "i",    "O",      "PyUnicode_FromOrdinal", builtin_return_type='str'),  # legacy Py2
    #('vars',      "",     "",      ""),
    #('zip',       "",     "",      ""),

    # Put in namespace append optimization.
    BuiltinFunction('__Pyx_PyObject_Append', "OO",  "O",     "__Pyx_PyObject_Append"),

    # This is conditionally looked up based on a compiler directive.
    BuiltinFunction('__Pyx_Globals',    "",     "O",     "__Pyx_Globals",
                    utility_code=globals_utility_code),
]


# Builtin types
#  bool
#  bytearray
#  bytes
#  classmethod
#  complex
#  dict
#  enumerate
#  float
#  frozenset
#  int
#  list
#  long
#  memoryview
#  object
#  property
#  range
#  set
#  slice
#  staticmethod
#  str
#  super
#  tuple
#  type

builtin_types_table = [

    ("type",    "&PyType_Type",     []),

    ("bool",   "&PyBool_Type",     []),

    ("int",     "&PyLong_Type",     []),
    ("float",   "&PyFloat_Type",   []),

    ("complex", "&PyComplex_Type", [BuiltinAttribute('cval', field_type_name = 'Py_complex'),
                                    BuiltinAttribute('real', 'cval.real', field_type = PyrexTypes.c_double_type),
                                    BuiltinAttribute('imag', 'cval.imag', field_type = PyrexTypes.c_double_type),
                                    ]),

    ("bytearray", "&PyByteArray_Type", [
                                    make_sequence_multiply_method("PyByteArray_Type"),
                                    ]),
    ("bytes",   "&PyBytes_Type",   [BuiltinMethod("join",  "TO",   "T", "__Pyx_PyBytes_Join",
                                                  utility_code=UtilityCode.load("StringJoin", "StringTools.c")),
                                    make_sequence_multiply_method("PyBytes_Type"),
                                    ]),
    ("str",     "&PyUnicode_Type", [BuiltinMethod("__contains__",  "TO",   "b", "PyUnicode_Contains"),
                                    BuiltinMethod("join",  "TO",   "T", "PyUnicode_Join"),
                                    make_sequence_multiply_method("PyUnicode_Type"),
                                    ]),

    ("tuple",  "&PyTuple_Type",    [make_sequence_multiply_method("PyTuple_Type"),
                                    ]),

    ("list",   "&PyList_Type",     [BuiltinMethod("insert",  "TzO",  "r", "PyList_Insert"),
                                    BuiltinMethod("reverse", "T",    "r", "PyList_Reverse"),
                                    BuiltinMethod("append",  "TO",   "r", "__Pyx_PyList_Append",
                                                  utility_code=UtilityCode.load("ListAppend", "Optimize.c")),
                                    BuiltinMethod("extend",  "TO",   "r", "__Pyx_PyList_Extend",
                                                  utility_code=UtilityCode.load("ListExtend", "Optimize.c")),
                                    make_sequence_multiply_method("PyList_Type"),
                                    ]),

    ("dict",   "&PyDict_Type",     [BuiltinMethod("__contains__",  "TO",   "b", "PyDict_Contains"),
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

    ("range",  "&PyRange_Type",    []),

    ("slice",  "&PySlice_Type",    [BuiltinProperty("start", PyrexTypes.py_object_type, '__Pyx_PySlice_Start',
                                                    utility_code=slice_accessor_utility_code),
                                    BuiltinProperty("stop", PyrexTypes.py_object_type, '__Pyx_PySlice_Stop',
                                                    utility_code=slice_accessor_utility_code),
                                    BuiltinProperty("step", PyrexTypes.py_object_type, '__Pyx_PySlice_Step',
                                                    utility_code=slice_accessor_utility_code),
                                    ]),

    ("set",      "&PySet_Type",    [BuiltinMethod("clear",   "T",  "r", "PySet_Clear"),
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
    ("frozenset", "&PyFrozenSet_Type", []),
    ("BaseException", "((PyTypeObject*)PyExc_BaseException)", []),
    ("Exception", "((PyTypeObject*)PyExc_Exception)", []),
    ("memoryview", "&PyMemoryView_Type", [
        # TODO - format would be nice, but hard to get
        # __len__ can be accessed through a direct lookup of the buffer (but probably in Optimize.c)
        # error checking would ideally be limited api only
        BuiltinProperty("ndim", PyrexTypes.c_int_type, '__Pyx_PyMemoryView_Get_ndim',
                        exception_value=-1, exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="ndim")
                        )
        ),
        BuiltinProperty("readonly", PyrexTypes.c_bint_type, '__Pyx_PyMemoryView_Get_readonly',
                        exception_value=-1, exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="readonly")
                        )
        ),
        BuiltinProperty("itemsize", PyrexTypes.c_py_ssize_t_type, '__Pyx_PyMemoryView_Get_itemsize',
                        exception_value=-1, exception_check=True,
                        utility_code=TempitaUtilityCode.load_cached(
                            "memoryview_get_from_buffer", "Builtins.c",
                            context=dict(name="itemsize")
                        )
        )]
    )
]


types_that_construct_their_instance = frozenset({
    # Some builtin types do not always return an instance of
    # themselves - these do:
    'type', 'bool', 'int', 'float', 'complex',
    'bytes', 'unicode', 'bytearray', 'str',
    'tuple', 'list', 'dict', 'set', 'frozenset',
    'memoryview', 'range',
    # All builtin exception types create their own instance.
    *filter(PyrexTypes.is_exception_type_name, KNOWN_PYTHON_BUILTINS),
})


# When updating this mapping, also update "unsafe_compile_time_methods" below
# if methods are added that are not safe to evaluate at compile time.
inferred_method_return_types = {
    'complex': dict(
        conjugate='complex',
    ),
    'int': dict(
        as_integer_ratio='tuple[int,int]',
        bit_count='T',
        bit_length='T',
        conjugate='T',
        from_bytes='T',  # classmethod
        is_integer='bint',
        to_bytes='bytes',
    ),
    'float': dict(
        as_integer_ratio='tuple[int,int]',
        conjugate='T',
        fromhex='T',  # classmethod
        hex='str',
        is_integer='bint',
    ),
    'list': dict(
        copy='T',
        count='Py_ssize_t',
        index='Py_ssize_t',
    ),
    'tuple': dict(
        count='Py_ssize_t',
        index='Py_ssize_t',
    ),
    'str': dict(
        capitalize='T',
        casefold='T',
        center='T',
        count='Py_ssize_t',
        encode='bytes',
        endswith='bint',
        expandtabs='T',
        find='Py_ssize_t',
        format='T',
        format_map='T',
        index='Py_ssize_t',
        isalnum='bint',
        isalpha='bint',
        isascii='bint',
        isdecimal='bint',
        isdigit='bint',
        isidentifier='bint',
        islower='bint',
        isnumeric='bint',
        isprintable='bint',
        isspace='bint',
        istitle='bint',
        isupper='bint',
        join='T',
        ljust='T',
        lower='T',
        lstrip='T',
        maketrans='dict[int,object]',  # staticmethod
        partition='tuple[T,T,T]',
        removeprefix='T',
        removesuffix='T',
        replace='T',
        rfind='Py_ssize_t',
        rindex='Py_ssize_t',
        rjust='T',
        rpartition='tuple[T,T,T]',
        rsplit='list[T]',
        rstrip='T',
        split='list[T]',
        splitlines='list[T]',
        startswith='bint',
        strip='T',
        swapcase='T',
        title='T',
        translate='T',
        upper='T',
        zfill='T',
    ),
    'bytes': dict(
        capitalize='T',
        center='T',
        count='Py_ssize_t',
        decode='str',
        endswith='bint',
        expandtabs='T',
        find='Py_ssize_t',
        fromhex='T',  # classmethod
        hex='str',
        index='Py_ssize_t',
        isalnum='bint',
        isalpha='bint',
        isascii='bint',
        isdigit='bint',
        islower='bint',
        isspace='bint',
        istitle='bint',
        isupper='bint',
        join='T',
        ljust='T',
        lower='T',
        lstrip='T',
        maketrans='bytes',  # staticmethod
        partition='tuple[T,T,T]',
        removeprefix='T',
        removesuffix='T',
        replace='T',
        rfind='Py_ssize_t',
        rindex='Py_ssize_t',
        rjust='T',
        rpartition='tuple[T,T,T]',
        rsplit='list[T]',
        rstrip='T',
        split='list[T]',
        splitlines='list[T]',
        startswith='bint',
        strip='T',
        swapcase='T',
        title='T',
        translate='T',
        upper='T',
        zfill='T',
    ),
    'bytearray': dict(
        # Inherited from 'bytes' below.
    ),
    'memoryview': dict(
        cast='T',
        hex='str',
        tobytes='bytes',
        tolist='list',
        toreadonly='T',
    ),
    'set': dict(
        copy='T',
        difference='T',
        intersection='T',
        isdisjoint='bint',
        issubset='bint',
        issuperset='bint',
        symmetric_difference='T',
        union='T',
    ),
    'frozenset': dict(
        # Inherited from 'set' below.
    ),
    'dict': dict(
        copy='T',
        fromkeys='T',  # classmethod
        popitem='tuple',
    ),
}

inferred_method_return_types['bytearray'].update(inferred_method_return_types['bytes'])
inferred_method_return_types['frozenset'].update(inferred_method_return_types['set'])


def find_return_type_of_builtin_method(builtin_type, method_name):
    type_name = builtin_type.name
    if type_name in inferred_method_return_types:
        methods = inferred_method_return_types[type_name]
        if method_name in methods:
            return_type_name = methods[method_name]
            if '[' in return_type_name:
                # TODO: Keep the "[...]" part when we add support for generics.
                return_type_name = return_type_name.partition('[')[0]
            if return_type_name == 'T':
                return builtin_type
            if 'T' in return_type_name:
                return_type_name = return_type_name.replace('T', builtin_type.name)
            if return_type_name == 'bint':
                return PyrexTypes.c_bint_type
            elif return_type_name == 'Py_ssize_t':
                return PyrexTypes.c_py_ssize_t_type
            return builtin_scope.lookup(return_type_name).type
    return PyrexTypes.py_object_type


unsafe_compile_time_methods = {
    # We name here only unsafe and non-portable methods if:
    # - the type has a literal representation, allowing for constant folding.
    # - the return type is not None (thus excluding modifier methods)
    #   and is listed in 'inferred_method_return_types' above.
    #
    # See the consistency check in TestBuiltin.py.
    #
    'complex': set(),
    'int': {
        'bit_count',  # Py3.10+
        'from_bytes',  # classmethod
        'is_integer',  # Py3.12+
        'to_bytes',  # changed in Py3.11
    },
    'float': {
        'fromhex',  # classmethod
    },
    'list': {
        'copy',
    },
    'tuple': set(),
    'str': {
        'replace',  # changed in Py3.13+
        'maketrans',  # staticmethod
        'removeprefix',  # Py3.9+
        'removesuffix',  # Py3.9+
    },
    'bytes': {
        'fromhex',  # classmethod
        'maketrans',  # staticmethod
        'removeprefix',  # Py3.9+
        'removesuffix',  # Py3.9+
    },
    'set': set(),
}


def is_safe_compile_time_method(builtin_type_name: str, method_name: str):
    unsafe_methods = unsafe_compile_time_methods.get(builtin_type_name)
    if unsafe_methods is None:
        # Not a literal type.
        return False
    if method_name in unsafe_methods:
        # Not a safe method.
        return False
    known_methods = inferred_method_return_types.get(builtin_type_name)
    if known_methods is None or method_name not in known_methods:
        # Not a known method.
        return False
    return True


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
        if name == 'frozenset':
            objstruct_cname = 'PySetObject'
        elif name == 'bytearray':
            objstruct_cname = 'PyByteArrayObject'
        elif name == 'int':
            objstruct_cname = 'PyLongObject'
        elif name == 'str':
            objstruct_cname = 'PyUnicodeObject'
        elif name == 'bool':
            objstruct_cname = 'PyLongObject'
        elif name == 'BaseException':
            objstruct_cname = "PyBaseExceptionObject"
        elif name == 'Exception':
            objstruct_cname = "PyBaseExceptionObject"
        else:
            objstruct_cname = 'Py%sObject' % name.capitalize()

        utility_code = None
        type_class = PyrexTypes.BuiltinObjectType
        if name in ['dict', 'list', 'set', 'frozenset']:
            type_class = PyrexTypes.BuiltinTypeConstructorObjectType
        elif name == 'tuple':
            type_class = PyrexTypes.PythonTupleTypeConstructor
        elif name == 'range':
            utility_code = range_utility_code
        the_type = builtin_scope.declare_builtin_type(
            name, cname, objstruct_cname=objstruct_cname, type_class=type_class, utility_code=utility_code)
        builtin_types[name] = the_type
        for method in methods:
            method.declare_in_type(the_type)


def init_builtin_exceptions():
    """Declare known builtin Python exceptions as types.
    """
    for name in KNOWN_PYTHON_BUILTINS:
        if name in uncachable_builtins:
            # Exclude builtins specific to later Python versions or platforms.
            continue
        if not PyrexTypes.is_exception_type_name(name):
            continue
        if builtin_scope.lookup_here(name) is not None:
            # Already declared as builtin type above in a more specialised way.
            continue
        utility_code = UtilityCode(
            proto=f"#define __Pyx_PyExc_{name}_Check(obj)  __Pyx_TypeCheck(obj, PyExc_{name})",
            name=f"Py{name}_Check",
        )
        builtin_types[name] = builtin_scope.declare_builtin_type(
            name, f"((PyTypeObject*)PyExc_{name})", utility_code=utility_code)


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
    init_builtin_exceptions()
    init_builtin_funcs()

    entry = builtin_scope.declare_var(
        '__debug__', PyrexTypes.c_const_type(PyrexTypes.c_bint_type),
        pos=None, cname='__pyx_assertions_enabled()', is_cdef=True)
    entry.utility_code = UtilityCode.load_cached("AssertionsEnabled", "Exceptions.c")

    global type_type, list_type, tuple_type, dict_type, set_type, frozenset_type
    global slice_type, range_type
    global bytes_type, unicode_type, bytearray_type
    global float_type, int_type, bool_type, complex_type
    global memoryview_type, py_buffer_type
    global sequence_types
    type_type  = builtin_scope.lookup('type').type
    list_type  = builtin_scope.lookup('list').type
    tuple_type = builtin_scope.lookup('tuple').type
    dict_type  = builtin_scope.lookup('dict').type
    set_type   = builtin_scope.lookup('set').type
    frozenset_type = builtin_scope.lookup('frozenset').type
    slice_type   = builtin_scope.lookup('slice').type
    range_type   = builtin_scope.lookup('range').type

    bytes_type = builtin_scope.lookup('bytes').type
    unicode_type = builtin_scope.lookup('str').type
    bytearray_type = builtin_scope.lookup('bytearray').type
    memoryview_type = builtin_scope.lookup('memoryview').type

    float_type = builtin_scope.lookup('float').type
    int_type = builtin_scope.lookup('int').type
    bool_type  = builtin_scope.lookup('bool').type
    complex_type  = builtin_scope.lookup('complex').type

    sequence_types = (
        list_type,
        tuple_type,
        bytes_type,
        unicode_type,
        bytearray_type,
        memoryview_type,
    )

    # Set up type inference links between equivalent Python/C types
    assert bool_type.name == 'bool', bool_type.name
    bool_type.equivalent_type = PyrexTypes.c_bint_type
    PyrexTypes.c_bint_type.equivalent_type = bool_type

    assert float_type.name == 'float', float_type.name
    float_type.equivalent_type = PyrexTypes.c_double_type
    PyrexTypes.c_double_type.equivalent_type = float_type

    assert complex_type.name == 'complex', complex_type.name
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

        for name in ['ClassVar', 'Optional', 'Union']:
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
