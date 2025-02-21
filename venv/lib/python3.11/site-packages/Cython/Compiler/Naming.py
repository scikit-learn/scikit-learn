#
#   C naming conventions
#
#
#   Prefixes for generating C names.
#   Collected here to facilitate ensuring uniqueness.
#
from .. import __version__

pyrex_prefix    = "__pyx_"
cyversion = __version__.replace('.', '_')


codewriter_temp_prefix = pyrex_prefix + "t_"

temp_prefix       = u"__cyt_"

pyunicode_identifier_prefix = pyrex_prefix + 'U'

builtin_prefix    = pyrex_prefix + "builtin_"
arg_prefix        = pyrex_prefix + "arg_"
genexpr_arg_prefix = pyrex_prefix + "genexpr_arg_"
funcdoc_prefix    = pyrex_prefix + "doc_"
enum_prefix       = pyrex_prefix + "e_"
func_prefix       = pyrex_prefix + "f_"
func_prefix_api   = pyrex_prefix + "api_f_"
pyfunc_prefix     = pyrex_prefix + "pf_"
pywrap_prefix     = pyrex_prefix + "pw_"
genbody_prefix    = pyrex_prefix + "gb_"
gstab_prefix      = pyrex_prefix + "getsets_"
prop_get_prefix   = pyrex_prefix + "getprop_"
const_prefix      = pyrex_prefix + "k_"
py_const_prefix   = pyrex_prefix + "kp_"
label_prefix      = pyrex_prefix + "L"
pymethdef_prefix  = pyrex_prefix + "mdef_"
method_wrapper_prefix = pyrex_prefix + "specialmethod_"
methtab_prefix    = pyrex_prefix + "methods_"
memtab_prefix     = pyrex_prefix + "members_"
objstruct_prefix  = pyrex_prefix + "obj_"
typeptr_prefix    = pyrex_prefix + "ptype_"
prop_set_prefix   = pyrex_prefix + "setprop_"
type_prefix       = pyrex_prefix + "t_"
typeobj_prefix    = pyrex_prefix + "type_"
var_prefix        = pyrex_prefix + "v_"
varptr_prefix     = pyrex_prefix + "vp_"
varptr_prefix_api = pyrex_prefix + "api_vp_"
wrapperbase_prefix= pyrex_prefix + "wrapperbase_"
pybuffernd_prefix   = pyrex_prefix + "pybuffernd_"
pybufferstruct_prefix  = pyrex_prefix + "pybuffer_"
vtable_prefix     = pyrex_prefix + "vtable_"
vtabptr_prefix    = pyrex_prefix + "vtabptr_"
vtabstruct_prefix = pyrex_prefix + "vtabstruct_"
unicode_vtabentry_prefix  = pyrex_prefix + "Uvtabentry_"
# vtab entries aren't normally mangled,
# but punycode names sometimes start with numbers leading to a C syntax error
unicode_structmember_prefix = pyrex_prefix + "Umember_"
# as above -
# not normally mangled but punycode names cause specific problems
opt_arg_prefix    = pyrex_prefix + "opt_args_"
convert_func_prefix = pyrex_prefix + "convert_"
closure_scope_prefix = pyrex_prefix + "scope_"
closure_class_prefix = pyrex_prefix + "scope_struct_"
lambda_func_prefix = pyrex_prefix + "lambda_"
module_is_main   = pyrex_prefix + "module_is_main"
defaults_struct_prefix = pyrex_prefix + "defaults"
dynamic_args_cname = pyrex_prefix + "dynamic_args"

interned_prefixes = {
    'str': pyrex_prefix + "n_",
    'int': pyrex_prefix + "int_",
    'float': pyrex_prefix + "float_",
    'tuple': pyrex_prefix + "tuple_",
    'codeobj': pyrex_prefix + "codeobj_",
    'slice': pyrex_prefix + "slice_",
    'ustring': pyrex_prefix + "ustring_",
    'umethod': pyrex_prefix + "umethod_",
}

ctuple_type_prefix = pyrex_prefix + "ctuple_"
args_cname       = pyrex_prefix + "args"
nargs_cname      = pyrex_prefix + "nargs"
kwvalues_cname   = pyrex_prefix + "kwvalues"
generator_cname  = pyrex_prefix + "generator"
sent_value_cname = pyrex_prefix + "sent_value"
pykwdlist_cname  = pyrex_prefix + "pyargnames"
obj_base_cname   = pyrex_prefix + "base"
builtins_cname   = pyrex_prefix + "b"
preimport_cname  = pyrex_prefix + "i"
moddict_cname    = pyrex_prefix + "d"
dummy_cname      = pyrex_prefix + "dummy"
filename_cname   = pyrex_prefix + "filename"
modulename_cname = pyrex_prefix + "modulename"
filetable_cname  = pyrex_prefix + "f"
intern_tab_cname = pyrex_prefix + "intern_tab"
kwds_cname       = pyrex_prefix + "kwds"
lineno_cname     = pyrex_prefix + "lineno"
clineno_cname    = pyrex_prefix + "clineno"
cfilenm_cname    = pyrex_prefix + "cfilenm"
local_tstate_cname = pyrex_prefix + "tstate"
module_cname     = pyrex_prefix + "m"
modulestate_cname = pyrex_prefix + "mstate"
modulestateglobal_cname = pyrex_prefix + "mstate_global"
moddoc_cname     = pyrex_prefix + "mdoc"
methtable_cname  = pyrex_prefix + "methods"
retval_cname     = pyrex_prefix + "r"
reqd_kwds_cname  = pyrex_prefix + "reqd_kwds"
self_cname       = pyrex_prefix + "self"
stringtab_cname  = pyrex_prefix + "string_tab"
vtabslot_cname   = pyrex_prefix + "vtab"
c_api_tab_cname  = pyrex_prefix + "c_api_tab"
gilstate_cname   = pyrex_prefix + "state"
skip_dispatch_cname = pyrex_prefix + "skip_dispatch"
empty_tuple      = pyrex_prefix + "empty_tuple"
empty_bytes      = pyrex_prefix + "empty_bytes"
empty_unicode    = pyrex_prefix + "empty_unicode"
print_function   = pyrex_prefix + "print"
print_function_kwargs   = pyrex_prefix + "print_kwargs"
cleanup_cname    = pyrex_prefix + "module_cleanup"
pymoduledef_cname = pyrex_prefix + "moduledef"
pymoduledef_slots_cname = pyrex_prefix + "moduledef_slots"
pymodinit_module_arg = pyrex_prefix + "pyinit_module"
pymodule_create_func_cname = pyrex_prefix + "pymod_create"
pymodule_exec_func_cname = pyrex_prefix + "pymod_exec"
optional_args_cname = pyrex_prefix + "optional_args"
import_star      = pyrex_prefix + "import_star"
import_star_set  = pyrex_prefix + "import_star_set"
outer_scope_cname= pyrex_prefix + "outer_scope"
cur_scope_cname  = pyrex_prefix + "cur_scope"
enc_scope_cname  = pyrex_prefix + "enc_scope"
frame_cname      = pyrex_prefix + "frame"
frame_code_cname = pyrex_prefix + "frame_code"
error_without_exception_cname = pyrex_prefix + "error_without_exception"
binding_cfunc    = pyrex_prefix + "binding_PyCFunctionType"
fused_func_prefix = pyrex_prefix + 'fuse_'
fused_dtype_prefix = pyrex_prefix + 'fused_dtype_'
quick_temp_cname = pyrex_prefix + "temp"  # temp variable for quick'n'dirty temping
tp_dict_version_temp = pyrex_prefix + "tp_dict_version"
obj_dict_version_temp = pyrex_prefix + "obj_dict_version"
type_dict_guard_temp = pyrex_prefix + "typedict_guard"
cython_runtime_cname   = pyrex_prefix + "cython_runtime"
cyfunction_type_cname = pyrex_prefix + "CyFunctionType"
fusedfunction_type_cname = pyrex_prefix + "FusedFunctionType"
# the name "dflt" was picked by analogy with the CPython dataclass module which stores
# the default values in variables named f"_dflt_{field.name}" in a hidden scope that's
# passed to the __init__ function. (The name is unimportant to the exact workings though)
dataclass_field_default_cname = pyrex_prefix + "dataclass_dflt"

global_code_object_cache_find = pyrex_prefix + 'find_code_object'
global_code_object_cache_insert = pyrex_prefix + 'insert_code_object'

genexpr_id_ref = 'genexpr'
freelist_name  = 'freelist'
freecount_name = 'freecount'

line_c_macro = "__LINE__"

file_c_macro = "__FILE__"

extern_c_macro  = pyrex_prefix.upper() + "EXTERN_C"

exc_type_name   = pyrex_prefix + "exc_type"
exc_value_name  = pyrex_prefix + "exc_value"
exc_tb_name     = pyrex_prefix + "exc_tb"
exc_lineno_name = pyrex_prefix + "exc_lineno"

parallel_exc_type = pyrex_prefix + "parallel_exc_type"
parallel_exc_value = pyrex_prefix + "parallel_exc_value"
parallel_exc_tb = pyrex_prefix + "parallel_exc_tb"
parallel_filename = pyrex_prefix + "parallel_filename"
parallel_lineno = pyrex_prefix + "parallel_lineno"
parallel_clineno = pyrex_prefix + "parallel_clineno"
parallel_why = pyrex_prefix + "parallel_why"

exc_vars = (exc_type_name, exc_value_name, exc_tb_name)

api_name        = pyrex_prefix + "capi__"

# the h and api guards get changed to:
#  __PYX_HAVE__FILENAME (for ascii filenames)
#  __PYX_HAVE_U_PUNYCODEFILENAME (for non-ascii filenames)
h_guard_prefix   = "__PYX_HAVE_"
api_guard_prefix = "__PYX_HAVE_API_"
api_func_guard   = "__PYX_HAVE_API_FUNC_"

PYX_NAN          = "__PYX_NAN()"

def py_version_hex(major, minor=0, micro=0, release_level=0, release_serial=0):
    return (major << 24) | (minor << 16) | (micro << 8) | (release_level << 4) | (release_serial)

# there's a few places where it's useful to iterate over all of these
used_types_and_macros = [
    (cyfunction_type_cname, '__Pyx_CyFunction_USED'),
    (fusedfunction_type_cname, '__Pyx_FusedFunction_USED'),
    ('__pyx_GeneratorType', '__Pyx_Generator_USED'),
    ('__pyx_IterableCoroutineType', '__Pyx_IterableCoroutine_USED'),
    ('__pyx_CoroutineAwaitType', '__Pyx_Coroutine_USED'),
    ('__pyx_CoroutineType', '__Pyx_Coroutine_USED'),
]
