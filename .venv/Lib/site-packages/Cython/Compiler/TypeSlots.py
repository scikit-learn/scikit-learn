#
#   Tables describing slots in the CPython type object
#   and associated know-how.
#


from . import Naming
from . import PyrexTypes
from .Errors import error, warn_once

import copy

invisible = ['__cinit__', '__dealloc__', '__richcmp__',
             '__nonzero__', '__bool__']

richcmp_special_methods = ['__eq__', '__ne__', '__lt__', '__gt__', '__le__', '__ge__']


class Signature:
    #  Method slot signature descriptor.
    #
    #  has_dummy_arg      boolean
    #  has_generic_args   boolean
    #  fixed_arg_format   string
    #  ret_format         string
    #  error_value        string
    #  use_fastcall       boolean
    #
    #  The formats are strings made up of the following
    #  characters:
    #
    #    'O'  Python object
    #    'T'  Python object of the type of 'self'
    #    'v'  void
    #    'p'  void *
    #    'P'  void **
    #    'i'  int
    #    'b'  bint
    #    'I'  int *
    #    'l'  long
    #    'f'  float
    #    'd'  double
    #    'h'  Py_hash_t
    #    'z'  Py_ssize_t
    #    'Z'  Py_ssize_t *
    #    's'  char *
    #    'S'  char **
    #    'r'  int used only to signal exception
    #    'B'  Py_buffer *
    #    '-'  dummy 'self' argument (not used)
    #    '*'  rest of args passed as generic Python
    #           arg tuple and kw dict (must be last
    #           char in format string)
    #    '?'  optional object arg (currently for pow only)

    format_map = {
        'O': PyrexTypes.py_object_type,
        'v': PyrexTypes.c_void_type,
        'p': PyrexTypes.c_void_ptr_type,
        'P': PyrexTypes.c_void_ptr_ptr_type,
        'i': PyrexTypes.c_int_type,
        'b': PyrexTypes.c_bint_type,
        'I': PyrexTypes.c_int_ptr_type,
        'l': PyrexTypes.c_long_type,
        'f': PyrexTypes.c_float_type,
        'd': PyrexTypes.c_double_type,
        'h': PyrexTypes.c_py_hash_t_type,
        'z': PyrexTypes.c_py_ssize_t_type,
        'Z': PyrexTypes.c_py_ssize_t_ptr_type,
        's': PyrexTypes.c_char_ptr_type,
        'S': PyrexTypes.c_char_ptr_ptr_type,
        'r': PyrexTypes.c_returncode_type,
        'B': PyrexTypes.c_py_buffer_ptr_type,
        '?': PyrexTypes.py_object_type
        # 'T', '-' and '*' are handled otherwise
        # and are not looked up in here
    }

    type_to_format_map = {type_: format_ for format_, type_ in format_map.items()}

    error_value_map = {
        'O': "NULL",
        'T': "NULL",
        'i': "-1",
        'b': "-1",
        'l': "-1",
        'r': "-1",
        'h': "-1",
        'z': "-1",
    }

    # Use METH_FASTCALL instead of METH_VARARGS
    use_fastcall = False

    def __init__(self, arg_format, ret_format, nogil=False):
        self.has_dummy_arg = False
        self.has_generic_args = False
        self.optional_object_arg_count = 0
        if arg_format[:1] == '-':
            self.has_dummy_arg = True
            arg_format = arg_format[1:]
        if arg_format[-1:] == '*':
            self.has_generic_args = True
            arg_format = arg_format[:-1]
        if arg_format[-1:] == '?':
            self.optional_object_arg_count += 1
        self.fixed_arg_format = arg_format
        self.ret_format = ret_format
        self.error_value = self.error_value_map.get(ret_format, None)
        self.exception_check = ret_format != 'r' and self.error_value is not None
        self.is_staticmethod = False
        self.nogil = nogil

    def __repr__(self):
        return '<Signature[%s(%s%s)]>' % (
            self.ret_format,
            ', '.join(self.fixed_arg_format),
            '*' if self.has_generic_args else '')

    def min_num_fixed_args(self):
        return self.max_num_fixed_args() - self.optional_object_arg_count

    def max_num_fixed_args(self):
        return len(self.fixed_arg_format)

    def is_self_arg(self, i):
        # argument is 'self' for methods or 'class' for classmethods
        return self.fixed_arg_format[i] == 'T'

    def returns_self_type(self):
        # return type is same as 'self' argument type
        return self.ret_format == 'T'

    def fixed_arg_type(self, i):
        return self.format_map[self.fixed_arg_format[i]]

    def return_type(self):
        return self.format_map[self.ret_format]

    def format_from_type(self, arg_type):
        if arg_type.is_pyobject:
            arg_type = PyrexTypes.py_object_type
        return self.type_to_format_map[arg_type]

    def exception_value(self):
        return self.error_value_map.get(self.ret_format)

    def function_type(self, self_arg_override=None):
        #  Construct a C function type descriptor for this signature
        args = []
        for i in range(self.max_num_fixed_args()):
            if self_arg_override is not None and self.is_self_arg(i):
                assert isinstance(self_arg_override, PyrexTypes.CFuncTypeArg)
                args.append(self_arg_override)
            else:
                arg_type = self.fixed_arg_type(i)
                args.append(PyrexTypes.CFuncTypeArg("", arg_type, None))
        if self_arg_override is not None and self.returns_self_type():
            ret_type = self_arg_override.type
        else:
            ret_type = self.return_type()
        exc_value = self.exception_value()
        return PyrexTypes.CFuncType(
            ret_type, args, exception_value=exc_value,
            exception_check=self.exception_check,
            nogil=self.nogil)

    def method_flags(self):
        if self.ret_format == "O":
            full_args = self.fixed_arg_format
            if self.has_dummy_arg:
                full_args = "O" + full_args
            if full_args in ["O", "T"]:
                if not self.has_generic_args:
                    return [method_noargs]
                elif self.use_fastcall:
                    return [method_fastcall, method_keywords]
                else:
                    return [method_varargs, method_keywords]
            elif full_args in ["OO", "TO"] and not self.has_generic_args:
                return [method_onearg]

            if self.is_staticmethod:
                if self.use_fastcall:
                    return [method_fastcall, method_keywords]
                else:
                    return [method_varargs, method_keywords]
        return None

    def method_function_type(self):
        # Return the C function type
        mflags = self.method_flags()
        kw = "WithKeywords" if (method_keywords in mflags) else ""
        for m in mflags:
            if m == method_noargs or m == method_onearg:
                return "PyCFunction"
            if m == method_varargs:
                return "PyCFunction" + kw
            if m == method_fastcall:
                return "__Pyx_PyCFunction_FastCall" + kw
        return None

    def with_fastcall(self):
        # Return a copy of this Signature with use_fastcall=True
        sig = copy.copy(self)
        sig.use_fastcall = True
        return sig

    @property
    def fastvar(self):
        # Used to select variants of functions, one dealing with METH_VARARGS
        # and one dealing with __Pyx_METH_FASTCALL
        if self.use_fastcall:
            return "FASTCALL"
        else:
            return "VARARGS"


class SlotDescriptor:
    #  Abstract base class for type slot descriptors.
    #
    #  slot_name    string           Member name of the slot in the type object
    #  is_initialised_dynamically    Is initialised by code in the module init function
    #  is_inherited                  Is inherited by subtypes (see PyType_Ready())
    #  ifdef                         Full #ifdef string that slot is wrapped in. Using this causes flags to be ignored.
    #  used_ifdef                    Full #ifdef string that the slot value is wrapped in (otherwise it is assigned NULL)
    #                                Unlike "ifdef" the slot is defined and this just controls if it receives a value

    def __init__(self, slot_name, dynamic=False, inherited=False,
                 ifdef=None, is_binop=False,
                 used_ifdef=None):
        self.slot_name = slot_name
        self.is_initialised_dynamically = dynamic
        self.is_inherited = inherited
        self.ifdef = ifdef
        self.used_ifdef = used_ifdef
        self.is_binop = is_binop

    def slot_code(self, scope):
        raise NotImplementedError()

    def spec_value(self, scope):
        return self.slot_code(scope)

    def preprocessor_guard_code(self):
        ifdef = self.ifdef
        guard = None
        if ifdef:
            guard = "#if %s" % ifdef
        return guard

    def generate_spec(self, scope, code):
        if self.is_initialised_dynamically:
            return
        value = self.spec_value(scope)
        if value == "0":
            return
        preprocessor_guard = self.preprocessor_guard_code()
        if not preprocessor_guard:
            if self.slot_name.startswith(('bf_', 'am_')):
                # The buffer protocol requires Limited API 3.11 and 'am_send' requires 3.10,
                # so check if the spec slots are available.
                preprocessor_guard = "#if defined(Py_%s)" % self.slot_name
        if preprocessor_guard:
            code.putln(preprocessor_guard)
        if self.used_ifdef:
            # different from preprocessor guard - this defines if we *want* to define it,
            # rather than if the slot exists
            code.putln(f"#if {self.used_ifdef}")
        code.putln("{Py_%s, (void *)%s}," % (self.slot_name, value))
        if self.used_ifdef:
            code.putln("#endif")
        if preprocessor_guard:
            code.putln("#endif")

    def generate(self, scope, code):
        preprocessor_guard = self.preprocessor_guard_code()
        if preprocessor_guard:
            code.putln(preprocessor_guard)

        end_pypy_guard = False
        if self.is_initialised_dynamically:
            value = "0"
        else:
            value = self.slot_code(scope)
            if value == "0" and self.is_inherited:
                # PyPy currently has a broken PyType_Ready() that fails to
                # inherit some slots.  To work around this, we explicitly
                # set inherited slots here, but only in PyPy since CPython
                # handles this better than we do (except for buffer slots in type specs).
                inherited_value = value
                current_scope = scope
                while (inherited_value == "0"
                       and current_scope.parent_type
                       and current_scope.parent_type.base_type
                       and current_scope.parent_type.base_type.scope):
                    current_scope = current_scope.parent_type.base_type.scope
                    inherited_value = self.slot_code(current_scope)
                if inherited_value != "0":
                    # we always need inherited buffer slots for the type spec
                    is_buffer_slot = int(self.slot_name in ("bf_getbuffer", "bf_releasebuffer"))
                    code.putln("#if CYTHON_COMPILING_IN_PYPY || %d" % is_buffer_slot)
                    code.putln("%s, /*%s*/" % (inherited_value, self.slot_name))
                    code.putln("#else")
                    end_pypy_guard = True

        if self.used_ifdef:
            code.putln("#if %s" % self.used_ifdef)
        code.putln("%s, /*%s*/" % (value, self.slot_name))
        if self.used_ifdef:
            code.putln("#else")
            code.putln("NULL, /*%s*/" % self.slot_name)
            code.putln("#endif")

        if end_pypy_guard:
            code.putln("#endif")

        if preprocessor_guard:
            code.putln("#endif")

    # Some C implementations have trouble statically
    # initialising a global with a pointer to an extern
    # function, so we initialise some of the type slots
    # in the module init function instead.

    def generate_dynamic_init_code(self, scope, code):
        if self.is_initialised_dynamically:
            self.generate_set_slot_code(
                self.slot_code(scope), scope, code)

    def generate_set_slot_code(self, value, scope, code):
        if value == "0":
            return

        if scope.parent_type.typeptr_cname:
            target = "%s->%s" % (
                code.typeptr_cname_in_module_state(scope.parent_type), self.slot_name)
        else:
            assert scope.parent_type.typeobj_cname
            target = "%s.%s" % (
                code.name_in_module_state(scope.parent_type.typeobj_cname), self.slot_name)

        code.putln("%s = %s;" % (target, value))


class FixedSlot(SlotDescriptor):
    #  Descriptor for a type slot with a fixed value.
    #
    #  value        string

    def __init__(self, slot_name, value, ifdef=None):
        SlotDescriptor.__init__(self, slot_name, ifdef=ifdef)
        self.value = value

    def slot_code(self, scope):
        return self.value


class EmptySlot(FixedSlot):
    #  Descriptor for a type slot whose value is always 0.

    def __init__(self, slot_name, ifdef=None):
        FixedSlot.__init__(self, slot_name, "0", ifdef=ifdef)


class MethodSlot(SlotDescriptor):
    #  Type slot descriptor for a user-definable method.
    #
    #  signature    Signature
    #  method_name  string           The __xxx__ name of the method
    #  alternatives [string]         Alternative list of __xxx__ names for the method

    def __init__(self, signature, slot_name, method_name, method_name_to_slot,
                 fallback=None, ifdef=None, inherited=True):
        SlotDescriptor.__init__(self, slot_name,
                                ifdef=ifdef, inherited=inherited)
        self.signature = signature
        self.slot_name = slot_name
        self.method_name = method_name
        self.alternatives = []
        method_name_to_slot[method_name] = self
        #
        if fallback:
            self.alternatives.append(fallback)

    def slot_code(self, scope):
        entry = scope.lookup_here(self.method_name)
        if entry and entry.is_special and entry.func_cname:
            for method_name in self.alternatives:
                alt_entry = scope.lookup_here(method_name)
                if alt_entry:
                    warn_once(alt_entry.pos,
                              f"{method_name} was removed in Python 3; ignoring it and using {self.method_name} instead",
                              2)
            return entry.func_cname
        for method_name in self.alternatives:
            entry = scope.lookup_here(method_name)
            if entry and entry.is_special and entry.func_cname:
                warn_once(entry.pos,
                          f"{method_name} was removed in Python 3; use {self.method_name} instead",
                          2)
                return entry.func_cname
        return "0"


class InternalMethodSlot(SlotDescriptor):
    #  Type slot descriptor for a method which is always
    #  synthesized by Cython.
    #
    #  slot_name    string           Member name of the slot in the type object

    def __init__(self, slot_name, **kargs):
        SlotDescriptor.__init__(self, slot_name, **kargs)

    def slot_code(self, scope):
        return scope.mangle_internal(self.slot_name)


class GCDependentSlot(InternalMethodSlot):
    #  Descriptor for a slot whose value depends on whether
    #  the type participates in GC.

    def __init__(self, slot_name, **kargs):
        InternalMethodSlot.__init__(self, slot_name, **kargs)

    def slot_code(self, scope):
        # We treat external types as needing gc, but don't generate a slot code
        # because we don't know it to be able to call it directly.
        if not scope.needs_gc() or scope.parent_type.is_external:
            return "0"
        if not scope.has_cyclic_pyobject_attrs:
            # if the type does not have GC relevant object attributes, it can
            # delegate GC methods to its parent - iff the parent functions
            # are defined in the same module
            parent_type_scope = scope.parent_type.base_type.scope
            if scope.parent_scope is parent_type_scope.parent_scope:
                entry = scope.parent_scope.lookup_here(scope.parent_type.base_type.name)
                if entry.visibility != 'extern':
                    return self.slot_code(parent_type_scope)
        return InternalMethodSlot.slot_code(self, scope)


class GCClearReferencesSlot(GCDependentSlot):

    def slot_code(self, scope):
        if scope.needs_tp_clear():
            return GCDependentSlot.slot_code(self, scope)
        return "0"


class ConstructorSlot(InternalMethodSlot):
    #  Descriptor for tp_new and tp_dealloc.

    def __init__(self, slot_name, method=None, **kargs):
        InternalMethodSlot.__init__(self, slot_name, **kargs)
        self.method = method

    def _needs_own(self, scope):
        if (scope.parent_type.base_type
                and not scope.has_pyobject_attrs
                and not scope.has_memoryview_attrs
                and not scope.has_explicitly_constructable_attrs
                and not (self.slot_name == 'tp_new' and scope.parent_type.vtabslot_cname)):
            entry = scope.lookup_here(self.method) if self.method else None
            if not (entry and entry.is_special):
                return False
        # Unless we can safely delegate to the parent, all types need a tp_new().
        return True

    def _parent_slot_function(self, scope):
        parent_type_scope = scope.parent_type.base_type.scope
        if scope.parent_scope is parent_type_scope.parent_scope:
            entry = scope.parent_scope.lookup_here(scope.parent_type.base_type.name)
            if entry.visibility != 'extern':
                return self.slot_code(parent_type_scope)
        return None

    def slot_code(self, scope):
        if not self._needs_own(scope):
            # if the type does not have object attributes, it can
            # delegate GC methods to its parent - iff the parent
            # functions are defined in the same module
            slot_code = self._parent_slot_function(scope)
            if slot_code is not None:
                return slot_code
        return InternalMethodSlot.slot_code(self, scope)

    def generate_dynamic_init_code(self, scope, code):
        if self.slot_code(scope) != '0':
            return
        # If we don't have our own slot function and don't know the
        # parent function statically, copy it dynamically.
        base_type = scope.parent_type.base_type
        if base_type.typeptr_cname:
            base_typeptr_cname = code.typeptr_cname_in_module_state(base_type)
            src = '%s->%s' % (base_typeptr_cname, self.slot_name)
        elif base_type.is_extension_type and base_type.typeobj_cname:
            src = '%s.%s' % (code.typeptr_cname_in_module_state(base_type), self.slot_name)
        else:
            return

        self.generate_set_slot_code(src, scope, code)


class SyntheticSlot(InternalMethodSlot):
    #  Type slot descriptor for a synthesized method which
    #  dispatches to one or more user-defined methods depending
    #  on its arguments. If none of the relevant methods are
    #  defined, the method will not be synthesized and an
    #  alternative default value will be placed in the type
    #  slot.

    def __init__(self, slot_name, user_methods, default_value, **kargs):
        InternalMethodSlot.__init__(self, slot_name, **kargs)
        self.user_methods = user_methods
        self.default_value = default_value

    def slot_code(self, scope):
        if scope.defines_any_special(self.user_methods):
            return InternalMethodSlot.slot_code(self, scope)
        else:
            return self.default_value

    def spec_value(self, scope):
        return self.slot_code(scope)


class BinopSlot(SyntheticSlot):
    def __init__(self, signature, slot_name, left_method, method_name_to_slot, **kargs):
        assert left_method.startswith('__')
        right_method = '__r' + left_method[2:]
        SyntheticSlot.__init__(
                self, slot_name, [left_method, right_method], "0", is_binop=True, **kargs)
        # MethodSlot causes special method registration.
        self.left_slot = MethodSlot(signature, "", left_method, method_name_to_slot, **kargs)
        self.right_slot = MethodSlot(signature, "", right_method, method_name_to_slot, **kargs)


class RichcmpSlot(MethodSlot):
    def slot_code(self, scope):
        entry = scope.lookup_here(self.method_name)
        if entry and entry.is_special and entry.func_cname:
            return entry.func_cname
        elif scope.defines_any_special(richcmp_special_methods):
            return scope.mangle_internal(self.slot_name)
        else:
            return "0"


class TypeFlagsSlot(SlotDescriptor):
    #  Descriptor for the type flags slot.

    def slot_code(self, scope):
        value = "Py_TPFLAGS_DEFAULT"
        if scope.directives['type_version_tag']:
            # No longer used since Py3.11.
            value += "|Py_TPFLAGS_HAVE_VERSION_TAG"
        else:
            # Used to be in 'Py_TPFLAGS_DEFAULT' up to Py3.10.
            value = f"({value}&~Py_TPFLAGS_HAVE_VERSION_TAG)"
        value += "|Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_HAVE_NEWBUFFER"
        if not scope.parent_type.is_final_type:
            value += "|Py_TPFLAGS_BASETYPE"
        if scope.needs_gc():
            value += "|Py_TPFLAGS_HAVE_GC"
        if scope.parent_type.has_sequence_flag:
            value += "|Py_TPFLAGS_SEQUENCE"
        return value

    def generate_spec(self, scope, code):
        # Flags are stored in the PyType_Spec, not in a PyType_Slot.
        return


class DocStringSlot(SlotDescriptor):
    #  Descriptor for the docstring slot.

    def slot_code(self, scope):
        doc = scope.doc
        if doc is None:
            return "0"
        if doc.is_unicode:
            doc = doc.as_utf8_string()
        return "PyDoc_STR(%s)" % doc.as_c_string_literal()


class SuiteSlot(SlotDescriptor):
    #  Descriptor for a substructure of the type object.
    #
    #  sub_slots   [SlotDescriptor]

    def __init__(self, sub_slots, slot_type, slot_name, substructures, ifdef=None, cast_cname=None):
        SlotDescriptor.__init__(self, slot_name, ifdef=ifdef)
        self.sub_slots = sub_slots
        self.slot_type = slot_type
        self.cast_cname = cast_cname
        substructures.append(self)

    def is_empty(self, scope):
        for slot in self.sub_slots:
            if slot.slot_code(scope) != "0":
                return False
        return True

    def substructure_cname(self, scope):
        return "%s%s_%s" % (Naming.pyrex_prefix, self.slot_name, scope.class_name)

    def slot_code(self, scope):
        if not self.is_empty(scope):
            cast = ""
            if self.cast_cname:
                cast = f"({self.cast_cname}*)"
            return f"{cast}&{self.substructure_cname(scope)}"
        return "0"

    def generate_substructure(self, scope, code):
        if not self.is_empty(scope):
            code.putln("")
            if self.ifdef:
                code.putln("#if %s" % self.ifdef)
            code.putln(
                "static %s %s = {" % (
                    self.slot_type,
                    self.substructure_cname(scope)))
            for slot in self.sub_slots:
                slot.generate(scope, code)
            code.putln("};")
            if self.ifdef:
                code.putln("#endif")

    def generate_spec(self, scope, code):
        for slot in self.sub_slots:
            slot.generate_spec(scope, code)

class MethodTableSlot(SlotDescriptor):
    #  Slot descriptor for the method table.

    def slot_code(self, scope):
        if scope.pyfunc_entries:
            return scope.method_table_cname
        else:
            return "0"


class MemberTableSlot(SlotDescriptor):
    #  Slot descriptor for the table of Python-accessible attributes.

    def slot_code(self, scope):
        # Only used in specs.
        return "0"

    def get_member_specs(self, scope):
        return [
            get_slot_by_name("tp_dictoffset", scope.directives).members_slot_value(scope),
            #get_slot_by_name("tp_weaklistoffset").spec_value(scope),
        ]

    def is_empty(self, scope):
        for member_entry in self.get_member_specs(scope):
            if member_entry:
                return False
        return True

    def substructure_cname(self, scope):
        return "%s%s_%s" % (Naming.pyrex_prefix, self.slot_name, scope.class_name)

    def generate_substructure_spec(self, scope, code):
        if self.is_empty(scope):
            return
        from .Code import UtilityCode
        code.globalstate.use_utility_code(UtilityCode.load_cached("IncludeStructmemberH", "ModuleSetupCode.c"))

        code.putln("static struct PyMemberDef %s[] = {" % self.substructure_cname(scope))
        for member_entry in self.get_member_specs(scope):
            if member_entry:
                code.putln(member_entry)
        code.putln("{NULL, 0, 0, 0, NULL}")
        code.putln("};")

    def spec_value(self, scope):
        if self.is_empty(scope):
            return "0"
        return self.substructure_cname(scope)


class GetSetSlot(SlotDescriptor):
    #  Slot descriptor for the table of attribute get & set methods.

    def slot_code(self, scope):
        if scope.property_entries:
            return scope.getset_table_cname
        else:
            return "0"


class BaseClassSlot(SlotDescriptor):
    #  Slot descriptor for the base class slot.

    def __init__(self, name):
        SlotDescriptor.__init__(self, name, dynamic=True)

    def generate_dynamic_init_code(self, scope, code):
        base_type = scope.parent_type.base_type
        if base_type:
            base_typeptr_cname = code.typeptr_cname_in_module_state(base_type)
            code.putln("%s->%s = %s;" % (
                code.typeptr_cname_in_module_state(scope.parent_type),
                self.slot_name,
                base_typeptr_cname))


class DictOffsetSlot(SlotDescriptor):
    #  Slot descriptor for a class' dict offset, for dynamic attributes.

    def slot_code(self, scope):
        dict_entry = scope.lookup_here("__dict__") if not scope.is_closure_class_scope else None
        if dict_entry and dict_entry.is_variable:
            from . import Builtin
            if dict_entry.type is not Builtin.dict_type:
                error(dict_entry.pos, "__dict__ slot must be of type 'dict'")
                return "0"
            type = scope.parent_type
            if type.typedef_flag:
                objstruct = type.objstruct_cname
            else:
                objstruct = "struct %s" % type.objstruct_cname
            return ("offsetof(%s, %s)" % (
                        objstruct,
                        dict_entry.cname))
        else:
            return "0"

    def members_slot_value(self, scope):
        dict_offset = self.slot_code(scope)
        if dict_offset == "0":
            return None
        return '{"__dictoffset__", T_PYSSIZET, %s, READONLY, NULL},' % dict_offset

## The following slots are (or could be) initialised with an
## extern function pointer.
#
#slots_initialised_from_extern = (
#    "tp_free",
#)

#------------------------------------------------------------------------------------------
#
#  Utility functions for accessing slot table data structures
#
#------------------------------------------------------------------------------------------


def get_property_accessor_signature(name):
    #  Return signature of accessor for an extension type
    #  property, else None.
    return property_accessor_signatures.get(name)


def get_base_slot_function(scope, slot):
    #  Returns the function implementing this slot in the baseclass.
    #  This is useful for enabling the compiler to optimize calls
    #  that recursively climb the class hierarchy.
    base_type = scope.parent_type.base_type
    if base_type and scope.parent_scope is base_type.scope.parent_scope:
        parent_slot = slot.slot_code(base_type.scope)
        if parent_slot != '0':
            entry = scope.parent_scope.lookup_here(scope.parent_type.base_type.name)
            if entry.visibility != 'extern':
                return parent_slot
    return None


def get_slot_function(scope, slot):
    #  Returns the function implementing this slot in the baseclass.
    #  This is useful for enabling the compiler to optimize calls
    #  that recursively climb the class hierarchy.
    slot_code = slot.slot_code(scope)
    if slot_code != '0':
        entry = scope.parent_scope.lookup_here(scope.parent_type.name)
        if entry.visibility != 'extern':
            return slot_code
    return None


def get_slot_by_name(slot_name, compiler_directives):
    # For now, only search the type struct, no referenced sub-structs.
    for slot in get_slot_table(compiler_directives).slot_table:
        if slot.slot_name == slot_name:
            return slot
    assert False, "Slot not found: %s" % slot_name


def get_slot_code_by_name(scope, slot_name):
    slot = get_slot_by_name(slot_name, scope.directives)
    return slot.slot_code(scope)

def is_binop_number_slot(name):
    """
    Tries to identify __add__/__radd__ and friends (so the METH_COEXIST flag can be applied).

    There's no great consequence if it inadvertently identifies a few other methods
    so just use a simple rule rather than an exact list.
    """
    slot_table = get_slot_table(None)
    for meth in get_slot_table(None).PyNumberMethods:
        if meth.is_binop and name in meth.user_methods:
            return True
    return False


#------------------------------------------------------------------------------------------
#
#  Signatures for generic Python functions and methods.
#
#------------------------------------------------------------------------------------------

pyfunction_signature = Signature("-*", "O")
pymethod_signature = Signature("T*", "O")

#------------------------------------------------------------------------------------------
#
#  Signatures for simple Python functions.
#
#------------------------------------------------------------------------------------------

pyfunction_noargs = Signature("-", "O")
pyfunction_onearg = Signature("-O", "O")

#------------------------------------------------------------------------------------------
#
#  Signatures for the various kinds of function that
#  can appear in the type object and its substructures.
#
#------------------------------------------------------------------------------------------

unaryfunc = Signature("T", "O")            # typedef PyObject * (*unaryfunc)(PyObject *);
binaryfunc = Signature("OO", "O")          # typedef PyObject * (*binaryfunc)(PyObject *, PyObject *);
ibinaryfunc = Signature("TO", "O")         # typedef PyObject * (*binaryfunc)(PyObject *, PyObject *);
powternaryfunc = Signature("OO?", "O")     # typedef PyObject * (*ternaryfunc)(PyObject *, PyObject *, PyObject *);
ipowternaryfunc = Signature("TO?", "O")    # typedef PyObject * (*ternaryfunc)(PyObject *, PyObject *, PyObject *);
callfunc = Signature("T*", "O")            # typedef PyObject * (*ternaryfunc)(PyObject *, PyObject *, PyObject *);
inquiry = Signature("T", "i")              # typedef int (*inquiry)(PyObject *);
lenfunc = Signature("T", "z")              # typedef Py_ssize_t (*lenfunc)(PyObject *);

                                           # typedef int (*coercion)(PyObject **, PyObject **);
intargfunc = Signature("Ti", "O")          # typedef PyObject *(*intargfunc)(PyObject *, int);
ssizeargfunc = Signature("Tz", "O")        # typedef PyObject *(*ssizeargfunc)(PyObject *, Py_ssize_t);
intintargfunc = Signature("Tii", "O")      # typedef PyObject *(*intintargfunc)(PyObject *, int, int);
ssizessizeargfunc = Signature("Tzz", "O")  # typedef PyObject *(*ssizessizeargfunc)(PyObject *, Py_ssize_t, Py_ssize_t);
intobjargproc = Signature("TiO", 'r')      # typedef int(*intobjargproc)(PyObject *, int, PyObject *);
ssizeobjargproc = Signature("TzO", 'r')    # typedef int(*ssizeobjargproc)(PyObject *, Py_ssize_t, PyObject *);
intintobjargproc = Signature("TiiO", 'r')  # typedef int(*intintobjargproc)(PyObject *, int, int, PyObject *);
ssizessizeobjargproc = Signature("TzzO", 'r')  # typedef int(*ssizessizeobjargproc)(PyObject *, Py_ssize_t, Py_ssize_t, PyObject *);

intintargproc = Signature("Tii", 'r')
ssizessizeargproc = Signature("Tzz", 'r')
objargfunc = Signature("TO", "O")
objobjargproc = Signature("TOO", 'r')      # typedef int (*objobjargproc)(PyObject *, PyObject *, PyObject *);
readbufferproc = Signature("TzP", "z")     # typedef Py_ssize_t (*readbufferproc)(PyObject *, Py_ssize_t, void **);
writebufferproc = Signature("TzP", "z")    # typedef Py_ssize_t (*writebufferproc)(PyObject *, Py_ssize_t, void **);
segcountproc = Signature("TZ", "z")        # typedef Py_ssize_t (*segcountproc)(PyObject *, Py_ssize_t *);
charbufferproc = Signature("TzS", "z")     # typedef Py_ssize_t (*charbufferproc)(PyObject *, Py_ssize_t, char **);
objargproc = Signature("TO", 'r')          # typedef int (*objobjproc)(PyObject *, PyObject *);
                                           # typedef int (*visitproc)(PyObject *, void *);
                                           # typedef int (*traverseproc)(PyObject *, visitproc, void *);

destructor = Signature("T", "v")           # typedef void (*destructor)(PyObject *);
# printfunc = Signature("TFi", 'r')        # typedef int (*printfunc)(PyObject *, FILE *, int);
                                           # typedef PyObject *(*getattrfunc)(PyObject *, char *);
getattrofunc = Signature("TO", "O")        # typedef PyObject *(*getattrofunc)(PyObject *, PyObject *);
                                           # typedef int (*setattrfunc)(PyObject *, char *, PyObject *);
setattrofunc = Signature("TOO", 'r')       # typedef int (*setattrofunc)(PyObject *, PyObject *, PyObject *);
delattrofunc = Signature("TO", 'r')
cmpfunc = Signature("TO", "i")             # typedef int (*cmpfunc)(PyObject *, PyObject *);
reprfunc = Signature("T", "O")             # typedef PyObject *(*reprfunc)(PyObject *);
hashfunc = Signature("T", "h")             # typedef Py_hash_t (*hashfunc)(PyObject *);
richcmpfunc = Signature("TOi", "O")        # typedef PyObject *(*richcmpfunc) (PyObject *, PyObject *, int);
getiterfunc = Signature("T", "O")          # typedef PyObject *(*getiterfunc) (PyObject *);
iternextfunc = Signature("T", "O")         # typedef PyObject *(*iternextfunc) (PyObject *);
descrgetfunc = Signature("TOO", "O")       # typedef PyObject *(*descrgetfunc) (PyObject *, PyObject *, PyObject *);
descrsetfunc = Signature("TOO", 'r')       # typedef int (*descrsetfunc) (PyObject *, PyObject *, PyObject *);
descrdelfunc = Signature("TO", 'r')
initproc = Signature("T*", 'r')            # typedef int (*initproc)(PyObject *, PyObject *, PyObject *);
                                           # typedef PyObject *(*newfunc)(struct _typeobject *, PyObject *, PyObject *);
                                           # typedef PyObject *(*allocfunc)(struct _typeobject *, int);

getbufferproc = Signature("TBi", "r")      # typedef int (*getbufferproc)(PyObject *, Py_buffer *, int);
releasebufferproc = Signature("TB", "v")   # typedef void (*releasebufferproc)(PyObject *, Py_buffer *);

# typedef PySendResult (*sendfunc)(PyObject* iter, PyObject* value, PyObject** result);
sendfunc = PyrexTypes.CPtrType(PyrexTypes.CFuncType(
    return_type=PyrexTypes.PySendResult_type,
    args=[
        PyrexTypes.CFuncTypeArg("iter", PyrexTypes.py_object_type),
        PyrexTypes.CFuncTypeArg("value", PyrexTypes.py_object_type),
        PyrexTypes.CFuncTypeArg("result", PyrexTypes.CPtrType(PyrexTypes.py_objptr_type)),
    ],
    exception_value="PYGEN_ERROR",
    exception_check=True,  # we allow returning PYGEN_ERROR without GeneratorExit / StopIteration
))


#------------------------------------------------------------------------------------------
#
#  Signatures for accessor methods of properties.
#
#------------------------------------------------------------------------------------------

property_accessor_signatures = {
    '__get__': Signature("T", "O"),
    '__set__': Signature("TO", 'r'),
    '__del__': Signature("T", 'r')
}

#------------------------------------------------------------------------------------------
#
#  The main slot table. This table contains descriptors for all the
#  top-level type slots, beginning with tp_dealloc, in the order they
#  appear in the type object.
#
# It depends on some compiler directives (currently c_api_binop_methods), so the
# slot tables for each set of compiler directives are generated lazily and put in
# the _slot_table_dict
#
#------------------------------------------------------------------------------------------

class SlotTable:
    def __init__(self, old_binops):
        # The following dictionary maps __xxx__ method names to slot descriptors.
        method_name_to_slot = {}
        self._get_slot_by_method_name = method_name_to_slot.get
        self.substructures = []   # List of all SuiteSlot instances

        bf = binaryfunc if old_binops else ibinaryfunc
        ptf = powternaryfunc if old_binops else ipowternaryfunc

        #  Descriptor tables for the slots of the various type object
        #  substructures, in the order they appear in the structure.
        self.PyNumberMethods = (
            BinopSlot(bf, "nb_add", "__add__", method_name_to_slot),
            BinopSlot(bf, "nb_subtract", "__sub__", method_name_to_slot),
            BinopSlot(bf, "nb_multiply", "__mul__", method_name_to_slot),
            BinopSlot(bf, "nb_remainder", "__mod__", method_name_to_slot),
            BinopSlot(bf, "nb_divmod", "__divmod__", method_name_to_slot),
            BinopSlot(ptf, "nb_power", "__pow__", method_name_to_slot),
            MethodSlot(unaryfunc, "nb_negative", "__neg__", method_name_to_slot),
            MethodSlot(unaryfunc, "nb_positive", "__pos__", method_name_to_slot),
            MethodSlot(unaryfunc, "nb_absolute", "__abs__", method_name_to_slot),
            MethodSlot(inquiry, "nb_bool", "__bool__", method_name_to_slot,
                       fallback="__nonzero__"),
            MethodSlot(unaryfunc, "nb_invert", "__invert__", method_name_to_slot),
            BinopSlot(bf, "nb_lshift", "__lshift__", method_name_to_slot),
            BinopSlot(bf, "nb_rshift", "__rshift__", method_name_to_slot),
            BinopSlot(bf, "nb_and", "__and__", method_name_to_slot),
            BinopSlot(bf, "nb_xor", "__xor__", method_name_to_slot),
            BinopSlot(bf, "nb_or", "__or__", method_name_to_slot),
            MethodSlot(unaryfunc, "nb_int", "__int__", method_name_to_slot, fallback="__long__"),
            EmptySlot("nb_long (reserved)"),
            MethodSlot(unaryfunc, "nb_float", "__float__", method_name_to_slot),

            # Added in release 2.0
            MethodSlot(ibinaryfunc, "nb_inplace_add", "__iadd__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_subtract", "__isub__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_multiply", "__imul__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_remainder", "__imod__", method_name_to_slot),
            MethodSlot(ptf, "nb_inplace_power", "__ipow__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_lshift", "__ilshift__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_rshift", "__irshift__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_and", "__iand__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_xor", "__ixor__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_or", "__ior__", method_name_to_slot),

            # Added in release 2.2
            # The following require the Py_TPFLAGS_HAVE_CLASS flag
            BinopSlot(bf, "nb_floor_divide", "__floordiv__", method_name_to_slot),
            BinopSlot(bf, "nb_true_divide", "__truediv__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_floor_divide", "__ifloordiv__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_true_divide", "__itruediv__", method_name_to_slot),

            # Added in release 2.5
            MethodSlot(unaryfunc, "nb_index", "__index__", method_name_to_slot),

            # Added in release 3.5
            BinopSlot(bf, "nb_matrix_multiply", "__matmul__", method_name_to_slot),
            MethodSlot(ibinaryfunc, "nb_inplace_matrix_multiply", "__imatmul__", method_name_to_slot),
        )

        self.PySequenceMethods = (
            MethodSlot(lenfunc, "sq_length", "__len__", method_name_to_slot),
            EmptySlot("sq_concat"),  # nb_add used instead
            EmptySlot("sq_repeat"),  # nb_multiply used instead
            SyntheticSlot("sq_item", ["__getitem__"], "0"),    #EmptySlot("sq_item"),   # mp_subscript used instead
            EmptySlot("sq_slice"),
            EmptySlot("sq_ass_item"),  # mp_ass_subscript used instead
            EmptySlot("sq_ass_slice"),
            MethodSlot(cmpfunc, "sq_contains", "__contains__", method_name_to_slot),
            EmptySlot("sq_inplace_concat"),  # nb_inplace_add used instead
            EmptySlot("sq_inplace_repeat"),  # nb_inplace_multiply used instead
        )

        self.PyMappingMethods = (
            MethodSlot(lenfunc, "mp_length", "__len__", method_name_to_slot),
            MethodSlot(objargfunc, "mp_subscript", "__getitem__", method_name_to_slot),
            SyntheticSlot("mp_ass_subscript", ["__setitem__", "__delitem__"], "0"),
        )

        self.PyBufferProcs = (
            MethodSlot(getbufferproc, "bf_getbuffer", "__getbuffer__", method_name_to_slot),
            MethodSlot(releasebufferproc, "bf_releasebuffer", "__releasebuffer__", method_name_to_slot)
        )

        self.PyAsyncMethods = (
            MethodSlot(unaryfunc, "am_await", "__await__", method_name_to_slot),
            MethodSlot(unaryfunc, "am_aiter", "__aiter__", method_name_to_slot),
            MethodSlot(unaryfunc, "am_anext", "__anext__", method_name_to_slot),
            # We should not map arbitrary .send() methods to an async slot.
            #MethodSlot(sendfunc, "am_send", "send", method_name_to_slot),
            EmptySlot("am_send"),
        )

        self.slot_table = (
            ConstructorSlot("tp_dealloc", '__dealloc__'),
            EmptySlot("tp_vectorcall_offset"),
            EmptySlot("tp_getattr"),
            EmptySlot("tp_setattr"),

            SuiteSlot(self. PyAsyncMethods, "__Pyx_PyAsyncMethodsStruct", "tp_as_async",
                      self.substructures, cast_cname="PyAsyncMethods"),

            MethodSlot(reprfunc, "tp_repr", "__repr__", method_name_to_slot),

            SuiteSlot(self.PyNumberMethods, "PyNumberMethods", "tp_as_number", self.substructures),
            SuiteSlot(self.PySequenceMethods, "PySequenceMethods", "tp_as_sequence", self.substructures),
            SuiteSlot(self.PyMappingMethods, "PyMappingMethods", "tp_as_mapping", self.substructures),

            MethodSlot(hashfunc, "tp_hash", "__hash__", method_name_to_slot,
                       inherited=False),    # Py3 checks for __richcmp__
            MethodSlot(callfunc, "tp_call", "__call__", method_name_to_slot),
            MethodSlot(reprfunc, "tp_str", "__str__", method_name_to_slot),

            SyntheticSlot("tp_getattro", ["__getattr__","__getattribute__"], "0"),  #"PyObject_GenericGetAttr"),
            SyntheticSlot("tp_setattro", ["__setattr__", "__delattr__"], "0"),  #"PyObject_GenericSetAttr"),

            SuiteSlot(self.PyBufferProcs, "PyBufferProcs", "tp_as_buffer", self.substructures),

            TypeFlagsSlot("tp_flags"),
            DocStringSlot("tp_doc"),

            GCDependentSlot("tp_traverse"),
            GCClearReferencesSlot("tp_clear"),

            RichcmpSlot(richcmpfunc, "tp_richcompare", "__richcmp__", method_name_to_slot,
                        inherited=False),  # Py3 checks for __hash__

            EmptySlot("tp_weaklistoffset"),

            MethodSlot(getiterfunc, "tp_iter", "__iter__", method_name_to_slot),
            MethodSlot(iternextfunc, "tp_iternext", "__next__", method_name_to_slot),

            MethodTableSlot("tp_methods"),
            MemberTableSlot("tp_members"),
            GetSetSlot("tp_getset"),

            BaseClassSlot("tp_base"),  #EmptySlot("tp_base"),
            EmptySlot("tp_dict"),

            SyntheticSlot("tp_descr_get", ["__get__"], "0"),
            SyntheticSlot("tp_descr_set", ["__set__", "__delete__"], "0"),

            DictOffsetSlot("tp_dictoffset", ifdef="!CYTHON_USE_TYPE_SPECS"),  # otherwise set via "__dictoffset__" member

            MethodSlot(initproc, "tp_init", "__init__", method_name_to_slot),
            EmptySlot("tp_alloc"),  #FixedSlot("tp_alloc", "PyType_GenericAlloc"),
            ConstructorSlot("tp_new", "__cinit__"),
            EmptySlot("tp_free"),

            EmptySlot("tp_is_gc"),
            EmptySlot("tp_bases"),
            EmptySlot("tp_mro"),
            EmptySlot("tp_cache"),
            EmptySlot("tp_subclasses"),
            EmptySlot("tp_weaklist"),
            EmptySlot("tp_del"),
            EmptySlot("tp_version_tag"),
            SyntheticSlot("tp_finalize", ["__del__"], "0",
                          used_ifdef="CYTHON_USE_TP_FINALIZE"),
            EmptySlot("tp_vectorcall", ifdef="!CYTHON_COMPILING_IN_PYPY || PYPY_VERSION_NUM >= 0x07030800"),
            EmptySlot("tp_print", ifdef="__PYX_NEED_TP_PRINT_SLOT == 1"),
            EmptySlot("tp_watched", ifdef="PY_VERSION_HEX >= 0x030C0000"),
            EmptySlot("tp_versions_used", ifdef="PY_VERSION_HEX >= 0x030d00A4"),
            # PyPy specific extension - only here to avoid C compiler warnings.
            EmptySlot("tp_pypy_flags", ifdef="CYTHON_COMPILING_IN_PYPY && PY_VERSION_HEX >= 0x03090000 && PY_VERSION_HEX < 0x030a0000"),
        )

        #------------------------------------------------------------------------------------------
        #
        #  Descriptors for special methods which don't appear directly
        #  in the type object or its substructures. These methods are
        #  called from slot functions synthesized by Cython.
        #
        #------------------------------------------------------------------------------------------

        MethodSlot(initproc, "", "__cinit__", method_name_to_slot)
        MethodSlot(destructor, "", "__dealloc__", method_name_to_slot)
        MethodSlot(destructor, "", "__del__", method_name_to_slot)
        MethodSlot(objobjargproc, "", "__setitem__", method_name_to_slot)
        MethodSlot(objargproc, "", "__delitem__", method_name_to_slot)
        MethodSlot(ssizessizeobjargproc, "", "__setslice__", method_name_to_slot)
        MethodSlot(ssizessizeargproc, "", "__delslice__", method_name_to_slot)
        MethodSlot(getattrofunc, "", "__getattr__", method_name_to_slot)
        MethodSlot(getattrofunc, "", "__getattribute__", method_name_to_slot)
        MethodSlot(setattrofunc, "", "__setattr__", method_name_to_slot)
        MethodSlot(delattrofunc, "", "__delattr__", method_name_to_slot)
        MethodSlot(descrgetfunc, "", "__get__", method_name_to_slot)
        MethodSlot(descrsetfunc, "", "__set__", method_name_to_slot)
        MethodSlot(descrdelfunc, "", "__delete__", method_name_to_slot)

        #-------------------------------------------------------------------------
        #
        # Legacy "fallback" Py2 slots. Don't appear in the generated slot table,
        # but match the "fallback" argument of a slot that does
        #
        #-------------------------------------------------------------------------
        MethodSlot(inquiry, "", "__nonzero__", method_name_to_slot)
        MethodSlot(unaryfunc, "", "__long__", method_name_to_slot)

    def get_special_method_signature(self, name):
        #  Given a method name, if it is a special method,
        #  return its signature, else return None.
        slot = self._get_slot_by_method_name(name)
        if slot:
            return slot.signature
        elif name in richcmp_special_methods:
            return ibinaryfunc
        else:
            return None

    def get_slot_by_method_name(self, method_name):
        # For now, only search the type struct, no referenced sub-structs.
        return self._get_slot_by_method_name(method_name)

    def __iter__(self):
        # make it easier to iterate over all the slots
        return iter(self.slot_table)


_slot_table_dict = {}

def get_slot_table(compiler_directives):
    if not compiler_directives:
        # fetch default directives here since the builtin type classes don't have
        # directives set
        from .Options import get_directive_defaults
        compiler_directives = get_directive_defaults()

    old_binops = compiler_directives['c_api_binop_methods']
    key = (old_binops,)
    if key not in _slot_table_dict:
        _slot_table_dict[key] = SlotTable(old_binops=old_binops)
    return _slot_table_dict[key]


# Populate "special_method_names" based on the default directives (so it can always be accessed quickly).
special_method_names = set(get_slot_table(compiler_directives=None))


# Method flags for python-exposed methods.

method_noargs   = "METH_NOARGS"
method_onearg   = "METH_O"
method_varargs  = "METH_VARARGS"
method_fastcall = "__Pyx_METH_FASTCALL"  # Actually VARARGS on versions < 3.7
method_keywords = "METH_KEYWORDS"
method_coexist  = "METH_COEXIST"
