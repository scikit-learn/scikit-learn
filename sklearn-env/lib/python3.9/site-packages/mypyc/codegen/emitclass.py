"""Code generation for native classes and related wrappers."""

from typing import Optional, List, Tuple, Dict, Callable, Mapping, Set

from mypy.backports import OrderedDict

from mypyc.common import PREFIX, NATIVE_PREFIX, REG_PREFIX, use_fastcall
from mypyc.codegen.emit import Emitter, HeaderDeclaration, ReturnHandler
from mypyc.codegen.emitfunc import native_function_header
from mypyc.codegen.emitwrapper import (
    generate_dunder_wrapper, generate_hash_wrapper, generate_richcompare_wrapper,
    generate_bool_wrapper, generate_get_wrapper, generate_len_wrapper,
    generate_set_del_item_wrapper, generate_contains_wrapper, generate_bin_op_wrapper
)
from mypyc.ir.rtypes import RType, RTuple, object_rprimitive
from mypyc.ir.func_ir import FuncIR, FuncDecl, FUNC_STATICMETHOD, FUNC_CLASSMETHOD
from mypyc.ir.class_ir import ClassIR, VTableEntries
from mypyc.sametype import is_same_type
from mypyc.namegen import NameGenerator


def native_slot(cl: ClassIR, fn: FuncIR, emitter: Emitter) -> str:
    return '{}{}'.format(NATIVE_PREFIX, fn.cname(emitter.names))


def wrapper_slot(cl: ClassIR, fn: FuncIR, emitter: Emitter) -> str:
    return '{}{}'.format(PREFIX, fn.cname(emitter.names))


# We maintain a table from dunder function names to struct slots they
# correspond to and functions that generate a wrapper (if necessary)
# and return the function name to stick in the slot.
# TODO: Add remaining dunder methods
SlotGenerator = Callable[[ClassIR, FuncIR, Emitter], str]
SlotTable = Mapping[str, Tuple[str, SlotGenerator]]

SLOT_DEFS: SlotTable = {
    '__init__': ('tp_init', lambda c, t, e: generate_init_for_class(c, t, e)),
    '__call__': ('tp_call', lambda c, t, e: generate_call_wrapper(c, t, e)),
    '__str__': ('tp_str', native_slot),
    '__repr__': ('tp_repr', native_slot),
    '__next__': ('tp_iternext', native_slot),
    '__iter__': ('tp_iter', native_slot),
    '__hash__': ('tp_hash', generate_hash_wrapper),
    '__get__': ('tp_descr_get', generate_get_wrapper),
}

AS_MAPPING_SLOT_DEFS: SlotTable = {
    '__getitem__': ('mp_subscript', generate_dunder_wrapper),
    '__setitem__': ('mp_ass_subscript', generate_set_del_item_wrapper),
    '__delitem__': ('mp_ass_subscript', generate_set_del_item_wrapper),
    '__len__': ('mp_length', generate_len_wrapper),
}

AS_SEQUENCE_SLOT_DEFS: SlotTable = {
    '__contains__': ('sq_contains', generate_contains_wrapper),
}

AS_NUMBER_SLOT_DEFS: SlotTable = {
    '__bool__': ('nb_bool', generate_bool_wrapper),
    '__neg__': ('nb_negative', generate_dunder_wrapper),
    '__invert__': ('nb_invert', generate_dunder_wrapper),
    '__int__': ('nb_int', generate_dunder_wrapper),
    '__float__': ('nb_float', generate_dunder_wrapper),
    '__add__': ('nb_add', generate_bin_op_wrapper),
    '__radd__': ('nb_add', generate_bin_op_wrapper),
    '__sub__': ('nb_subtract', generate_bin_op_wrapper),
    '__rsub__': ('nb_subtract', generate_bin_op_wrapper),
    '__mul__': ('nb_multiply', generate_bin_op_wrapper),
    '__rmul__': ('nb_multiply', generate_bin_op_wrapper),
    '__mod__': ('nb_remainder', generate_bin_op_wrapper),
    '__rmod__': ('nb_remainder', generate_bin_op_wrapper),
    '__truediv__': ('nb_true_divide', generate_bin_op_wrapper),
    '__rtruediv__': ('nb_true_divide', generate_bin_op_wrapper),
    '__floordiv__': ('nb_floor_divide', generate_bin_op_wrapper),
    '__rfloordiv__': ('nb_floor_divide', generate_bin_op_wrapper),
    '__lshift__': ('nb_lshift', generate_bin_op_wrapper),
    '__rlshift__': ('nb_lshift', generate_bin_op_wrapper),
    '__rshift__': ('nb_rshift', generate_bin_op_wrapper),
    '__rrshift__': ('nb_rshift', generate_bin_op_wrapper),
    '__and__': ('nb_and', generate_bin_op_wrapper),
    '__rand__': ('nb_and', generate_bin_op_wrapper),
    '__or__': ('nb_or', generate_bin_op_wrapper),
    '__ror__': ('nb_or', generate_bin_op_wrapper),
    '__xor__': ('nb_xor', generate_bin_op_wrapper),
    '__rxor__': ('nb_xor', generate_bin_op_wrapper),
    '__matmul__': ('nb_matrix_multiply', generate_bin_op_wrapper),
    '__rmatmul__': ('nb_matrix_multiply', generate_bin_op_wrapper),
    '__iadd__': ('nb_inplace_add', generate_dunder_wrapper),
    '__isub__': ('nb_inplace_subtract', generate_dunder_wrapper),
    '__imul__': ('nb_inplace_multiply', generate_dunder_wrapper),
    '__imod__': ('nb_inplace_remainder', generate_dunder_wrapper),
    '__itruediv__': ('nb_inplace_true_divide', generate_dunder_wrapper),
    '__ifloordiv__': ('nb_inplace_floor_divide', generate_dunder_wrapper),
    '__ilshift__': ('nb_inplace_lshift', generate_dunder_wrapper),
    '__irshift__': ('nb_inplace_rshift', generate_dunder_wrapper),
    '__iand__': ('nb_inplace_and', generate_dunder_wrapper),
    '__ior__': ('nb_inplace_or', generate_dunder_wrapper),
    '__ixor__': ('nb_inplace_xor', generate_dunder_wrapper),
    '__imatmul__': ('nb_inplace_matrix_multiply', generate_dunder_wrapper),
}

AS_ASYNC_SLOT_DEFS: SlotTable = {
    '__await__': ('am_await', native_slot),
    '__aiter__': ('am_aiter', native_slot),
    '__anext__': ('am_anext', native_slot),
}

SIDE_TABLES = [
    ('as_mapping', 'PyMappingMethods', AS_MAPPING_SLOT_DEFS),
    ('as_sequence', 'PySequenceMethods', AS_SEQUENCE_SLOT_DEFS),
    ('as_number', 'PyNumberMethods', AS_NUMBER_SLOT_DEFS),
    ('as_async', 'PyAsyncMethods', AS_ASYNC_SLOT_DEFS),
]

# Slots that need to always be filled in because they don't get
# inherited right.
ALWAYS_FILL = {
    '__hash__',
}


def generate_call_wrapper(cl: ClassIR, fn: FuncIR, emitter: Emitter) -> str:
    if emitter.use_vectorcall():
        # Use vectorcall wrapper if supported (PEP 590).
        return 'PyVectorcall_Call'
    else:
        # On older Pythons use the legacy wrapper.
        return wrapper_slot(cl, fn, emitter)


def slot_key(attr: str) -> str:
    """Map dunder method name to sort key.

    Sort reverse operator methods and __delitem__ after others ('x' > '_').
    """
    if (attr.startswith('__r') and attr != '__rshift__') or attr == '__delitem__':
        return 'x' + attr
    return attr


def generate_slots(cl: ClassIR, table: SlotTable, emitter: Emitter) -> Dict[str, str]:
    fields: Dict[str, str] = OrderedDict()
    generated: Dict[str, str] = {}
    # Sort for determinism on Python 3.5
    for name, (slot, generator) in sorted(table.items(), key=lambda x: slot_key(x[0])):
        method_cls = cl.get_method_and_class(name)
        if method_cls and (method_cls[1] == cl or name in ALWAYS_FILL):
            if slot in generated:
                # Reuse previously generated wrapper.
                fields[slot] = generated[slot]
            else:
                # Generate new wrapper.
                name = generator(cl, method_cls[0], emitter)
                fields[slot] = name
                generated[slot] = name

    return fields


def generate_class_type_decl(cl: ClassIR, c_emitter: Emitter,
                             external_emitter: Emitter,
                             emitter: Emitter) -> None:
    context = c_emitter.context
    name = emitter.type_struct_name(cl)
    context.declarations[name] = HeaderDeclaration(
        'PyTypeObject *{};'.format(emitter.type_struct_name(cl)),
        needs_export=True)

    # If this is a non-extension class, all we want is the type object decl.
    if not cl.is_ext_class:
        return

    generate_object_struct(cl, external_emitter)
    generate_full = not cl.is_trait and not cl.builtin_base
    if generate_full:
        context.declarations[emitter.native_function_name(cl.ctor)] = HeaderDeclaration(
            '{};'.format(native_function_header(cl.ctor, emitter)),
            needs_export=True,
        )


def generate_class(cl: ClassIR, module: str, emitter: Emitter) -> None:
    """Generate C code for a class.

    This is the main entry point to the module.
    """
    name = cl.name
    name_prefix = cl.name_prefix(emitter.names)

    setup_name = '{}_setup'.format(name_prefix)
    new_name = '{}_new'.format(name_prefix)
    members_name = '{}_members'.format(name_prefix)
    getseters_name = '{}_getseters'.format(name_prefix)
    vtable_name = '{}_vtable'.format(name_prefix)
    traverse_name = '{}_traverse'.format(name_prefix)
    clear_name = '{}_clear'.format(name_prefix)
    dealloc_name = '{}_dealloc'.format(name_prefix)
    methods_name = '{}_methods'.format(name_prefix)
    vtable_setup_name = '{}_trait_vtable_setup'.format(name_prefix)

    fields: Dict[str, str] = OrderedDict()
    fields['tp_name'] = '"{}"'.format(name)

    generate_full = not cl.is_trait and not cl.builtin_base
    needs_getseters = cl.needs_getseters or not cl.is_generated

    if not cl.builtin_base:
        fields['tp_new'] = new_name

    if generate_full:
        fields['tp_dealloc'] = '(destructor){}_dealloc'.format(name_prefix)
        fields['tp_traverse'] = '(traverseproc){}_traverse'.format(name_prefix)
        fields['tp_clear'] = '(inquiry){}_clear'.format(name_prefix)
    if needs_getseters:
        fields['tp_getset'] = getseters_name
    fields['tp_methods'] = methods_name

    def emit_line() -> None:
        emitter.emit_line()

    emit_line()

    # If the class has a method to initialize default attribute
    # values, we need to call it during initialization.
    defaults_fn = cl.get_method('__mypyc_defaults_setup')

    # If there is a __init__ method, we'll use it in the native constructor.
    init_fn = cl.get_method('__init__')

    # Fill out slots in the type object from dunder methods.
    fields.update(generate_slots(cl, SLOT_DEFS, emitter))

    # Fill out dunder methods that live in tables hanging off the side.
    for table_name, type, slot_defs in SIDE_TABLES:
        slots = generate_slots(cl, slot_defs, emitter)
        if slots:
            table_struct_name = generate_side_table_for_class(cl, table_name, type, slots, emitter)
            fields['tp_{}'.format(table_name)] = '&{}'.format(table_struct_name)

    richcompare_name = generate_richcompare_wrapper(cl, emitter)
    if richcompare_name:
        fields['tp_richcompare'] = richcompare_name

    # If the class inherits from python, make space for a __dict__
    struct_name = cl.struct_name(emitter.names)
    if cl.builtin_base:
        base_size = 'sizeof({})'.format(cl.builtin_base)
    elif cl.is_trait:
        base_size = 'sizeof(PyObject)'
    else:
        base_size = 'sizeof({})'.format(struct_name)
    # Since our types aren't allocated using type() we need to
    # populate these fields ourselves if we want them to have correct
    # values. PyType_Ready will inherit the offsets from tp_base but
    # that isn't what we want.

    # XXX: there is no reason for the __weakref__ stuff to be mixed up with __dict__
    if cl.has_dict:
        # __dict__ lives right after the struct and __weakref__ lives right after that
        # TODO: They should get members in the struct instead of doing this nonsense.
        weak_offset = '{} + sizeof(PyObject *)'.format(base_size)
        emitter.emit_lines(
            'PyMemberDef {}[] = {{'.format(members_name),
            '{{"__dict__", T_OBJECT_EX, {}, 0, NULL}},'.format(base_size),
            '{{"__weakref__", T_OBJECT_EX, {}, 0, NULL}},'.format(weak_offset),
            '{0}',
            '};',
        )

        fields['tp_members'] = members_name
        fields['tp_basicsize'] = '{} + 2*sizeof(PyObject *)'.format(base_size)
        fields['tp_dictoffset'] = base_size
        fields['tp_weaklistoffset'] = weak_offset
    else:
        fields['tp_basicsize'] = base_size

    if generate_full:
        # Declare setup method that allocates and initializes an object. type is the
        # type of the class being initialized, which could be another class if there
        # is an interpreted subclass.
        emitter.emit_line('static PyObject *{}(PyTypeObject *type);'.format(setup_name))
        assert cl.ctor is not None
        emitter.emit_line(native_function_header(cl.ctor, emitter) + ';')

        emit_line()
        generate_new_for_class(cl, new_name, vtable_name, setup_name, emitter)
        emit_line()
        generate_traverse_for_class(cl, traverse_name, emitter)
        emit_line()
        generate_clear_for_class(cl, clear_name, emitter)
        emit_line()
        generate_dealloc_for_class(cl, dealloc_name, clear_name, emitter)
        emit_line()

        if cl.allow_interpreted_subclasses:
            shadow_vtable_name: Optional[str] = generate_vtables(
                cl, vtable_setup_name + "_shadow", vtable_name + "_shadow", emitter, shadow=True
            )
            emit_line()
        else:
            shadow_vtable_name = None
        vtable_name = generate_vtables(cl, vtable_setup_name, vtable_name, emitter, shadow=False)
        emit_line()
    if needs_getseters:
        generate_getseter_declarations(cl, emitter)
        emit_line()
        generate_getseters_table(cl, getseters_name, emitter)
        emit_line()

    if cl.is_trait:
        generate_new_for_trait(cl, new_name, emitter)

    generate_methods_table(cl, methods_name, emitter)
    emit_line()

    flags = ['Py_TPFLAGS_DEFAULT', 'Py_TPFLAGS_HEAPTYPE', 'Py_TPFLAGS_BASETYPE']
    if generate_full:
        flags.append('Py_TPFLAGS_HAVE_GC')
    if cl.has_method('__call__') and emitter.use_vectorcall():
        fields['tp_vectorcall_offset'] = 'offsetof({}, vectorcall)'.format(
            cl.struct_name(emitter.names))
        flags.append('_Py_TPFLAGS_HAVE_VECTORCALL')
    fields['tp_flags'] = ' | '.join(flags)

    emitter.emit_line("static PyTypeObject {}_template_ = {{".format(emitter.type_struct_name(cl)))
    emitter.emit_line("PyVarObject_HEAD_INIT(NULL, 0)")
    for field, value in fields.items():
        emitter.emit_line(".{} = {},".format(field, value))
    emitter.emit_line("};")
    emitter.emit_line("static PyTypeObject *{t}_template = &{t}_template_;".format(
        t=emitter.type_struct_name(cl)))

    emitter.emit_line()
    if generate_full:
        generate_setup_for_class(
            cl, setup_name, defaults_fn, vtable_name, shadow_vtable_name, emitter)
        emitter.emit_line()
        generate_constructor_for_class(
            cl, cl.ctor, init_fn, setup_name, vtable_name, emitter)
        emitter.emit_line()
    if needs_getseters:
        generate_getseters(cl, emitter)


def getter_name(cl: ClassIR, attribute: str, names: NameGenerator) -> str:
    return names.private_name(cl.module_name, '{}_get{}'.format(cl.name, attribute))


def setter_name(cl: ClassIR, attribute: str, names: NameGenerator) -> str:
    return names.private_name(cl.module_name, '{}_set{}'.format(cl.name, attribute))


def generate_object_struct(cl: ClassIR, emitter: Emitter) -> None:
    seen_attrs: Set[Tuple[str, RType]] = set()
    lines: List[str] = []
    lines += ['typedef struct {',
              'PyObject_HEAD',
              'CPyVTableItem *vtable;']
    if cl.has_method('__call__') and emitter.use_vectorcall():
        lines.append('vectorcallfunc vectorcall;')
    for base in reversed(cl.base_mro):
        if not base.is_trait:
            for attr, rtype in base.attributes.items():
                if (attr, rtype) not in seen_attrs:
                    lines.append('{}{};'.format(emitter.ctype_spaced(rtype),
                                                emitter.attr(attr)))
                    seen_attrs.add((attr, rtype))

                    if isinstance(rtype, RTuple):
                        emitter.declare_tuple_struct(rtype)

    lines.append('}} {};'.format(cl.struct_name(emitter.names)))
    lines.append('')
    emitter.context.declarations[cl.struct_name(emitter.names)] = HeaderDeclaration(
        lines,
        is_type=True
    )


def generate_vtables(base: ClassIR,
                     vtable_setup_name: str,
                     vtable_name: str,
                     emitter: Emitter,
                     shadow: bool) -> str:
    """Emit the vtables and vtable setup functions for a class.

    This includes both the primary vtable and any trait implementation vtables.
    The trait vtables go before the main vtable, and have the following layout:
        {
            CPyType_T1,         // pointer to type object
            C_T1_trait_vtable,  // pointer to array of method pointers
            C_T1_offset_table,  // pointer to array of attribute offsets
            CPyType_T2,
            C_T2_trait_vtable,
            C_T2_offset_table,
            ...
        }
    The method implementations are calculated at the end of IR pass, attribute
    offsets are {offsetof(native__C, _x1), offsetof(native__C, _y1), ...}.

    To account for both dynamic loading and dynamic class creation,
    vtables are populated dynamically at class creation time, so we
    emit empty array definitions to store the vtables and a function to
    populate them.

    If shadow is True, generate "shadow vtables" that point to the
    shadow glue methods (which should dispatch via the Python C-API).

    Returns the expression to use to refer to the vtable, which might be
    different than the name, if there are trait vtables.
    """

    def trait_vtable_name(trait: ClassIR) -> str:
        return '{}_{}_trait_vtable{}'.format(
            base.name_prefix(emitter.names), trait.name_prefix(emitter.names),
            '_shadow' if shadow else '')

    def trait_offset_table_name(trait: ClassIR) -> str:
        return '{}_{}_offset_table'.format(
            base.name_prefix(emitter.names), trait.name_prefix(emitter.names)
        )

    # Emit array definitions with enough space for all the entries
    emitter.emit_line('static CPyVTableItem {}[{}];'.format(
        vtable_name,
        max(1, len(base.vtable_entries) + 3 * len(base.trait_vtables))))

    for trait, vtable in base.trait_vtables.items():
        # Trait methods entry (vtable index -> method implementation).
        emitter.emit_line('static CPyVTableItem {}[{}];'.format(
            trait_vtable_name(trait),
            max(1, len(vtable))))
        # Trait attributes entry (attribute number in trait -> offset in actual struct).
        emitter.emit_line('static size_t {}[{}];'.format(
            trait_offset_table_name(trait),
            max(1, len(trait.attributes)))
        )

    # Emit vtable setup function
    emitter.emit_line('static bool')
    emitter.emit_line('{}{}(void)'.format(NATIVE_PREFIX, vtable_setup_name))
    emitter.emit_line('{')

    if base.allow_interpreted_subclasses and not shadow:
        emitter.emit_line('{}{}_shadow();'.format(NATIVE_PREFIX, vtable_setup_name))

    subtables = []
    for trait, vtable in base.trait_vtables.items():
        name = trait_vtable_name(trait)
        offset_name = trait_offset_table_name(trait)
        generate_vtable(vtable, name, emitter, [], shadow)
        generate_offset_table(offset_name, emitter, trait, base)
        subtables.append((trait, name, offset_name))

    generate_vtable(base.vtable_entries, vtable_name, emitter, subtables, shadow)

    emitter.emit_line('return 1;')
    emitter.emit_line('}')

    return vtable_name if not subtables else "{} + {}".format(vtable_name, len(subtables) * 3)


def generate_offset_table(trait_offset_table_name: str,
                          emitter: Emitter,
                          trait: ClassIR,
                          cl: ClassIR) -> None:
    """Generate attribute offset row of a trait vtable."""
    emitter.emit_line('size_t {}_scratch[] = {{'.format(trait_offset_table_name))
    for attr in trait.attributes:
        emitter.emit_line('offsetof({}, {}),'.format(
            cl.struct_name(emitter.names), emitter.attr(attr)
        ))
    if not trait.attributes:
        # This is for msvc.
        emitter.emit_line('0')
    emitter.emit_line('};')
    emitter.emit_line('memcpy({name}, {name}_scratch, sizeof({name}));'.format(
        name=trait_offset_table_name)
    )


def generate_vtable(entries: VTableEntries,
                    vtable_name: str,
                    emitter: Emitter,
                    subtables: List[Tuple[ClassIR, str, str]],
                    shadow: bool) -> None:
    emitter.emit_line('CPyVTableItem {}_scratch[] = {{'.format(vtable_name))
    if subtables:
        emitter.emit_line('/* Array of trait vtables */')
        for trait, table, offset_table in subtables:
            emitter.emit_line(
                '(CPyVTableItem){}, (CPyVTableItem){}, (CPyVTableItem){},'.format(
                    emitter.type_struct_name(trait), table, offset_table))
        emitter.emit_line('/* Start of real vtable */')

    for entry in entries:
        method = entry.shadow_method if shadow and entry.shadow_method else entry.method
        emitter.emit_line('(CPyVTableItem){}{}{},'.format(
            emitter.get_group_prefix(entry.method.decl),
            NATIVE_PREFIX,
            method.cname(emitter.names)))

    # msvc doesn't allow empty arrays; maybe allowing them at all is an extension?
    if not entries:
        emitter.emit_line('NULL')
    emitter.emit_line('};')
    emitter.emit_line('memcpy({name}, {name}_scratch, sizeof({name}));'.format(name=vtable_name))


def generate_setup_for_class(cl: ClassIR,
                             func_name: str,
                             defaults_fn: Optional[FuncIR],
                             vtable_name: str,
                             shadow_vtable_name: Optional[str],
                             emitter: Emitter) -> None:
    """Generate a native function that allocates an instance of a class."""
    emitter.emit_line('static PyObject *')
    emitter.emit_line('{}(PyTypeObject *type)'.format(func_name))
    emitter.emit_line('{')
    emitter.emit_line('{} *self;'.format(cl.struct_name(emitter.names)))
    emitter.emit_line('self = ({struct} *)type->tp_alloc(type, 0);'.format(
        struct=cl.struct_name(emitter.names)))
    emitter.emit_line('if (self == NULL)')
    emitter.emit_line('    return NULL;')

    if shadow_vtable_name:
        emitter.emit_line('if (type != {}) {{'.format(emitter.type_struct_name(cl)))
        emitter.emit_line('self->vtable = {};'.format(shadow_vtable_name))
        emitter.emit_line('} else {')
        emitter.emit_line('self->vtable = {};'.format(vtable_name))
        emitter.emit_line('}')
    else:
        emitter.emit_line('self->vtable = {};'.format(vtable_name))

    if cl.has_method('__call__') and emitter.use_vectorcall():
        name = cl.method_decl('__call__').cname(emitter.names)
        emitter.emit_line('self->vectorcall = {}{};'.format(PREFIX, name))

    for base in reversed(cl.base_mro):
        for attr, rtype in base.attributes.items():
            emitter.emit_line('self->{} = {};'.format(
                emitter.attr(attr), emitter.c_undefined_value(rtype)))

    # Initialize attributes to default values, if necessary
    if defaults_fn is not None:
        emitter.emit_lines(
            'if ({}{}((PyObject *)self) == 0) {{'.format(
                NATIVE_PREFIX, defaults_fn.cname(emitter.names)),
            'Py_DECREF(self);',
            'return NULL;',
            '}')

    emitter.emit_line('return (PyObject *)self;')
    emitter.emit_line('}')


def generate_constructor_for_class(cl: ClassIR,
                                   fn: FuncDecl,
                                   init_fn: Optional[FuncIR],
                                   setup_name: str,
                                   vtable_name: str,
                                   emitter: Emitter) -> None:
    """Generate a native function that allocates and initializes an instance of a class."""
    emitter.emit_line('{}'.format(native_function_header(fn, emitter)))
    emitter.emit_line('{')
    emitter.emit_line('PyObject *self = {}({});'.format(setup_name, emitter.type_struct_name(cl)))
    emitter.emit_line('if (self == NULL)')
    emitter.emit_line('    return NULL;')
    args = ', '.join(['self'] + [REG_PREFIX + arg.name for arg in fn.sig.args])
    if init_fn is not None:
        emitter.emit_line('char res = {}{}{}({});'.format(
            emitter.get_group_prefix(init_fn.decl),
            NATIVE_PREFIX, init_fn.cname(emitter.names), args))
        emitter.emit_line('if (res == 2) {')
        emitter.emit_line('Py_DECREF(self);')
        emitter.emit_line('return NULL;')
        emitter.emit_line('}')

    # If there is a nontrivial ctor that we didn't define, invoke it via tp_init
    elif len(fn.sig.args) > 1:
        emitter.emit_line(
            'int res = {}->tp_init({});'.format(
                emitter.type_struct_name(cl),
                args))

        emitter.emit_line('if (res < 0) {')
        emitter.emit_line('Py_DECREF(self);')
        emitter.emit_line('return NULL;')
        emitter.emit_line('}')

    emitter.emit_line('return self;')
    emitter.emit_line('}')


def generate_init_for_class(cl: ClassIR,
                            init_fn: FuncIR,
                            emitter: Emitter) -> str:
    """Generate an init function suitable for use as tp_init.

    tp_init needs to be a function that returns an int, and our
    __init__ methods return a PyObject. Translate NULL to -1,
    everything else to 0.
    """
    func_name = '{}_init'.format(cl.name_prefix(emitter.names))

    emitter.emit_line('static int')
    emitter.emit_line(
        '{}(PyObject *self, PyObject *args, PyObject *kwds)'.format(func_name))
    emitter.emit_line('{')
    emitter.emit_line('return {}{}(self, args, kwds) != NULL ? 0 : -1;'.format(
        PREFIX, init_fn.cname(emitter.names)))
    emitter.emit_line('}')

    return func_name


def generate_new_for_class(cl: ClassIR,
                           func_name: str,
                           vtable_name: str,
                           setup_name: str,
                           emitter: Emitter) -> None:
    emitter.emit_line('static PyObject *')
    emitter.emit_line(
        '{}(PyTypeObject *type, PyObject *args, PyObject *kwds)'.format(func_name))
    emitter.emit_line('{')
    # TODO: Check and unbox arguments
    if not cl.allow_interpreted_subclasses:
        emitter.emit_line('if (type != {}) {{'.format(emitter.type_struct_name(cl)))
        emitter.emit_line(
            'PyErr_SetString(PyExc_TypeError, "interpreted classes cannot inherit from compiled");'
        )
        emitter.emit_line('return NULL;')
        emitter.emit_line('}')

    emitter.emit_line('return {}(type);'.format(setup_name))
    emitter.emit_line('}')


def generate_new_for_trait(cl: ClassIR,
                           func_name: str,
                           emitter: Emitter) -> None:
    emitter.emit_line('static PyObject *')
    emitter.emit_line(
        '{}(PyTypeObject *type, PyObject *args, PyObject *kwds)'.format(func_name))
    emitter.emit_line('{')
    emitter.emit_line('if (type != {}) {{'.format(emitter.type_struct_name(cl)))
    emitter.emit_line(
        'PyErr_SetString(PyExc_TypeError, '
        '"interpreted classes cannot inherit from compiled traits");'
    )
    emitter.emit_line('} else {')
    emitter.emit_line(
        'PyErr_SetString(PyExc_TypeError, "traits may not be directly created");'
    )
    emitter.emit_line('}')
    emitter.emit_line('return NULL;')
    emitter.emit_line('}')


def generate_traverse_for_class(cl: ClassIR,
                                func_name: str,
                                emitter: Emitter) -> None:
    """Emit function that performs cycle GC traversal of an instance."""
    emitter.emit_line('static int')
    emitter.emit_line('{}({} *self, visitproc visit, void *arg)'.format(
        func_name,
        cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    for base in reversed(cl.base_mro):
        for attr, rtype in base.attributes.items():
            emitter.emit_gc_visit('self->{}'.format(emitter.attr(attr)), rtype)
    if cl.has_dict:
        struct_name = cl.struct_name(emitter.names)
        # __dict__ lives right after the struct and __weakref__ lives right after that
        emitter.emit_gc_visit('*((PyObject **)((char *)self + sizeof({})))'.format(
            struct_name), object_rprimitive)
        emitter.emit_gc_visit(
            '*((PyObject **)((char *)self + sizeof(PyObject *) + sizeof({})))'.format(
                struct_name),
            object_rprimitive)
    emitter.emit_line('return 0;')
    emitter.emit_line('}')


def generate_clear_for_class(cl: ClassIR,
                             func_name: str,
                             emitter: Emitter) -> None:
    emitter.emit_line('static int')
    emitter.emit_line('{}({} *self)'.format(func_name, cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    for base in reversed(cl.base_mro):
        for attr, rtype in base.attributes.items():
            emitter.emit_gc_clear('self->{}'.format(emitter.attr(attr)), rtype)
    if cl.has_dict:
        struct_name = cl.struct_name(emitter.names)
        # __dict__ lives right after the struct and __weakref__ lives right after that
        emitter.emit_gc_clear('*((PyObject **)((char *)self + sizeof({})))'.format(
            struct_name), object_rprimitive)
        emitter.emit_gc_clear(
            '*((PyObject **)((char *)self + sizeof(PyObject *) + sizeof({})))'.format(
                struct_name),
            object_rprimitive)
    emitter.emit_line('return 0;')
    emitter.emit_line('}')


def generate_dealloc_for_class(cl: ClassIR,
                               dealloc_func_name: str,
                               clear_func_name: str,
                               emitter: Emitter) -> None:
    emitter.emit_line('static void')
    emitter.emit_line('{}({} *self)'.format(dealloc_func_name, cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    emitter.emit_line('PyObject_GC_UnTrack(self);')
    # The trashcan is needed to handle deep recursive deallocations
    emitter.emit_line('CPy_TRASHCAN_BEGIN(self, {})'.format(dealloc_func_name))
    emitter.emit_line('{}(self);'.format(clear_func_name))
    emitter.emit_line('Py_TYPE(self)->tp_free((PyObject *)self);')
    emitter.emit_line('CPy_TRASHCAN_END(self)')
    emitter.emit_line('}')


def generate_methods_table(cl: ClassIR,
                           name: str,
                           emitter: Emitter) -> None:
    emitter.emit_line('static PyMethodDef {}[] = {{'.format(name))
    for fn in cl.methods.values():
        if fn.decl.is_prop_setter or fn.decl.is_prop_getter:
            continue
        emitter.emit_line('{{"{}",'.format(fn.name))
        emitter.emit_line(' (PyCFunction){}{},'.format(PREFIX, fn.cname(emitter.names)))
        if use_fastcall(emitter.capi_version):
            flags = ['METH_FASTCALL']
        else:
            flags = ['METH_VARARGS']
        flags.append('METH_KEYWORDS')
        if fn.decl.kind == FUNC_STATICMETHOD:
            flags.append('METH_STATIC')
        elif fn.decl.kind == FUNC_CLASSMETHOD:
            flags.append('METH_CLASS')

        emitter.emit_line(' {}, NULL}},'.format(' | '.join(flags)))

    # Provide a default __getstate__ and __setstate__
    if not cl.has_method('__setstate__') and not cl.has_method('__getstate__'):
        emitter.emit_lines(
            '{"__setstate__", (PyCFunction)CPyPickle_SetState, METH_O, NULL},',
            '{"__getstate__", (PyCFunction)CPyPickle_GetState, METH_NOARGS, NULL},',
        )

    emitter.emit_line('{NULL}  /* Sentinel */')
    emitter.emit_line('};')


def generate_side_table_for_class(cl: ClassIR,
                                  name: str,
                                  type: str,
                                  slots: Dict[str, str],
                                  emitter: Emitter) -> Optional[str]:
    name = '{}_{}'.format(cl.name_prefix(emitter.names), name)
    emitter.emit_line('static {} {} = {{'.format(type, name))
    for field, value in slots.items():
        emitter.emit_line(".{} = {},".format(field, value))
    emitter.emit_line("};")
    return name


def generate_getseter_declarations(cl: ClassIR, emitter: Emitter) -> None:
    if not cl.is_trait:
        for attr in cl.attributes:
            emitter.emit_line('static PyObject *')
            emitter.emit_line('{}({} *self, void *closure);'.format(
                getter_name(cl, attr, emitter.names),
                cl.struct_name(emitter.names)))
            emitter.emit_line('static int')
            emitter.emit_line('{}({} *self, PyObject *value, void *closure);'.format(
                setter_name(cl, attr, emitter.names),
                cl.struct_name(emitter.names)))

    for prop in cl.properties:
        # Generate getter declaration
        emitter.emit_line('static PyObject *')
        emitter.emit_line('{}({} *self, void *closure);'.format(
            getter_name(cl, prop, emitter.names),
            cl.struct_name(emitter.names)))

        # Generate property setter declaration if a setter exists
        if cl.properties[prop][1]:
            emitter.emit_line('static int')
            emitter.emit_line('{}({} *self, PyObject *value, void *closure);'.format(
                setter_name(cl, prop, emitter.names),
                cl.struct_name(emitter.names)))


def generate_getseters_table(cl: ClassIR,
                             name: str,
                             emitter: Emitter) -> None:
    emitter.emit_line('static PyGetSetDef {}[] = {{'.format(name))
    if not cl.is_trait:
        for attr in cl.attributes:
            emitter.emit_line('{{"{}",'.format(attr))
            emitter.emit_line(' (getter){}, (setter){},'.format(
                getter_name(cl, attr, emitter.names), setter_name(cl, attr, emitter.names)))
            emitter.emit_line(' NULL, NULL},')
    for prop in cl.properties:
        emitter.emit_line('{{"{}",'.format(prop))
        emitter.emit_line(' (getter){},'.format(getter_name(cl, prop, emitter.names)))

        setter = cl.properties[prop][1]
        if setter:
            emitter.emit_line(' (setter){},'.format(setter_name(cl, prop, emitter.names)))
            emitter.emit_line('NULL, NULL},')
        else:
            emitter.emit_line('NULL, NULL, NULL},')

    emitter.emit_line('{NULL}  /* Sentinel */')
    emitter.emit_line('};')


def generate_getseters(cl: ClassIR, emitter: Emitter) -> None:
    if not cl.is_trait:
        for i, (attr, rtype) in enumerate(cl.attributes.items()):
            generate_getter(cl, attr, rtype, emitter)
            emitter.emit_line('')
            generate_setter(cl, attr, rtype, emitter)
            if i < len(cl.attributes) - 1:
                emitter.emit_line('')
    for prop, (getter, setter) in cl.properties.items():
        rtype = getter.sig.ret_type
        emitter.emit_line('')
        generate_readonly_getter(cl, prop, rtype, getter, emitter)
        if setter:
            arg_type = setter.sig.args[1].type
            emitter.emit_line('')
            generate_property_setter(cl, prop, arg_type, setter, emitter)


def generate_getter(cl: ClassIR,
                    attr: str,
                    rtype: RType,
                    emitter: Emitter) -> None:
    attr_field = emitter.attr(attr)
    emitter.emit_line('static PyObject *')
    emitter.emit_line('{}({} *self, void *closure)'.format(getter_name(cl, attr, emitter.names),
                                                           cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    attr_expr = 'self->{}'.format(attr_field)
    emitter.emit_undefined_attr_check(rtype, attr_expr, '==', unlikely=True)
    emitter.emit_line('PyErr_SetString(PyExc_AttributeError,')
    emitter.emit_line('    "attribute {} of {} undefined");'.format(repr(attr),
                                                                    repr(cl.name)))
    emitter.emit_line('return NULL;')
    emitter.emit_line('}')
    emitter.emit_inc_ref('self->{}'.format(attr_field), rtype)
    emitter.emit_box('self->{}'.format(attr_field), 'retval', rtype, declare_dest=True)
    emitter.emit_line('return retval;')
    emitter.emit_line('}')


def generate_setter(cl: ClassIR,
                    attr: str,
                    rtype: RType,
                    emitter: Emitter) -> None:
    attr_field = emitter.attr(attr)
    emitter.emit_line('static int')
    emitter.emit_line('{}({} *self, PyObject *value, void *closure)'.format(
        setter_name(cl, attr, emitter.names),
        cl.struct_name(emitter.names)))
    emitter.emit_line('{')

    deletable = cl.is_deletable(attr)
    if not deletable:
        emitter.emit_line('if (value == NULL) {')
        emitter.emit_line('PyErr_SetString(PyExc_AttributeError,')
        emitter.emit_line('    "{} object attribute {} cannot be deleted");'.format(repr(cl.name),
                                                                                    repr(attr)))
        emitter.emit_line('return -1;')
        emitter.emit_line('}')

    if rtype.is_refcounted:
        attr_expr = 'self->{}'.format(attr_field)
        emitter.emit_undefined_attr_check(rtype, attr_expr, '!=')
        emitter.emit_dec_ref('self->{}'.format(attr_field), rtype)
        emitter.emit_line('}')

    if deletable:
        emitter.emit_line('if (value != NULL) {')
    if rtype.is_unboxed:
        emitter.emit_unbox('value', 'tmp', rtype, error=ReturnHandler('-1'), declare_dest=True)
    elif is_same_type(rtype, object_rprimitive):
        emitter.emit_line('PyObject *tmp = value;')
    else:
        emitter.emit_cast('value', 'tmp', rtype, declare_dest=True)
        emitter.emit_lines('if (!tmp)',
                           '    return -1;')
    emitter.emit_inc_ref('tmp', rtype)
    emitter.emit_line('self->{} = tmp;'.format(attr_field))
    if deletable:
        emitter.emit_line('} else')
        emitter.emit_line('    self->{} = {};'.format(attr_field,
                                                      emitter.c_undefined_value(rtype)))
    emitter.emit_line('return 0;')
    emitter.emit_line('}')


def generate_readonly_getter(cl: ClassIR,
                             attr: str,
                             rtype: RType,
                             func_ir: FuncIR,
                             emitter: Emitter) -> None:
    emitter.emit_line('static PyObject *')
    emitter.emit_line('{}({} *self, void *closure)'.format(getter_name(cl, attr, emitter.names),
                                                           cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    if rtype.is_unboxed:
        emitter.emit_line('{}retval = {}{}((PyObject *) self);'.format(
            emitter.ctype_spaced(rtype), NATIVE_PREFIX, func_ir.cname(emitter.names)))
        emitter.emit_box('retval', 'retbox', rtype, declare_dest=True)
        emitter.emit_line('return retbox;')
    else:
        emitter.emit_line('return {}{}((PyObject *) self);'.format(NATIVE_PREFIX,
                                                                   func_ir.cname(emitter.names)))
    emitter.emit_line('}')


def generate_property_setter(cl: ClassIR,
                             attr: str,
                             arg_type: RType,
                             func_ir: FuncIR,
                             emitter: Emitter) -> None:

    emitter.emit_line('static int')
    emitter.emit_line('{}({} *self, PyObject *value, void *closure)'.format(
        setter_name(cl, attr, emitter.names),
        cl.struct_name(emitter.names)))
    emitter.emit_line('{')
    if arg_type.is_unboxed:
        emitter.emit_unbox('value', 'tmp', arg_type, error=ReturnHandler('-1'),
                           declare_dest=True)
        emitter.emit_line('{}{}((PyObject *) self, tmp);'.format(
                          NATIVE_PREFIX,
                          func_ir.cname(emitter.names)))
    else:
        emitter.emit_line('{}{}((PyObject *) self, value);'.format(
                          NATIVE_PREFIX,
                          func_ir.cname(emitter.names)))
    emitter.emit_line('return 0;')
    emitter.emit_line('}')
