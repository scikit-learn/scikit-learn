import copy
import hashlib

from . import (ExprNodes, PyrexTypes,
               ParseTreeTransforms, StringEncoding, Errors,
               Naming)
from .ExprNodes import CloneNode, CodeObjectNode, ProxyNode, TupleNode
from .Nodes import FuncDefNode, StatListNode, DefNode
from ..Utils import OrderedSet
from .Errors import error, CannotSpecialize


class FusedCFuncDefNode(StatListNode):
    """
    This node replaces a function with fused arguments. It deep-copies the
    function for every permutation of fused types, and allocates a new local
    scope for it. It keeps track of the original function in self.node, and
    the entry of the original function in the symbol table is given the
    'fused_cfunction' attribute which points back to us.
    Then when a function lookup occurs (to e.g. call it), the call can be
    dispatched to the right function.

    node    FuncDefNode    the original function
    nodes   [FuncDefNode]  list of copies of node with different specific types
    py_func DefNode        the fused python function subscriptable from
                           Python space
    __signatures__         A DictNode mapping signature specialization strings
                           to PyCFunction nodes
    resulting_fused_function  PyCFunction for the fused DefNode that delegates
                              to specializations
    fused_func_assignment   Assignment of the fused function to the function name
    defaults_tuple          TupleNode of defaults (letting PyCFunctionNode build
                            defaults would result in many different tuples)
    specialized_pycfuncs    List of synthesized pycfunction nodes for the
                            specializations

    fused_compound_types    All fused (compound) types (e.g. floating[:])
    """

    __signatures__ = None
    resulting_fused_function = None
    fused_func_assignment = None
    py_func = None
    defaults_tuple = None
    decorators = None

    child_attrs = StatListNode.child_attrs + [
        '__signatures__', 'resulting_fused_function', 'fused_func_assignment']

    def __init__(self, node, env):
        super().__init__(node.pos)

        self.nodes = []
        self.node = node

        is_def = isinstance(self.node, DefNode)
        if is_def:
            # self.node.decorators = []
            self.copy_def(env)
        else:
            self.copy_cdef(env)

        # Perform some sanity checks. If anything fails, it's a bug
        for n in self.nodes:
            assert not n.entry.type.is_fused
            assert not n.local_scope.return_type.is_fused
            if node.return_type.is_fused:
                assert not n.return_type.is_fused

            if not is_def and n.cfunc_declarator.optional_arg_count:
                assert n.type.op_arg_struct

        node.entry.fused_cfunction = self
        # Copy the nodes as AnalyseDeclarationsTransform will prepend
        # self.py_func to self.stats, as we only want specialized
        # CFuncDefNodes in self.nodes
        self.stats = self.nodes[:]

    def copy_def(self, env):
        """
        Create a copy of the original def or lambda function for specialized
        versions.
        """
        fused_compound_types = PyrexTypes.unique(
            [arg.type for arg in self.node.args if arg.type.is_fused])
        fused_types = self._get_fused_base_types(fused_compound_types)
        permutations = PyrexTypes.get_all_specialized_permutations(fused_types)

        self.fused_compound_types = fused_compound_types

        if self.node.entry in env.pyfunc_entries:
            env.pyfunc_entries.remove(self.node.entry)

        for cname, fused_to_specific in permutations:
            copied_node = copy.deepcopy(self.node)
            # keep signature object identity for special casing in DefNode.analyse_declarations()
            copied_node.entry.signature = self.node.entry.signature

            self._specialize_function_args(copied_node.args, fused_to_specific)
            copied_node.return_type = self.node.return_type.specialize(
                                                    fused_to_specific)
            copied_node.code_object = CodeObjectNode(copied_node)
            copied_node.analyse_declarations(env)
            # copied_node.is_staticmethod = self.node.is_staticmethod
            # copied_node.is_classmethod = self.node.is_classmethod
            self.create_new_local_scope(copied_node, env, fused_to_specific)
            self.specialize_copied_def(copied_node, cname, self.node.entry,
                                       fused_to_specific, fused_compound_types)

            PyrexTypes.specialize_entry(copied_node.entry, cname)
            copied_node.entry.used = True
            env.entries[copied_node.entry.name] = copied_node.entry

            specialised_type_names = [
                sarg.type.declaration_code('', for_display=True)
                for (farg, sarg) in zip(self.node.args, copied_node.args)
                if farg.type.is_fused
            ]
            copied_node.name = StringEncoding.EncodedString(f"{copied_node.name}[{','.join(specialised_type_names)}]")

            if not self.replace_fused_typechecks(copied_node):
                break

        self.orig_py_func = self.node
        self.py_func = self.make_fused_cpdef(self.node, env, is_def=True)

    def copy_cdef(self, env):
        """
        Create a copy of the original c(p)def function for all specialized
        versions.
        """
        permutations = self.node.type.get_all_specialized_permutations()
        # print 'Node %s has %d specializations:' % (self.node.entry.name,
        #                                            len(permutations))
        # import pprint; pprint.pprint([d for cname, d in permutations])

        # Prevent copying of the python function
        self.orig_py_func = orig_py_func = self.node.py_func
        self.node.py_func = None
        if orig_py_func:
            env.pyfunc_entries.remove(orig_py_func.entry)

        fused_types = self.node.type.get_fused_types()
        self.fused_compound_types = fused_types

        new_cfunc_entries = []
        for cname, fused_to_specific in permutations:
            copied_node = copy.deepcopy(self.node)

            # Make the types in our CFuncType specific.
            try:
                type = copied_node.type.specialize(fused_to_specific)
            except CannotSpecialize:
                # unlike for the argument types, specializing the return type can fail
                error(copied_node.pos, "Return type is a fused type that cannot "
                      "be determined from the function arguments")
                self.py_func = None  # this is just to let the compiler exit gracefully
                return
            entry = copied_node.entry
            type.specialize_entry(entry, cname)

            # Reuse existing Entries (e.g. from .pxd files).
            for orig_entry in env.cfunc_entries:
                if entry.cname == orig_entry.cname and type.same_as_resolved_type(orig_entry.type):
                    copied_node.entry = orig_entry
                    if not copied_node.entry.func_cname:
                        copied_node.entry.func_cname = entry.func_cname
                    entry = orig_entry
                    type = orig_entry.type
                    break
            else:
                new_cfunc_entries.append(entry)

            copied_node.type = type
            entry.type, type.entry = type, entry

            entry.used = (entry.used or
                          self.node.entry.defined_in_pxd or
                          env.is_c_class_scope or
                          entry.is_cmethod)

            if self.node.cfunc_declarator.optional_arg_count:
                self.node.cfunc_declarator.declare_optional_arg_struct(
                                           type, env, fused_cname=cname)

            copied_node.return_type = type.return_type
            self.create_new_local_scope(copied_node, env, fused_to_specific)

            # Make the argument types in the CFuncDeclarator specific
            self._specialize_function_args(copied_node.cfunc_declarator.args,
                                           fused_to_specific)

            # If a cpdef, declare all specialized cpdefs (this
            # also calls analyse_declarations)
            copied_node.declare_cpdef_wrapper(env)
            if copied_node.py_func:
                env.pyfunc_entries.remove(copied_node.py_func.entry)

                self.specialize_copied_def(
                        copied_node.py_func, cname, self.node.entry.as_variable,
                        fused_to_specific, fused_types)

            if not self.replace_fused_typechecks(copied_node):
                break

        # replace old entry with new entries
        if self.node.entry in env.cfunc_entries:
            cindex = env.cfunc_entries.index(self.node.entry)
            env.cfunc_entries[cindex:cindex+1] = new_cfunc_entries
        else:
            env.cfunc_entries.extend(new_cfunc_entries)

        if orig_py_func:
            self.py_func = self.make_fused_cpdef(orig_py_func, env,
                                                 is_def=False)
        else:
            self.py_func = orig_py_func

    def _get_fused_base_types(self, fused_compound_types):
        """
        Get a list of unique basic fused types, from a list of
        (possibly) compound fused types.
        """
        base_types = []
        seen = set()
        for fused_type in fused_compound_types:
            fused_type.get_fused_types(result=base_types, seen=seen)
        return base_types

    def _specialize_function_args(self, args, fused_to_specific):
        for arg in args:
            if arg.type.is_fused:
                arg.type = arg.type.specialize(fused_to_specific)
                if arg.type.is_memoryviewslice:
                    arg.type.validate_memslice_dtype(arg.pos)
                if arg.annotation:
                    # TODO might be nice if annotations were specialized instead?
                    # (Or might be hard to do reliably)
                    arg.annotation.untyped = True

    def create_new_local_scope(self, node, env, f2s):
        """
        Create a new local scope for the copied node and append it to
        self.nodes. A new local scope is needed because the arguments with the
        fused types are already in the local scope, and we need the specialized
        entries created after analyse_declarations on each specialized version
        of the (CFunc)DefNode.
        f2s is a dict mapping each fused type to its specialized version
        """
        node.create_local_scope(env)
        node.local_scope.fused_to_specific = f2s

        # This is copied from the original function, set it to false to
        # stop recursion
        node.has_fused_arguments = False
        self.nodes.append(node)

    def specialize_copied_def(self, node, cname, py_entry, f2s, fused_compound_types):
        """Specialize the copy of a DefNode given the copied node,
        the specialization cname and the original DefNode entry"""
        fused_types = self._get_fused_base_types(fused_compound_types)
        type_strings = [
            PyrexTypes.specialization_signature_string(fused_type, f2s)
                for fused_type in fused_types
        ]

        node.specialized_signature_string = '|'.join(type_strings)

        node.entry.pymethdef_cname = PyrexTypes.get_fused_cname(
                                        cname, node.entry.pymethdef_cname)
        node.entry.doc = py_entry.doc
        node.entry.doc_cname = py_entry.doc_cname

    def replace_fused_typechecks(self, copied_node):
        """
        Branch-prune fused type checks like

            if fused_t is int:
                ...

        Returns whether an error was issued and whether we should stop in
        in order to prevent a flood of errors.
        """
        num_errors = Errors.get_errors_count()
        transform = ParseTreeTransforms.ReplaceFusedTypeChecks(
                                       copied_node.local_scope)
        transform(copied_node)

        if Errors.get_errors_count() > num_errors:
            return False

        return True

    def _fused_instance_checks(self, normal_types, pyx_code, env):
        """
        Generate Cython code for instance checks, matching an object to
        specialized types.
        """
        for specialized_type in normal_types:
            # all_numeric = all_numeric and specialized_type.is_numeric
            py_type_name = specialized_type.py_type_name()
            pyx_code.put_chunk(
                f"""
                    if isinstance(arg, {py_type_name}):
                        return '{specialized_type.specialization_string}'
                """
            )

    def _dtype_name(self, dtype):
        name = str(dtype).replace('_', '__').replace(' ', '_')
        if dtype.is_typedef:
            name = Naming.fused_dtype_prefix + name
        return name

    def _dtype_type(self, dtype):
        if dtype.is_typedef:
            return self._dtype_name(dtype)
        return str(dtype)

    def _sizeof_dtype(self, dtype):
        if dtype.is_pyobject:
            return 'sizeof(void *)'
        else:
            return f"sizeof({self._dtype_type(dtype)})"

    def _buffer_check_numpy_dtype_setup_cases(self, pyx_code):
        "Setup some common cases to match dtypes against specializations"
        with pyx_code.indenter("if kind in u'iu':"):
            pyx_code.putln("pass")
            pyx_code.named_insertion_point("dtype_int")

        with pyx_code.indenter("elif kind == u'f':"):
            pyx_code.putln("pass")
            pyx_code.named_insertion_point("dtype_float")

        with pyx_code.indenter("elif kind == u'c':"):
            pyx_code.putln("pass")
            pyx_code.named_insertion_point("dtype_complex")

    def _buffer_check_numpy_dtype(self, pyx_code, specialized_buffer_types, pythran_types):
        """
        Match a numpy dtype object to the individual specializations.
        """
        self._buffer_check_numpy_dtype_setup_cases(pyx_code)

        for specialized_type in pythran_types+specialized_buffer_types:
            final_type = specialized_type
            if specialized_type.is_pythran_expr:
                specialized_type = specialized_type.org_buffer
            dtype = specialized_type.dtype

            itemsize_match = self._sizeof_dtype(dtype) + " == itemsize"
            signed_match = f" and not ({self._dtype_name(dtype)}_is_signed ^ dtype_signed)"

            dtypes = [
                (dtype.is_int, pyx_code['dtype_int']),
                (dtype.is_float, pyx_code['dtype_float']),
                (dtype.is_complex, pyx_code['dtype_complex'])
            ]

            for dtype_category, codewriter in dtypes:
                if not dtype_category:
                    continue

                cond = f'{itemsize_match} and (<Py_ssize_t>arg.ndim) == {specialized_type.ndim}'
                if dtype.is_int:
                    cond += signed_match
                if final_type.is_pythran_expr:
                    cond += ' and arg_is_pythran_compatible'

                with codewriter.indenter(f"if {cond}:"):
                    #codewriter.putln("print 'buffer match found based on numpy dtype'")
                    codewriter.putln(f"return '{final_type.specialization_string}'")

    def _buffer_parse_format_string_check(self, pyx_code, decl_code, specialized_type, env):
        """
        For each specialized type, try to coerce the object to a memoryview
        slice of that type. This means obtaining a buffer and parsing the
        format string.
        TODO: separate buffer acquisition from format parsing
        """
        dtype = specialized_type.dtype
        if specialized_type.is_buffer:
            axes = [('direct', 'strided')] * specialized_type.ndim
        else:
            axes = specialized_type.axes

        memslice_type = PyrexTypes.MemoryViewSliceType(dtype, axes)
        memslice_type.create_from_py_utility_code(env)
        coerce_from_py_func = memslice_type.from_py_function

        decl_code.putln(
            f"{Naming.memviewslice_cname} {coerce_from_py_func}(object, int)")

        match = specialized_type.specialization_string
        sizeof_dtype = self._sizeof_dtype(dtype)
        ndim_dtype = specialized_type.ndim

        # Use the memoryview object to check itemsize and ndim.
        # In principle it could check more, but these are the easiest to do quickly.
        pyx_code.put_chunk(
            f"""
                # try {dtype}
                if (((itemsize == -1 and arg_as_memoryview.itemsize == {sizeof_dtype})
                        or itemsize == {sizeof_dtype})
                        and arg_as_memoryview.ndim == {ndim_dtype}):
                    memslice = {coerce_from_py_func}(arg_as_memoryview, 0)
                    if memslice.memview:
                        __PYX_XCLEAR_MEMVIEW(&memslice, 1)
                        # print 'found a match for the buffer through format parsing'
                        return '{match}'
                    else:
                        __pyx_PyErr_Clear()
            """
        )

    def _buffer_checks(self, buffer_types, pythran_types, pyx_code, decl_code, accept_none, env):
        """
        Generate Cython code to match objects to buffer specializations.
        First try to get a numpy dtype object and match it against the individual
        specializations. If that fails, try naively to coerce the object
        to each specialization, which obtains the buffer each time and tries
        to match the format string.
        """
        # The first thing to find a match in this loop breaks out of the loop
        pyx_code.put_chunk(
            """
                """ + ("arg_is_pythran_compatible = False" if pythran_types else "") + """
                if ndarray is not None:
                    if isinstance(arg, ndarray):
                        dtype = arg.dtype
                        """ + ("arg_is_pythran_compatible = True" if pythran_types else "") + """
                    elif __pyx_memoryview_check(arg):
                        arg_base = arg.base
                        if isinstance(arg_base, ndarray):
                            dtype = arg_base.dtype
                        else:
                            dtype = None
                    else:
                        dtype = None

                    itemsize = -1
                    if dtype is not None:
                        itemsize = dtype.itemsize
                        kind = ord(dtype.kind)
                        dtype_signed = kind == u'i'
            """)
        pyx_code.indent(2)
        if pythran_types:
            pyx_code.put_chunk(
                """
                        # Pythran only supports the endianness of the current compiler
                        byteorder = dtype.byteorder
                        if byteorder == "<" and not __Pyx_Is_Little_Endian():
                            arg_is_pythran_compatible = False
                        elif byteorder == ">" and __Pyx_Is_Little_Endian():
                            arg_is_pythran_compatible = False
                        if arg_is_pythran_compatible:
                            cur_stride = itemsize
                            shape = arg.shape
                            strides = arg.strides
                            for i in range(arg.ndim-1, -1, -1):
                                if (<Py_ssize_t>strides[i]) != cur_stride:
                                    arg_is_pythran_compatible = False
                                    break
                                cur_stride *= <Py_ssize_t> shape[i]
                            else:
                                arg_is_pythran_compatible = not (arg.flags.f_contiguous and (<Py_ssize_t>arg.ndim) > 1)
                """)
        self._buffer_check_numpy_dtype(pyx_code, buffer_types, pythran_types)
        pyx_code.dedent(2)

        if accept_none:
            # If None is acceptable, then Cython <3.0 matched None with the
            # first type. This behaviour isn't ideal, but keep it for backwards
            # compatibility. Better behaviour would be to see if subsequent
            # arguments give a stronger match.
            pyx_code.put_chunk(
                f"""
                if arg is None:
                    return '{buffer_types[0].specialization_string}'
                """
            )

        # creating a Cython memoryview from a Python memoryview avoids the
        # need to get the buffer multiple times, and we can
        # also use it to check itemsizes etc
        pyx_code.put_chunk(
            """
            try:
                arg_as_memoryview = memoryview(arg)
            except (ValueError, TypeError):
                pass
            """)
        with pyx_code.indenter("else:"):
            for specialized_type in buffer_types:
                self._buffer_parse_format_string_check(
                        pyx_code, decl_code, specialized_type, env)

    def _buffer_declarations(self, pyx_code, decl_code, all_buffer_types, pythran_types):
        """
        If we have any buffer specializations, write out some variable and type declarations.
        """
        decl_code.put_chunk(
            f"""
                ctypedef struct {Naming.memviewslice_cname}:
                    void *memview

                void __PYX_XCLEAR_MEMVIEW({Naming.memviewslice_cname} *, int have_gil)
                bint __pyx_memoryview_check(object)
            """)

        pyx_code['local_variable_declarations'].put_chunk(
            f"""
                cdef {Naming.memviewslice_cname} memslice
                cdef Py_ssize_t itemsize
                cdef bint dtype_signed
                cdef Py_UCS4 kind

                itemsize = -1
            """)

        if pythran_types:
            pyx_code['local_variable_declarations'].put_chunk("""
                cdef bint arg_is_pythran_compatible
                cdef Py_ssize_t cur_stride
                cdef Py_ssize_t i
            """)

        pyx_code['local_variable_declarations'].put_chunk(
            """
                cdef memoryview arg_as_memoryview
            """
        )

        seen_typedefs = set()
        seen_int_dtypes = set()
        seen_structs = set()
        for buffer_type in all_buffer_types:
            dtype = buffer_type.dtype
            dtype_name = self._dtype_name(dtype)
            if dtype.is_struct_or_union:
                if dtype_name not in seen_structs:
                    seen_structs.add(dtype_name)
                    decl_code.putln(
                        f'ctypedef {dtype.kind} {dtype_name} "{dtype.empty_declaration_code()}": pass')

            elif dtype.is_typedef:
                if dtype_name not in seen_typedefs:
                    seen_typedefs.add(dtype_name)
                    decl_code.putln(
                        f'ctypedef {dtype.resolve()} {dtype_name} "{dtype.empty_declaration_code()}"')

            # 'is_signed' is also needed for typedefs.
            if dtype.is_int:
                if str(dtype) not in seen_int_dtypes:
                    seen_int_dtypes.add(str(dtype))
                    dtype_type = self._dtype_type(dtype)
                    pyx_code['local_variable_declarations'].put_chunk(
                        f"""
                            cdef bint {dtype_name}_is_signed
                            {dtype_name}_is_signed = not (<{dtype_type}> -1 > 0)
                        """)

    def _split_fused_types(self, arg):
        """
        Specialize fused types and split into normal types and buffer types.
        """
        specialized_types = PyrexTypes.get_specialized_types(arg.type)

        # Prefer long over int, etc by sorting (see type classes in PyrexTypes.py)
        specialized_types.sort()

        seen_py_type_names = set()
        normal_types, buffer_types, pythran_types = [], [], []
        has_object_fallback = False
        for specialized_type in specialized_types:
            py_type_name = specialized_type.py_type_name()
            if py_type_name:
                if py_type_name in seen_py_type_names:
                    continue
                seen_py_type_names.add(py_type_name)
                if py_type_name == 'object':
                    has_object_fallback = True
                else:
                    normal_types.append(specialized_type)
            elif specialized_type.is_pythran_expr:
                pythran_types.append(specialized_type)
            elif specialized_type.is_buffer or specialized_type.is_memoryviewslice:
                buffer_types.append(specialized_type)

        return normal_types, buffer_types, pythran_types, has_object_fallback

    def _unpack_argument(self, pyx_code, arg, arg_tuple_idx, min_positional_args, default_idx):
        pyx_code.put_chunk(
            f"""
                # PROCESSING ARGUMENT {arg_tuple_idx}
                if {arg_tuple_idx} < len(<tuple>args):
                    arg = (<tuple>args)[{arg_tuple_idx}]
                elif kwargs is not None and '{arg.name}' in <dict>kwargs:
                    arg = (<dict>kwargs)['{arg.name}']
                else:
            """
        )
        pyx_code.indent()
        if arg.default:
            pyx_code.putln(
                f"arg = (<tuple>defaults)[{default_idx}]")
        elif arg_tuple_idx < min_positional_args:
            pyx_code.putln(
                'raise TypeError("Expected at least %d argument%s, got %d" % ('
                f'''{min_positional_args}, {'"s"' if min_positional_args != 1 else '""'}, len(<tuple>args)))'''
            )
        else:
            pyx_code.putln(f"""raise TypeError("Missing keyword-only argument: '%s'" % "{arg.name}")""")
        pyx_code.dedent()

    def make_fused_cpdef(self, orig_py_func, env, is_def):
        """
        This creates the function that is indexable from Python and does
        runtime dispatch based on the argument types. The function gets the
        arg tuple and kwargs dict (or None) and the defaults tuple
        as arguments from the Binding Fused Function's tp_call.
        """
        from . import TreeFragment, Code, UtilityCode

        min_positional_args = (
            self.node.num_required_args - self.node.num_required_kw_args
            if is_def else
            sum(1 for arg in self.node.args if arg.default is None)
        )

        pyx_code = Code.PyxCodeWriter()
        decl_code = Code.PyxCodeWriter()
        type_mapper = Code.PyxCodeWriter()
        decl_code.put_chunk(
            """
                cdef extern from *:
                    type __Pyx_ImportNumPyArrayTypeIfAvailable()

                    # from FusedFunction utility code
                    object __pyx_ff_match_signatures_single(dict signatures, dest_type)
                    object __pyx_ff_match_signatures(dict signatures, tuple dest_sig, dict sigindex)
            """)
        decl_code.indent()

        pyx_code.put_chunk(
            """
                def __pyx_fused_cpdef(signatures, args, kwargs, defaults, _fused_sigindex={}):
                    # FIXME: use a typed signature - currently fails badly because
                    #        default arguments inherit the types we specify here!

                    if kwargs is not None and not kwargs:
                        kwargs = None

                    # instance check body
            """)

        pyx_code.indent()  # indent following code to function body
        pyx_code.named_insertion_point("imports")

        fused_index = 0
        default_idx = 0
        all_buffer_types = OrderedSet()
        seen_fused_types = set()
        for i, arg in enumerate(self.node.args):
            if arg.type.is_fused:
                arg_fused_types = arg.type.get_fused_types()
                if len(arg_fused_types) > 1:
                    raise NotImplementedError("Determination of more than one fused base "
                                              "type per argument is not implemented.")
                fused_type = arg_fused_types[0]

            if arg.type.is_fused and fused_type not in seen_fused_types:
                seen_fused_types.add(fused_type)

                normal_types, buffer_types, pythran_types, has_object_fallback = self._split_fused_types(arg)
                self._unpack_argument(pyx_code, arg, i, min_positional_args, default_idx)

                mapper_arg_types = ['object', 'type']
                mapper_arg_names = ['arg']
                if buffer_types or pythran_types:
                    mapper_arg_names.append('ndarray')

                mapper_sig = ', '.join(f"{atype} {aname}" for atype, aname in zip(mapper_arg_types, mapper_arg_names))
                mapper_args = ', '.join(mapper_arg_names)

                mapper_decl_code = type_mapper.insertion_point()
                mapper_decl_code.put_chunk(
                    """
                    cdef extern from *:
                        void __pyx_PyErr_Clear "PyErr_Clear" ()
                        int __Pyx_Is_Little_Endian()
                    """
                )
                mapper_decl_code.indent()

                type_mapper.putln('')
                type_mapper.putln("@TYPE_MAPPER_CNAME_PLACEHOLDER")
                with type_mapper.indenter(f"cdef str map_fused_type({mapper_sig}):"):

                    type_mapper.named_insertion_point("local_variable_declarations")

                    if normal_types:
                        self._fused_instance_checks(normal_types, type_mapper, env)

                    if buffer_types or pythran_types:
                        mapper_buffer_types = OrderedSet()
                        mapper_buffer_types.update(buffer_types)
                        mapper_buffer_types.update(ty.org_buffer for ty in pythran_types)

                        self._buffer_declarations(type_mapper, mapper_decl_code, mapper_buffer_types, pythran_types)

                        self._buffer_checks(
                            buffer_types, pythran_types, type_mapper, mapper_decl_code,
                            arg.accept_none, env)

                    type_mapper.putln("return 'object'" if has_object_fallback else "return None")

                type_mapper_impl = type_mapper.getvalue()
                type_mapper.reset()

                # Generate a unique name for the mapper function based on type declarations and implementation.
                impl_hash = hashlib.sha256(type_mapper_impl.encode('utf-8')).hexdigest()
                type_mapper_cname = (
                    f"__pyx_ff_map_fused_{impl_hash[:6]}"
                    f"_{len(mapper_arg_names)}_{len(fused_type.types)}"
                    f"_{PyrexTypes.type_list_identifier(fused_type.types)}"
                )
                type_mapper_impl = type_mapper_impl.replace(
                    "\n@TYPE_MAPPER_CNAME_PLACEHOLDER\n", f"\n@cname('{type_mapper_cname}')\n")
                # print(''.join(f"{i:3d}  {line}" for i, line in enumerate(type_mapper_impl.splitlines(keepends=True))))

                env.use_utility_code(
                    UtilityCode.CythonUtilityCode(type_mapper_impl, name=type_mapper_cname))

                decl_code.putln(f"str {type_mapper_cname}({mapper_sig})")
                pyx_code.putln(f"dest_sig{fused_index} = {type_mapper_cname}({mapper_args})")

                fused_index += 1
                all_buffer_types.update(buffer_types)
                all_buffer_types.update(ty.org_buffer for ty in pythran_types)

            if arg.default:
                default_idx += 1

        if all_buffer_types:
            env.use_utility_code(
                Code.UtilityCode.load_cached("IsLittleEndian", "ModuleSetupCode.c"))
            env.use_utility_code(
                Code.UtilityCode.load_cached("Import", "ImportExport.c"))
            env.use_utility_code(
                Code.UtilityCode.load_cached("ImportNumPyArray", "ImportExport.c"))

            pyx_code['imports'].put_chunk(
                """
                    cdef type ndarray
                    ndarray = __Pyx_ImportNumPyArrayTypeIfAvailable()
                """)

        if len(seen_fused_types) == 1:
            # Fast and common case: a single fused type across all arguments.
            env.use_utility_code(
                UtilityCode.CythonUtilityCode.load("match_signatures_single", "FusedFunction.pyx"))
            pyx_code.put_chunk(
                """
                return __pyx_ff_match_signatures_single(<dict> signatures, dest_sig0)
                """
            )
        else:
            env.use_utility_code(
                UtilityCode.CythonUtilityCode.load("match_signatures", "FusedFunction.pyx"))
            dest_sig_tuple = ', '.join(f'dest_sig{i}' for i in range(len(seen_fused_types)))
            pyx_code.put_chunk(
                f"""
                return __pyx_ff_match_signatures(<dict> signatures, ({dest_sig_tuple}), <dict> _fused_sigindex)
                """
            )

        fragment_code = pyx_code.getvalue()
        # print(decl_code.getvalue())
        # print(fragment_code)
        # print(''.join(f"{i:3d}  {line}" for i, line in enumerate(fragment_code.splitlines(keepends=True))))
        from .Optimize import ConstantFolding
        fragment = TreeFragment.TreeFragment(
            fragment_code, level='module', pipeline=[ConstantFolding()])
        ast = TreeFragment.SetPosTransform(self.node.pos)(fragment.root)
        UtilityCode.declare_declarations_in_scope(
            decl_code.getvalue(), env.global_scope())
        ast.scope = env
        # FIXME: for static methods of cdef classes, we build the wrong signature here: first arg becomes 'self'
        ast.analyse_declarations(env)
        py_func = ast.stats[-1]  # the DefNode
        self.fragment_scope = ast.scope

        if isinstance(self.node, DefNode):
            py_func.specialized_cpdefs = self.nodes[:]
        else:
            py_func.specialized_cpdefs = [n.py_func for n in self.nodes]

        return py_func

    def update_fused_defnode_entry(self, env):
        copy_attributes = (
            'name', 'pos', 'cname', 'func_cname', 'pyfunc_cname',
            'pymethdef_cname', 'doc', 'doc_cname', 'is_member',
            'scope'
        )

        entry = self.py_func.entry

        for attr in copy_attributes:
            setattr(entry, attr,
                    getattr(self.orig_py_func.entry, attr))

        self.py_func.name = self.orig_py_func.name
        self.py_func.doc = self.orig_py_func.doc

        env.entries.pop('__pyx_fused_cpdef', None)
        if isinstance(self.node, DefNode):
            env.entries[entry.name] = entry
        else:
            env.entries[entry.name].as_variable = entry

        env.pyfunc_entries.append(entry)

        self.py_func.entry.fused_cfunction = self
        def_nodes = []
        for node in self.nodes:
            if isinstance(self.node, DefNode):
                def_nodes.append(node)
                node.fused_py_func = self.py_func
            else:
                def_nodes.append(node.py_func)
                node.py_func.fused_py_func = self.py_func
                node.entry.as_variable = entry

        self.synthesize_defnodes(def_nodes)

    def analyse_expressions(self, env):
        """
        Analyse the expressions. Take care to only evaluate default arguments
        once and clone the result for all specializations
        """
        for fused_compound_type in self.fused_compound_types:
            for fused_type in fused_compound_type.get_fused_types():
                for specialization_type in fused_type.types:
                    if specialization_type.is_complex:
                        specialization_type.create_declaration_utility_code(env)

        if self.py_func:
            self.__signatures__ = self.__signatures__.analyse_expressions(env)
            self.py_func = self.py_func.analyse_expressions(env)
            self.resulting_fused_function = self.resulting_fused_function.analyse_expressions(env)
            self.fused_func_assignment = self.fused_func_assignment.analyse_expressions(env)

        self.defaults = defaults = []

        for arg in self.node.args:
            if arg.default:
                arg.default = arg.default.analyse_expressions(env)
                if arg.default.is_literal:
                    defaults.append(copy.copy(arg.default))
                else:
                    # coerce the argument to temp since CloneNode really requires a temp
                    defaults.append(ProxyNode(arg.default.coerce_to_temp(env)))
            else:
                defaults.append(None)

        for i, stat in enumerate(self.stats):
            stat = self.stats[i] = stat.analyse_expressions(env)
            if isinstance(stat, FuncDefNode) and stat is not self.py_func:
                # the dispatcher specifically doesn't want its defaults overriding
                for arg, default in zip(stat.args, defaults):
                    if default is not None:
                        if default.is_literal:
                            arg.default = default.coerce_to(arg.type, env)
                        else:
                            arg.default = CloneNode(default).analyse_expressions(env).coerce_to(arg.type, env)

        if self.py_func:
            args = [CloneNode(default) for default in defaults if default]
            self.defaults_tuple = TupleNode(self.pos, args=args)
            self.defaults_tuple = self.defaults_tuple.analyse_types(env, skip_children=True).coerce_to_pyobject(env)
            self.defaults_tuple = ProxyNode(self.defaults_tuple)

            fused_func = self.resulting_fused_function.arg
            fused_func.defaults_tuple = CloneNode(self.defaults_tuple)

            for i, pycfunc in enumerate(self.specialized_pycfuncs):
                pycfunc = self.specialized_pycfuncs[i] = pycfunc.analyse_types(env)
                pycfunc.defaults_tuple = CloneNode(self.defaults_tuple)
        return self

    def synthesize_defnodes(self, nodes):
        """
        Create the __signatures__ dict of PyCFunctionNode specializations.
        """
        # For the moment, fused functions do not support METH_FASTCALL
        for node in nodes:
            node.entry.signature.use_fastcall = False

        signatures = [StringEncoding.EncodedString(node.specialized_signature_string)
                      for node in nodes]
        keys = [ExprNodes.UnicodeNode(node.pos, value=sig)
                for node, sig in zip(nodes, signatures)]
        values = [ExprNodes.PyCFunctionNode.from_defnode(node, binding=True)
                  for node in nodes]

        self.__signatures__ = ExprNodes.DictNode.from_pairs(self.pos, zip(keys, values))

        self.specialized_pycfuncs = values
        for pycfuncnode in values:
            pycfuncnode.is_specialization = True
        # use code object from first defnode to get as close to a correct signature as possible
        self.py_func.code_object = CodeObjectNode(nodes[0])

    def generate_function_definitions(self, env, code):
        if self.py_func:
            self.py_func.pymethdef_required = True
            self.fused_func_assignment.generate_function_definitions(env, code)

        from . import Options
        for stat in self.stats:
            if isinstance(stat, FuncDefNode) and (
                    stat.entry.used or
                    (Options.cimport_from_pyx and not stat.entry.visibility == 'extern')):
                code.mark_pos(stat.pos)
                stat.generate_function_definitions(env, code)

    def generate_execution_code(self, code):
        # Note: all def function specialization are wrapped in PyCFunction
        # nodes in the self.__signatures__ dictnode.
        for default in self.defaults:
            if default is not None:
                default.generate_evaluation_code(code)

        if self.py_func:
            self.defaults_tuple.generate_evaluation_code(code)

        super().generate_execution_code(code)

        if self.__signatures__:
            signatures = self.__signatures__
            signatures.generate_evaluation_code(code)
            fused_func = self.resulting_fused_function
            fused_func.generate_evaluation_code(code)

            code.putln(
                f"((__pyx_FusedFunctionObject *) {fused_func.result()})->__signatures__ = {signatures.result()};")

            signatures.generate_giveref(code)
            signatures.generate_post_assignment_code(code)
            signatures.free_temps(code)

            self.fused_func_assignment.generate_execution_code(code)

            # Dispose of results
            fused_func.generate_disposal_code(code)
            fused_func.free_temps(code)
            self.defaults_tuple.generate_disposal_code(code)
            self.defaults_tuple.free_temps(code)

        for default in self.defaults:
            if default is not None:
                default.generate_disposal_code(code)
                default.free_temps(code)

    def annotate(self, code):
        for stat in self.stats:
            stat.annotate(code)
