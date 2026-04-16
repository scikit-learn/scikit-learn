cimport cython
from ..StringIOTree cimport StringIOTree


cdef class AbstractUtilityCode:
    pass


cdef class UtilityCodeBase(AbstractUtilityCode):
    cpdef format_code(self, code_string, replace_empty_lines=*)


cdef class UtilityCode(UtilityCodeBase):
    cdef public object name
    cdef public object proto
    cdef public object export
    cdef public object impl
    cdef public object init
    cdef public object cleanup
    cdef object proto_block
    cdef readonly object module_state_decls
    cdef readonly object module_state_traverse
    cdef readonly object module_state_clear
    cdef public object requires
    cdef dict _cache
    cdef list specialize_list
    cdef public object file
    cdef readonly tuple _parts_tuple
    cdef list shared_utility_functions

    cpdef none_or_sub(self, s, context)
    # TODO - Signature not compatible with previous declaration
    #@cython.final
    #cdef bint _put_code_section(self, writer, code_type: str) except -1


@cython.final
cdef class FunctionState:
    cdef set names_taken
    cdef public object owner
    cdef public object scope

    cdef public object error_label
    cdef public size_t label_counter
    cdef public set labels_used
    cdef public object return_label
    cdef public object continue_label
    cdef public object break_label
    cdef readonly list yield_labels

    cdef public object return_from_error_cleanup_label # not used in __init__ ?

    cdef public object exc_vars
    cdef public object current_except
    cdef public bint can_trace
    cdef public bint gil_owned

    cdef list temps_allocated
    cdef dict temps_free
    cdef dict temps_used_type
    cdef set zombie_temps
    cdef size_t temp_counter
    cdef list collect_temps_stack

    cdef readonly object closure_temps
    cdef bint should_declare_error_indicator
    cdef public bint uses_error_indicator
    cdef public bint error_without_exception

    cdef public bint needs_refnanny

    cpdef new_label(self, name=*)
    cpdef tuple get_loop_labels(self)
    cpdef set_loop_labels(self, labels)
    cpdef tuple get_all_labels(self)
    cpdef set_all_labels(self, labels)
    cpdef start_collecting_temps(self)
    cpdef stop_collecting_temps(self)

    cpdef list temps_in_use(self)


@cython.final
cdef class PyObjectConst:
    cdef readonly object cname
    cdef readonly object type


@cython.final
cdef class StringConst:
    cdef readonly object cname
    cdef readonly object text
    cdef readonly object escaped_value
    cdef readonly dict py_strings
    cdef public bint c_used

    cpdef get_py_string_const(self, encoding, identifier=*)


@cython.final
cdef class PyStringConst:
    cdef readonly object cname
    cdef readonly object encoding
    cdef readonly bint is_unicode
    cdef readonly bint intern


#class GlobalState(object):

#def funccontext_property(name):

cdef class CCodeWriter(object):
    cdef readonly StringIOTree buffer
    cdef readonly list pyclass_stack
    cdef readonly object globalstate
    cdef readonly FunctionState funcstate
    cdef object code_config
    cdef tuple last_pos
    cdef tuple last_marked_pos
    cdef Py_ssize_t level
    cdef public Py_ssize_t call_level  # debug-only, see Nodes.py
    cdef bint bol

    cpdef write(self, s)
    cdef _write_lines(self, s)
    cpdef _write_to_buffer(self, s)
    cdef put_multilines(self, code)
    cpdef put(self, code)
    cpdef put_safe(self, code)
    cpdef putln(self, code=*, bint safe=*)
    cdef emit_marker(self)
    cdef _build_marker(self, tuple pos)
    cdef increase_indent(self)
    cdef decrease_indent(self)
    cdef indent(self)


cdef class PyrexCodeWriter:
    cdef readonly object f
    cdef readonly Py_ssize_t level


cdef class PyxCodeWriter:
    cdef public StringIOTree buffer
    cdef public object context
    cdef object encoding
    cdef Py_ssize_t level
    cdef Py_ssize_t original_level
    cdef dict _insertion_points
