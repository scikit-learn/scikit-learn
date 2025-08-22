import cython

from ..Plex.Scanners cimport Scanner

cdef unicode any_string_prefix, IDENT

cdef get_lexicon()
cdef initial_compile_time_env()

## methods commented with '##' out are used by Parsing.py when compiled.

@cython.final
cdef class CompileTimeScope:
    cdef public dict entries
    cdef public CompileTimeScope outer
    ##cdef declare(self, name, value)
    ##cdef lookup_here(self, name)
    ##cpdef lookup(self, name)

@cython.final
cdef class PyrexScanner(Scanner):
    cdef public context
    cdef public list included_files
    cdef public CompileTimeScope compile_time_env
    cdef public bint compile_time_eval
    cdef public bint compile_time_expr
    cdef public bint parse_comments
    cdef public bint in_python_file
    cdef public source_encoding
    cdef dict keywords
    cdef public list indentation_stack
    cdef public Py_UCS4 indentation_char
    cdef public int bracket_nesting_level
    cdef readonly bint async_enabled
    cdef public unicode sy
    cdef public systring  # EncodedString
    cdef public list put_back_on_failure

    cdef Py_ssize_t current_level(self)
    cdef int error_at_scanpos(self, str message) except -1
