##################### UFuncDefinition ######################

cdef extern from *:
    ctypedef int npy_intp
    struct PyObject
    PyObject* __Pyx_NewRef(object)
    {{inline_func_declaration}}

# variable names have to come from tempita to avoid duplication
@cname("{{func_cname}}")
cdef void {{func_cname}}(char **args, const npy_intp *dimensions, const npy_intp* steps, void* data) except * {{"nogil" if will_be_called_without_gil else ""}}:
    cdef npy_intp i
    cdef npy_intp n = dimensions[0]
    {{for idx, tp in enumerate(in_types)}}
    cdef char* in_{{idx}} = args[{{idx}}]
    cdef {{tp.empty_declaration_code(pyrex=True)}} cast_in_{{idx}}
    {{endfor}}
    {{for idx, tp in enumerate(out_types)}}
    cdef char* out_{{idx}} = args[{{idx+len(in_types)}}]
    cdef {{tp.empty_declaration_code(pyrex=True)}} cast_out_{{idx}}
    {{endfor}}
    {{for idx in range(len(out_types)+len(in_types))}}
    cdef npy_intp step_{{idx}} = steps[{{idx}}]
    {{endfor}}

    {{"with gil" if (not nogil and will_be_called_without_gil) else "if True"}}:
        for i in range(n):
            {{for idx, tp in enumerate(in_types)}}
            {{if tp.is_pyobject}}
            cast_in_{{idx}} = (<{{tp.empty_declaration_code(pyrex=True)}}>(<void**>in_{{idx}})[0])
            {{else}}
            cast_in_{{idx}} = (<{{tp.empty_declaration_code(pyrex=True)}}*>in_{{idx}})[0]
            {{endif}}
            {{endfor}}

            {{", ".join("cast_out_{}".format(idx) for idx in range(len(out_types)))}} = \
                {{inline_func_call}}({{", ".join("cast_in_{}".format(idx) for idx in range(len(in_types)))}})

            {{for idx, tp in enumerate(out_types)}}
            {{if tp.is_pyobject}}
            (<void**>out_{{idx}})[0] = <void*>__Pyx_NewRef(cast_out_{{idx}})
            {{else}}
            (<{{tp.empty_declaration_code(pyrex=True)}}*>out_{{idx}})[0] = cast_out_{{idx}}
            {{endif}}
            {{endfor}}
            {{for idx in range(len(in_types))}}
            in_{{idx}} += step_{{idx}}
            {{endfor}}
            {{for idx in range(len(out_types))}}
            out_{{idx}} += step_{{idx+len(in_types)}}
            {{endfor}}
