#################### match_signatures_single ####################

@cname("__pyx_ff_match_signatures_single")
cdef object match_signatures_single(signatures: dict, dest_type):
    found_match = signatures.get(dest_type)
    if found_match is None:
        raise TypeError("No matching signature found")
    return found_match


#################### match_signatures ####################

cimport cython

@cname("__pyx_ff_index_signature")
cdef object index_signature(dest_sig: tuple, signatures: dict, sigindex: dict):
    matched_function = None
    for sig, function in signatures.items():
        types = (<str> sig).strip('()').split('|')

        for type_name, dest_type in zip(types, dest_sig):
            if dest_type is None:
                continue
            if dest_type != type_name:
                break
        else:
            if matched_function is not None:
                raise TypeError("Function call with ambiguous argument types")
            matched_function = function

    if matched_function is None:
        raise TypeError("No matching signature found")

    sigindex[dest_sig] = matched_function
    return matched_function


@cname("__pyx_ff_match_signatures")
cdef object match_signatures(signatures: dict, dest_sig: tuple, sigindex: dict):
    found_match = sigindex.get(dest_sig)
    if found_match is None:
        # Unknown, ambiguous or not initialised yet.
        found_match = index_signature(dest_sig, signatures, sigindex)
    return found_match
