cimport murmurhash


cpdef unsigned int murmurhash3_int(int key, unsigned int seed):
    cdef int out
    MurmurHash3_x86_32(&key, 1, seed, &out)
    return 0


cpdef unsigned int murmurhash3_bytes(bytes key, unsigned int seed):
    # TODO
    return 0


def murmurhash3(key, seed):
    if isinstance(key, bytes):
        return murmurhash3_bytes(key, seed)
    elif isinstance(key, unicode):
        return murmurhash3_bytes(key.encode('utf-8'), seed)
    elif isinstance(key, int):
        return murmurhash3_int(key, seed)
    else:
        raise ValueError(
            "key %r with type %s is not supported. "
            "Explicit conversion to bytes is required" % (key, type(key)))
