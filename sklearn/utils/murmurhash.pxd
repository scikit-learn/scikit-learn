"""Export fast murmurhash C/C++ routines + cython wrappers"""

cdef extern from "MurmurHash3.h":
    void MurmurHash3_x86_32(void* key, int len, unsigned int seed,
                            void* out)

    void MurmurHash3_x86_128(void* key, int len, unsigned int seed,
                             void* out)

    void MurmurHash3_x64_128(void* key, int len, unsigned int seed,
                             void* out)


cpdef unsigned int murmurhash3_int_uint(int key, unsigned int seed)


cpdef int murmurhash3_int_int(int key, unsigned int seed)


cpdef unsigned int murmurhash3_bytes_uint(bytes key, unsigned int seed)


cpdef int murmurhash3_bytes_int(bytes key, unsigned int seed)
