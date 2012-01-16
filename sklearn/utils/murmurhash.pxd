cdef extern from "MurmurHash3.h":
    void MurmurHash3_x86_32(void* key, int len, unsigned int seed,
                            void* out)
