"""Wrapper for newrand.h"""

cdef extern from "newrand.h":
    void set_seed(unsigned int)
    unsigned int bounded_rand_int(unsigned int)


def set_seed_wrap(unsigned int custom_seed):
    set_seed(custom_seed)


def bounded_rand_int_wrap(unsigned int range_):
    return bounded_rand_int(range_)
