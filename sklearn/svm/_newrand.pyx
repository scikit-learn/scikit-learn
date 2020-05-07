"""
Wrapper for newrand.h

"""

cdef extern from "newrand-cache.h":
	void set_seed(int)
	int bounded_rand_int(int)

def set_seed_wrap(unsigned custom_seed):
	set_seed(custom_seed)

def bounded_rand_int_wrap(int orig_range):
	return bounded_rand_int(orig_range)
