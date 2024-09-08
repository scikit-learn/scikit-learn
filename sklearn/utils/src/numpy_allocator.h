/* See https://docs.python.org/3/extending/extending.html */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

/* We need at least numpy 1.22 */
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION
#include "numpy/arrayobject.h"

#include <stdlib.h>  /* memory allocation */
#include <string.h>  /* memset */


static inline void *
aligned_64_malloc(void *ctx, size_t size)
{
    size_t rest = size % 64;
    if (rest > 0) {
        size += 64 - rest;  /* allocate a bit more than we need */
    }
    return aligned_alloc(64, size);
}


static inline void *
aligned_64_calloc(void *ctx, size_t nelem, size_t elsize)
{
    size_t size = nelem * elsize;
    void *x = aligned_64_malloc(ctx, size);
    if (x != NULL) {
        size_t rest = size % 64;
        if (rest > 0) {
            size += 64 - rest;
        }
        memset(x, 0, size);
    }
    return x;
}


static inline void *
std_realloc(void *ctx, void *ptr, size_t new_size)
{
    /* TODO: Do we need aligned realloc? */
    return realloc(ptr, new_size);
}

static inline void
std_free(void *ctx, void *ptr, size_t size)
{
    /* TODO: Deal with MSCV incompatibility of aligned_alloc. */
    free(ptr);
}


static PyDataMem_Handler aligned_handler = {
    "aligned_allocator",    /* name */
    1,                      /* version */
    {
        NULL,               /* ctx */
        aligned_64_malloc,  /* malloc */
        aligned_64_calloc,  /* calloc */
        std_realloc,        /* realloc */
        std_free            /* free */
    }
};


static inline void
c_set_aligned_allocation()
{
    /* Name must be "mem_handler". */
    PyObject *aligned_handler_obj = PyCapsule_New(&aligned_handler, "mem_handler", NULL);
    if (aligned_handler_obj == NULL) {
        return;
    }
    PyObject *old = PyDataMem_SetHandler(aligned_handler_obj);
    Py_DECREF(aligned_handler_obj);
};


static inline void
c_set_numpy_default_allocation()
{
    PyDataMem_SetHandler(NULL);
};
