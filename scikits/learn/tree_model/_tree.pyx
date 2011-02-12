# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
#
# Modified from Milk
# Copyright (C) 2010, Luis Pedro Coelho <lpc@cmu.edu>
# License: MIT

import numpy as np
import sys
from time import time

cimport numpy as np
cimport cython

#def myfunc(np.ndarray[np.float64_t, ndim=2] A):
def set_entropy(np.ndarray[np.int_t, ndim=1] data):
    '''Set entropy function
    '''
    cdef Py_ssize_t j
    cdef unsigned int N, clen
    cdef int value
    N = len(data)
    clen = np.max(data)+1
    cdef np.ndarray[np.int_t, ndim=1] counts = np.zeros((clen), dtype=np.int)
    for j from 0 <= j < N:
        value = data[j]
        if (value >= clen):
            raise RuntimeError("tree_model._tree.set_entropy: label value too large. aborting");
        counts[value] += 1.;

    '''
    Here is the formula:
    
    H = - \sum px \log(px)
      = - \sum (cx/N) \log( cx / N)
      = - 1/N \sum { cx [ \log cx - \log N ] }
      = - 1/N { (\sum cx \log cx ) - ( \sum cx \log N ) }
      = - 1/N { (\sum cx \log cx ) - N \log N }
      = ( - 1/N \sum cx \log cx ) + \log N
    '''
    
    cdef double entropy = 0.
    cdef double cx
    for j from 0 <= j < clen:
        cx = counts[j]
        if cx:
            entropy += cx * np.log(cx);
    entropy /= -N;
    entropy += np.log(N);
    return entropy, counts

'''
double set_entropy(PyArrayObject* points, double* counts, long clen) {
    const int N = PyArray_SIZE(points);
#define HANDLE(type_marker, type) \
    if (PyArray_TYPE(points) == type_marker) return set_entropy<type>(reinterpret_cast<const type*>(PyArray_DATA(points)), N, counts, clen);

    HANDLE(NPY_INT, int)
    HANDLE(NPY_LONG, long)
#undef HANDLE
    assert(false);
    return 0.;
}

PyObject* py_set_entropy(PyObject* self, PyObject* args) {
    const char* errmsg = "Arguments were not what was expected for set_entropy.\n"
                        "This is an internal function: Do not call directly unless you know exactly what you're doing.\n";
    PyArrayObject* labels;
    PyArrayObject* counts;
    if (!PyArg_ParseTuple(args, "OO", &labels, &counts)) {
        PyErr_SetString(PyExc_RuntimeError,errmsg);
        return 0;
    }
    if (!PyArray_Check(labels) || (PyArray_TYPE(labels) != NPY_INT && PyArray_TYPE(labels) != NPY_LONG) || !PyArray_ISCONTIGUOUS(labels) ||
        !PyArray_Check(counts) || PyArray_TYPE(counts) != NPY_DOUBLE || !PyArray_ISCONTIGUOUS(counts)) {
        PyErr_SetString(PyExc_RuntimeError,errmsg);
        return 0;
    }
    double* cdata = reinterpret_cast<double*>(PyArray_DATA(counts));
    const int nlabels = PyArray_DIM(counts, 0);
    double res = set_entropy(labels, cdata, nlabels);
    return PyFloat_FromDouble(res);
}

PyObject* py_information_gain(PyObject* self, PyObject* args) {
    const char* errmsg = "Arguments were not what was expected for information_gain.\n"
                        "This is an internal function: Do not call directly unless you know exactly what you're doing.\n";
    PyArrayObject* labels0;
    PyArrayObject* labels1;
    if (!PyArg_ParseTuple(args, "OO", &labels0, &labels1)) {
        PyErr_SetString(PyExc_RuntimeError,errmsg);
        return 0;
    }
    if (!PyArray_Check(labels0) || (PyArray_TYPE(labels0) != NPY_INT && PyArray_TYPE(labels0) != NPY_LONG) || !PyArray_ISCONTIGUOUS(labels0) ||
        !PyArray_Check(labels1) || (PyArray_TYPE(labels1) != NPY_INT && PyArray_TYPE(labels1) != NPY_LONG) || !PyArray_ISCONTIGUOUS(labels1)) {
        PyErr_SetString(PyExc_RuntimeError,errmsg);
        std::cout << PyArray_TYPE(labels0) << " " << PyArray_TYPE(labels1) << '\n';
        return 0;
    }
    double* counts;
    double counts_array[8];
    long clen = 0;
#define GET_MAX(index) \
        const int N ## index = PyArray_DIM(labels ## index, 0); \
        { \
            for (int i = 0; i != N ## index; ++i) { \
                long val = 0; \
                if (PyArray_TYPE(labels ## index) == NPY_INT) val = *reinterpret_cast<const int*>(PyArray_GETPTR1(labels ## index, i)); \
                else if (PyArray_TYPE(labels ## index) == NPY_LONG) val = *reinterpret_cast<const int*>(PyArray_GETPTR1(labels ## index, i)); \
                if (val > clen) clen = val; \
            } \
        }
    GET_MAX(0);
    GET_MAX(1);
    ++clen;
    if (clen > 8) {
        counts = new(std::nothrow) double[clen];
        if (!counts) {
            PyErr_NoMemory();
            return 0;
        }
    } else {
        counts = counts_array;
    }
    const double N = N0 + N1;
    double H = - N0/N * set_entropy(labels0, counts, clen) - N1/N * set_entropy(labels1, counts, clen);
    if (clen > 8) delete [] counts;
    return PyFloat_FromDouble(H);
}

PyMethodDef methods[] = {
  {"set_entropy", py_set_entropy, METH_VARARGS , "Do NOT call directly.\n" },
  {"information_gain", py_information_gain, METH_VARARGS , "Do NOT call directly.\n" },
  {NULL, NULL,0,NULL},
};

const char  * module_doc =
    "Internal Module.\n"
    "\n"
    "Do NOT use directly!\n";

} // namespace

extern "C"
void init_tree()
  {
    import_array();
    (void)Py_InitModule3("_tree", methods, module_doc);
  }

'''

