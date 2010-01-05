/*
 * Last Change: Wed Dec 06 08:00 PM 2006 J
 */
#include <Python.h>
#include <numpy/arrayobject.h>

#include <math.h>

PyObject* compute_ss_frame_1d_py(PyObject* dum, PyObject *arg);
PyObject* compute_em_frame_1d_py(PyObject* dum, PyObject *arg);

/*
 * Pure C methods
 */
static int compute_ss_frame_1d(double sample, int nc, double* w, double* mu, double *va,
        double *cx, double *cxx, double nu);
static int update_em_frame_1d(const double* w, const double* cx, const double *cxx,
        int nc, double *cmu, double *cva);

static PyMethodDef mymethods[] = {
    {"compute_ss_frame_1d", compute_ss_frame_1d_py, METH_VARARGS, ""},
    {"compute_em_frame_1d", compute_em_frame_1d_py, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/*
 * function table
 */
PyMODINIT_FUNC init_rawden(void)
{
    (void)Py_InitModule("_rawden", mymethods);
    import_array();
}

PyObject* compute_ss_frame_1d_py(PyObject* self, PyObject *args)
{
    PyObject    *w, *mu, *va, *cx, *cxx;
    PyObject    *w_a, *mu_a, *va_a, *cx_a, *cxx_a;
    double      nu, sample;
    npy_intp    ndim, *dims, len;
    double      *mu_ca, *va_ca, *w_ca, *cx_ca, *cxx_ca;

    const npy_intp  mrank   = 1;

    int         st;

    /*
     * Init python object holder to NULL so that we can use Py_XDECREF on all 
     * the objets when something is woring
     */
    w   = NULL;
    mu  = NULL;
    va  = NULL;
    cx  = NULL;
    cxx = NULL;

    w_a     = NULL;
    mu_a    = NULL;
    va_a    = NULL;
    cx_a    = NULL;
    cxx_a   = NULL;
    /*
     * Parsing of args: w, cx and cxx are inout
     */
    if (!PyArg_ParseTuple(args, "dOOOOOd",  &sample, &w, &mu, &va, &cx, &cxx, &nu)){
        return NULL;
    }
    if ( (nu > 1) | (nu <= 0) ) {
        PyErr_SetString(PyExc_TypeError, "nu should be between 0 and 1");
        return NULL;
    }

    /* inout entries */
    w_a     = PyArray_FROM_OTF(w, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (w_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "w array not convertible");
        return NULL;
    }

    cx_a    = PyArray_FROM_OTF(cx, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (cx_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "cx array not convertible");
        goto fail;
    }

    cxx_a   = PyArray_FROM_OTF(cxx, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (cxx_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "cxx array not convertible");
        goto fail;
    }

    /* in only entries */
    mu_a    = PyArray_FROM_OTF(mu, NPY_DOUBLE, NPY_IN_ARRAY);
    if (mu_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "mu array not convertible");
        goto fail;
    }

    va_a    = PyArray_FROM_OTF(va, NPY_DOUBLE, NPY_IN_ARRAY);
    if (va_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "va array not convertible");
        goto fail;
    }

    /*
     * Check that in and out have same size and same rank
     */
    ndim    = PyArray_NDIM(w_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "w rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(cx_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "cx rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(cxx_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "cxx rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(mu_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "mu rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(va_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "va rank should be 1");
        goto fail;
    }

    dims    = PyArray_DIMS(w_a);
    len     = dims[0];
    //fprintf(stderr, "%s:%s, len is %d\n", __FILE__, __func__, len);
    dims    = PyArray_DIMS(cx_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "cx shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(cxx_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "cxx shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(mu_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "mu_a shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(va_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "va_a shape should match !");
        goto fail;
    }

    /* 
     * Get pointer to the data
     */
    w_ca    = PyArray_DATA(w_a);
    if (w_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for w_ca");
        goto fail;
    }
    cx_ca   = PyArray_DATA(cx_a);
    if (cx_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for cx_ca");
        goto fail;
    }
    cxx_ca  = PyArray_DATA(cxx_a);
    if (w_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for cxx_ca");
        goto fail;
    }
    mu_ca   = PyArray_DATA(mu_a);
    if (mu_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for mu_ca");
        goto fail;
    }
    va_ca   = PyArray_DATA(va_a);
    if (va_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for va_ca");
        goto fail;
    }
    /*
     * Call actual implementation
     */
    st  = compute_ss_frame_1d(sample, len, w_ca, mu_ca, va_ca, cx_ca, cxx_ca, nu);
    if (st) {
        PyErr_SetString(PyExc_TypeError, "Error while calling multi_gauss....");
        goto fail;
    }

    Py_DECREF(w_a);
    Py_DECREF(cx_a);
    Py_DECREF(cxx_a);
    Py_DECREF(mu_a);
    Py_DECREF(va_a);

    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(w_a);
    Py_XDECREF(cx_a);
    Py_XDECREF(cxx_a);
    Py_XDECREF(mu_a);
    Py_XDECREF(va_a);
    return NULL;
}

PyObject* compute_em_frame_1d_py(PyObject* self, PyObject *args)
{
    PyObject    *w, *mu, *va, *cx, *cxx;
    PyObject    *w_a, *mu_a, *va_a, *cx_a, *cxx_a;
    double      nu, sample;
    npy_intp    ndim, *dims, len;
    double      *mu_ca, *va_ca, *w_ca, *cx_ca, *cxx_ca;

    const npy_intp  mrank   = 1;

    int         st;

    /*
     * Init python object holder to NULL so that we can use Py_XDECREF on all 
     * the objets when something is woring
     */
    w   = NULL;
    mu  = NULL;
    va  = NULL;
    cx  = NULL;
    cxx = NULL;

    w_a     = NULL;
    mu_a    = NULL;
    va_a    = NULL;
    cx_a    = NULL;
    cxx_a   = NULL;
    /*
     * Parsing of args: w, cx and cxx are inout
     */
    if (!PyArg_ParseTuple(args, "dOOOOOd",  &sample, &w, &mu, &va, &cx, &cxx, &nu)){
        return NULL;
    }
    if ( (nu > 1) | (nu <= 0) ) {
        PyErr_SetString(PyExc_TypeError, "nu should be between 0 and 1");
        return NULL;
    }

    /* inout entries */
    w_a     = PyArray_FROM_OTF(w, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (w_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "w array not convertible");
        return NULL;
    }

    cx_a    = PyArray_FROM_OTF(cx, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (cx_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "cx array not convertible");
        goto fail;
    }

    cxx_a   = PyArray_FROM_OTF(cxx, NPY_DOUBLE, NPY_INOUT_ARRAY);
    if (cxx_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "cxx array not convertible");
        goto fail;
    }

    /* in only entries */
    mu_a    = PyArray_FROM_OTF(mu, NPY_DOUBLE, NPY_IN_ARRAY);
    if (mu_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "mu array not convertible");
        goto fail;
    }

    va_a    = PyArray_FROM_OTF(va, NPY_DOUBLE, NPY_IN_ARRAY);
    if (va_a == NULL) {
        PyErr_SetString(PyExc_TypeError, "va array not convertible");
        goto fail;
    }

    /*
     * Check that in and out have same size and same rank
     */
    ndim    = PyArray_NDIM(w_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "w rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(cx_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "cx rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(cxx_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "cxx rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(mu_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "mu rank should be 1");
        goto fail;
    }
    ndim    = PyArray_NDIM(va_a);
    if(ndim != mrank) {
        PyErr_SetString(PyExc_TypeError, "va rank should be 1");
        goto fail;
    }

    dims    = PyArray_DIMS(w_a);
    len     = dims[0];
    //fprintf(stderr, "%s:%s, len is %d\n", __FILE__, __func__, len);
    dims    = PyArray_DIMS(cx_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "cx shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(cxx_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "cxx shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(mu_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "mu_a shape should match !");
        goto fail;
    }
    dims    = PyArray_DIMS(va_a);
    if(dims[0] != len) {
        PyErr_SetString(PyExc_TypeError, "va_a shape should match !");
        goto fail;
    }

    /* 
     * Get pointer to the data
     */
    w_ca    = PyArray_DATA(w_a);
    if (w_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for w_ca");
        goto fail;
    }
    cx_ca   = PyArray_DATA(cx_a);
    if (cx_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for cx_ca");
        goto fail;
    }
    cxx_ca  = PyArray_DATA(cxx_a);
    if (w_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for cxx_ca");
        goto fail;
    }
    mu_ca   = PyArray_DATA(mu_a);
    if (mu_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for mu_ca");
        goto fail;
    }
    va_ca   = PyArray_DATA(va_a);
    if (va_ca == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unknown Error for va_ca");
        goto fail;
    }
    /*
     * Call actual implementation
     */
    st  = compute_ss_frame_1d(sample, len, w_ca, mu_ca, va_ca, cx_ca, cxx_ca, nu);
    if (st) {
        PyErr_SetString(PyExc_TypeError, "Error while calling multi_gauss....");
        goto fail;
    }
    st  = update_em_frame_1d(w_ca, cx_ca, cxx_ca, len, mu_ca, va_ca);
    if (st) {
        PyErr_SetString(PyExc_TypeError, "Error while calling update_em_frame_1d....");
        goto fail;
    }

    Py_DECREF(w_a);
    Py_DECREF(cx_a);
    Py_DECREF(cxx_a);
    Py_DECREF(mu_a);
    Py_DECREF(va_a);

    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(w_a);
    Py_XDECREF(cx_a);
    Py_XDECREF(cxx_a);
    Py_XDECREF(mu_a);
    Py_XDECREF(va_a);
    return NULL;
}

int compute_ss_frame_1d(double sample, int nc, double* w, double* mu, double *va,
        double *cx, double *cxx, double nu)
{
    /*
     * TODO: check va division
     */
    int     i;
    double  inva, fac, *gam, acc;

    gam = malloc(sizeof(*gam) * nc);
    if (gam == NULL) {
        return -1;
    }

    /* Compute gamma */
    acc = 0;
    for (i = 0; i < nc; ++i) {
        inva        = 1/va[i];
        fac         = 1 / sqrt(2 * M_PI * va[i]);
        gam[i]      = fac * exp( -0.5 * inva * (sample - mu[i]) * (sample - mu[i]));
        gam[i]      *= w[i];
        acc         += gam[i];
    }
    /* Normalize gamma */
    for (i = 0; i < nc; ++i) {
        gam[i]  /= acc;
    }

    /* Compute new SS from EM (cx and cxx) */
    for (i = 0; i < nc; ++i) {
        w[i]  *= (1 - nu);
        w[i]  += nu * gam[i];
        cx[i]   = (1 - nu) * cx[i] + nu * sample * gam[i]; 
        cxx[i]  = (1 - nu) * cxx[i] + nu * sample * sample * gam[i]; 
    }

    free(gam);

    return 0;
}

/*
 * update mu and va from SS w, cx and cxx. Only mu and va are modified
 * all arrays have same length nc
 */
int update_em_frame_1d(const double* cw, const double* cx, const double *cxx,
        int nc, double *cmu, double *cva)
{
    /*
     * TODO: check va division
     */
    int     i;
    double  invw;

    /* Compute new SS from EM (cx and cxx) */
    for (i = 0; i < nc; ++i) {
        invw    = 1/cw[i];
        cmu[i]  = cx[i] * invw;
        cva[i]  = cxx[i] * invw - cmu[i] * cmu[i]; 
    }

    return 0;
}
