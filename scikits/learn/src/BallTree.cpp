#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <iostream>
#include <structmember.h>
#include "BallTree.h"
#include "BallTreePoint.h"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

/* define as extern functions that must be called from python */
extern "C" {
#ifdef IS_PY3K
    PyObject * PyInit_ball_tree (void);
#else
    PyMODINIT_FUNC initball_tree(void);
#endif
}

/************************************************************
 * BallTree.cpp
 *  the following code defines a BallTree python class which
 *  wraps BallTree<double> in BallTree.h
 *
 *  Currently, all python arrays are converted to aligned
 *  double arrays.  If the input object is not of type double
 *  or not aligned, then an internal copy is made.
 ************************************************************/

//define the BallTree object
typedef struct{
    PyObject_HEAD
        BallTree<BallTree_Point>* tree;
    int size;
    int dim;
    std::vector<BallTree_Point*> *Points;
    PyObject* data;
} BallTreeObject;

//deallocation of BallTree Object
static void
    BallTree_dealloc(BallTreeObject* self)
{
    if( self->Points != NULL ){
        for(size_t i=0; i<self->Points->size(); i++)
            delete self->Points->at(i);
        self->Points->resize(0);
        delete self->Points;
    }
    if(self->tree != NULL)
        delete self->tree;

    if( self->data != NULL)
        Py_DECREF(self->data);

    self->Points = 0;
    self->tree = 0;
}

//BallTree.__new__ will use the default

//initialization of BallTree object
// argument is a single array of size [D,N]
static int
BallTree_init(BallTreeObject *self, PyObject *args, PyObject *kwds){
  //we use goto statements : all variables should be declared up front
    PyObject *arg=NULL;
    PyObject *arr=NULL;
    long int leaf_size=20;

    static char *kwlist[] = {"x", "leafsize", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|l", kwlist,
        &arg,&leaf_size) )
        goto fail;

    if(leaf_size <= 0){
        PyErr_SetString(PyExc_ValueError,
            "BallTree : leaf size must be greater than zero");
        goto fail;
    }

  //view this object as an array of doubles,
    arr = PyArray_FROM_OTF(arg,NPY_DOUBLE,0);
    if(arr==NULL)
        goto fail;

  //check that it is 2D
    if( PyArray_NDIM(arr) != 2)
        goto fail;

    if(self != NULL){
    //create the list of points
        self->size = PyArray_DIMS(arr)[0];
        self->dim = PyArray_DIMS(arr)[1];

        int inc = PyArray_STRIDES(arr)[1]/PyArray_DESCR(arr)->elsize;
        self->Points = new std::vector<BallTree_Point*>(self->size);
        for(int i=0;i<self->size;i++)
            self->Points->at(i) = new BallTree_Point(arr,
            (double*)PyArray_GETPTR2(arr,i,0),
            inc, PyArray_DIM(arr,1));
        self->tree = new BallTree<BallTree_Point>(*(self->Points),
            leaf_size);
    }
    self->data = arr;
  //Py_DECREF(arr);
    return 0;

    fail:
    Py_XDECREF(arr);
    return -1;
}

//query the ball tree.  Arguments are the array of search points
// and the number of nearest neighbors, k (optional)
static PyObject *
BallTree_query(BallTreeObject *self, PyObject *args, PyObject *kwds){
  //we use goto statements : all variables should be declared up front
    int return_distance = 1;
    long int k = 1;
    PyObject *arg = NULL;
    PyObject *arr = NULL;
    PyObject *nbrs = NULL;
    PyObject *dist = NULL;
    PyArrayIterObject *arr_iter = NULL;
    PyArrayIterObject *nbrs_iter = NULL;
    PyArrayIterObject *dist_iter = NULL;
    long int* nbrs_data;
    double* dist_data;
    static char *kwlist[] = {"x", "k", "return_distance", NULL};

    int nd, pt_size, pt_inc;
    npy_intp* dim;

  //parse arguments.  If k is not provided, the default is 1
    if(!PyArg_ParseTupleAndKeywords(args,kwds,"O|li",kwlist,
    &arg,&k,&return_distance)){
        goto fail;
    }

  //check value of k
    if(k < 1){
        PyErr_SetString(PyExc_ValueError,
            "k must be positive");
        goto fail;
    }

    if(k > self->size){
        PyErr_SetString(PyExc_ValueError,
            "k must not be greater than number of points");
        goto fail;
    }

  //get the array object from the first argument
    arr = PyArray_FROM_OTF(arg,NPY_DOUBLE,NPY_ALIGNED);

  //check that the array was properly constructed
    if(arr==NULL){
        PyErr_SetString(PyExc_ValueError,
            "pt must be convertable to array");
        goto fail;
    }

    nd = PyArray_NDIM(arr);
    if(nd == 0){
        PyErr_SetString(PyExc_ValueError,
            "pt cannot be zero-sized array");
        goto fail;
    }
    pt_size = PyArray_DIM(arr,nd-1);
    if( pt_size != self->tree->PointSize() ){
        PyErr_SetString(PyExc_ValueError,
            "points are incorrect dimension");
        goto fail;
    }

  //create a neighbors array and distance array
    dim = new npy_intp[nd];
    for(int i=0; i<nd-1;i++)
        dim[i] = PyArray_DIM(arr,i);
    dim[nd-1] = k;
    nbrs = (PyObject*)PyArray_SimpleNew(nd,dim,PyArray_LONG);
    if(return_distance)
        dist = (PyObject*)PyArray_SimpleNew(nd,dim,PyArray_DOUBLE);
    delete[] dim;

    if(nbrs==NULL)
        goto fail;

  //create iterators to cycle through points
    nd-=1;
    arr_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(arr,&nd);
    nbrs_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(nbrs,&nd);
    if(return_distance)
        dist_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(dist,&nd);
    nd+=1;

    if( arr_iter==NULL || nbrs_iter==NULL ||
        (arr_iter->size != nbrs_iter->size) ||
        (return_distance &&
    (dist_iter==NULL || (arr_iter->size != dist_iter->size))) ){
        PyErr_SetString(PyExc_ValueError,
            "failure constructing iterators");
        goto fail;
    }

    pt_inc = PyArray_STRIDES(arr)[nd-1] / PyArray_DESCR(arr)->elsize;
    if( (PyArray_STRIDES(nbrs)[nd-1] != PyArray_DESCR(nbrs)->elsize ) ||
        (return_distance &&
    (PyArray_STRIDES(dist)[nd-1] != PyArray_DESCR(dist)->elsize )) ){
        PyErr_SetString(PyExc_ValueError,
            "nbrs & dist not allocated as a C-array");
        goto fail;
    }

  //iterate through points and determine neighbors
  //warning: if nbrs is not a C-array, or if we're not iterating
  // over the last dimension, this may cause a seg fault.
    if(return_distance){
        while(arr_iter->index < arr_iter->size){
            BallTree_Point pt(arr,
                (double*)PyArray_ITER_DATA(arr_iter),
                pt_inc,pt_size);
            nbrs_data = (long int*)(PyArray_ITER_DATA(nbrs_iter));
            dist_data = (double*)(PyArray_ITER_DATA(dist_iter));
            self->tree->query(pt,k,nbrs_data,dist_data);
            PyArray_ITER_NEXT(arr_iter);
            PyArray_ITER_NEXT(nbrs_iter);
            PyArray_ITER_NEXT(dist_iter);
        }
    }else{
        while(arr_iter->index < arr_iter->size){
            BallTree_Point pt(arr,
                (double*)PyArray_ITER_DATA(arr_iter),
                pt_inc,pt_size);
            nbrs_data = (long int*)(PyArray_ITER_DATA(nbrs_iter));
            self->tree->query(pt,k,nbrs_data);
            PyArray_ITER_NEXT(arr_iter);
            PyArray_ITER_NEXT(nbrs_iter);
        }
    }

    if(return_distance){
        Py_DECREF(arr_iter);
        Py_DECREF(nbrs_iter);
        Py_DECREF(dist_iter);
        Py_DECREF(arr);

        arr = Py_BuildValue("(OO)",dist,nbrs);
        Py_DECREF(nbrs);
        Py_DECREF(dist);
        return arr;

    }else{
        Py_DECREF(arr_iter);
        Py_DECREF(nbrs_iter);
        Py_DECREF(arr);
        return nbrs;
    }

    fail:
    Py_XDECREF(arr);
    Py_XDECREF(nbrs);
    Py_XDECREF(dist);
    Py_XDECREF(arr_iter);
    Py_XDECREF(nbrs_iter);
    Py_XDECREF(dist_iter);
    return NULL;
}

//query the ball tree.  Arguments are the array of search points
// and the radius around each point to search
static PyObject *
BallTree_queryball(BallTreeObject *self, PyObject *args, PyObject *kwds){
  //we use goto statements : all variables should be declared up front
    int count_only = 0;
    double r;
    PyObject *arg = NULL;
    PyObject *arr = NULL;
    PyObject *nbrs = NULL;
    PyArrayIterObject* arr_iter = NULL;
    PyArrayIterObject* nbrs_iter = NULL;
    int nd, pt_size, pt_inc;
    static char *kwlist[] = {"x", "r", "count_only", NULL};

  //parse arguments.  If kmax is not provided, the default is 20
    if(!PyArg_ParseTupleAndKeywords(args,kwds,"Od|i",
                                    kwlist,&arg,&r,&count_only)){
        goto fail;
    }

  //check value of r
    if(r < 0){
        PyErr_SetString(PyExc_ValueError,
            "r must not be negative");
        goto fail;
    }

  //get the array object from the first argument
    arr = PyArray_FROM_OTF(arg,NPY_DOUBLE,NPY_ALIGNED);

  //check that the array was properly constructed
    if(arr==NULL){
        PyErr_SetString(PyExc_ValueError,
            "pt must be convertable to array");
        goto fail;
    }

    nd = PyArray_NDIM(arr);
    if(nd == 0){
        PyErr_SetString(PyExc_ValueError,
            "pt cannot be zero-sized array");
        goto fail;
    }
    pt_size = PyArray_DIM(arr,nd-1);
    if( pt_size != self->tree->PointSize() ){
        PyErr_SetString(PyExc_ValueError,
            "points are incorrect dimension");
        goto fail;
    }


  // Case 1: return arrays of all neighbors for each point
  //
    if(!count_only){
    //create a neighbors array.  This is an array of python objects.
    // each of which will be a numpy array of neighbors
        nbrs = (PyObject*)PyArray_SimpleNew(nd-1,PyArray_DIMS(arr),
                                            PyArray_OBJECT);

        if(nbrs==NULL){
            goto fail;
        }

    //create iterators to cycle through points
        --nd;
        arr_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(arr,&nd);
        nbrs_iter = (PyArrayIterObject*)PyArray_IterNew(nbrs);
        ++nd;

        if( arr_iter==NULL || nbrs_iter==NULL ||
        (arr_iter->size != nbrs_iter->size)){
            PyErr_SetString(PyExc_ValueError,
                "unable to construct iterator");
            goto fail;
        }


        pt_inc = PyArray_STRIDES(arr)[nd-1] / PyArray_DESCR(arr)->elsize;
        if(PyArray_NDIM(nbrs)==0){
            BallTree_Point pt(arr,
                (double*)PyArray_ITER_DATA(arr_iter),
                pt_inc,pt_size);
            std::vector<long int> nbrs_vec;
            self->tree->query_ball(pt,r,nbrs_vec);
            npy_intp N_nbrs = nbrs_vec.size();
            PyObject* nbrs_obj = PyArray_SimpleNew(1, &N_nbrs, PyArray_LONG);
            long int* data = (long int*)PyArray_DATA(nbrs_obj);
            for(int i=0; i<N_nbrs; i++)
                data[i] = nbrs_vec[i];
            PyObject* tmp = nbrs;
            nbrs = nbrs_obj;
            Py_DECREF(tmp);
        }else{
            while(arr_iter->index < arr_iter->size){
                BallTree_Point pt(arr,
                    (double*)PyArray_ITER_DATA(arr_iter),
                    pt_inc,pt_size);
                std::vector<long int> nbrs_vec;
                self->tree->query_ball(pt,r,nbrs_vec);

                npy_intp N_nbrs = nbrs_vec.size();

                PyObject* nbrs_obj = PyArray_SimpleNew(1, &N_nbrs,
                                                       PyArray_LONG);
                long int* data = (long int*)PyArray_DATA(nbrs_obj);
                for(int i=0; i<N_nbrs; i++)
                    data[i] = nbrs_vec[i];

                PyObject** nbrs_data = (PyObject**)PyArray_ITER_DATA(nbrs_iter);
                PyObject* tmp = nbrs_data[0];
                nbrs_data[0] = nbrs_obj;
                Py_XDECREF(tmp);

                PyArray_ITER_NEXT(arr_iter);
                PyArray_ITER_NEXT(nbrs_iter);
            }
        }
    }

  // Case 2 : return number of neighbors for each point
    else{
    //create an array to keep track of the count
        nbrs = (PyObject*)PyArray_SimpleNew(nd-1,PyArray_DIMS(arr),
                                            PyArray_LONG);
        if(nbrs==NULL){
            goto fail;
        }

    //create iterators to cycle through points
        --nd;
        arr_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(arr,&nd);
        nbrs_iter = (PyArrayIterObject*)PyArray_IterNew(nbrs);
        ++nd;

        if( arr_iter==NULL || nbrs_iter==NULL ||
            (arr_iter->size != nbrs_iter->size)){
            PyErr_SetString(PyExc_ValueError,
                "unable to construct iterator");
            goto fail;
        }

    //go through points and call BallTree::query_ball to count neighbors
        pt_inc = PyArray_STRIDES(arr)[nd-1] / PyArray_DESCR(arr)->elsize;
        while(arr_iter->index < arr_iter->size){
            BallTree_Point pt(arr,
                (double*)PyArray_ITER_DATA(arr_iter),
                pt_inc,pt_size);

            long int* nbrs_count = (long int*)PyArray_ITER_DATA(nbrs_iter);
            *nbrs_count = self->tree->query_ball(pt,r);

            PyArray_ITER_NEXT(arr_iter);
            PyArray_ITER_NEXT(nbrs_iter);
        }
    }
    Py_DECREF(nbrs_iter);
    Py_DECREF(arr_iter);
    Py_DECREF(arr);
    return nbrs;

    fail:
    Py_XDECREF(nbrs_iter);
    Py_XDECREF(arr_iter);
    Py_XDECREF(arr);
    Py_XDECREF(nbrs);
    return NULL;
}


//define the data members of BallTree
static PyMemberDef BallTree_members[] = {
    {"size", T_INT, offsetof(BallTreeObject,size),0,
    "Number of points in the Ball Tree"},
    {"dim", T_INT, offsetof(BallTreeObject,dim),0,
    "Dimension of the Ball Tree"},
    {"data", T_OBJECT_EX, offsetof(BallTreeObject,data),0,
    "View of data making up the Ball Tree"},
    {NULL} /* Sentinel */
};


//define the methods of BallTree
static PyMethodDef BallTree_methods[] = {
    {"query", (PyCFunction)BallTree_query, METH_VARARGS|METH_KEYWORDS,
     "query(x, k=1, return_distance=True)                  \n"
     "                                                     \n"
     "   query the Ball Tree for the k nearest neighbors   \n"
     "                                                     \n"
     "Parameters                                           \n"
     "----------                                           \n"
     "x : array-like, last dimension self.dim              \n"
     "      An array of points to query                    \n"
     "k : integer  (default = 1)                           \n"
     "      The number of nearest neighbors to return      \n"
     "return_distance : boolean (default = True)           \n"
     "      if True, return a tuple (d,i)                  \n"
     "      if False, return array i                       \n"
     "                                                     \n"
     "Returns                                              \n"
     "-------                                              \n"
     "i    : if return_distance == False                   \n"
     "(d,i) : if return_distance == True                   \n"
     "                                                     \n"
     "d : array of doubles - shape: x.shape[:-1] + (k,)    \n"
     "    each entry gives the list of distances to the    \n"
     "    neighbors of the corresponding point             \n"
     "    (note that distances are not sorted)             \n"
     "                                                     \n"
     "i : array of integers - shape: x.shape[:-1] + (k,)   \n"
     "    each entry gives the list of indices of          \n"
     "    neighbors of the corresponding point             \n"
     "    (note that neighbors are not sorted)             \n"
    },
    {"query_ball", (PyCFunction)BallTree_queryball, METH_VARARGS|METH_KEYWORDS,
     "query_ball(x,r,count_only = False)                   \n"
     "                                                     \n"
     "    query the Ball Tree for the k nearest neighbors  \n"
     "                                                     \n"
     " Parameters                                          \n"
     " ----------                                          \n"
     " x : array-like, last dimension self.dim             \n"
     "       An array of points to query                   \n"
     " r : floating-point value                            \n"
     "       Radius around each point within which all     \n"
     "       neighbors are returned                        \n"
     " count_only : boolean (default = False)              \n"
     "                if True,  return count of neighbors  \n"
     "                           for each point            \n"
     "                if False, return full list of        \n"
     "                           neighbors for each point  \n"
     "                                                     \n"
     " Returns                                             \n"
     " -------                                             \n"
     " i : array of integers, shape: x.shape[:-1]          \n"
     "     if count_only is False each entry gives the     \n"
     "     list of neighbors of the corresponding point    \n"
     "     (note that neighbors are not sorted). Otherwise \n"
     "     return only the number of neighbors.            \n"
    },
    {NULL}  /* Sentinel */
};

//construct the BallTree type from the preceeding functions
static PyTypeObject BallTreeType = {
#if defined (IS_PY3K)
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
#endif
    "BallTree.BallTree",       /*tp_name*/
    sizeof(BallTreeObject),  /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)BallTree_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
#if defined (IS_PY3K)
    (void *) 0,                /* tp_reserved */
#else
    0,                         /*tp_compare*/
#endif
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Ball Tree for fast nearest-neighbor searches :           \n"
    "                                                         \n"
    " BallTree(M, leafsize=20)                                \n"
    "                                                         \n"
    "                                                         \n"
    "   Parameters                                            \n"
    "   ----------                                            \n"
    "   M : array-like, shape = [N,D]                         \n"
    "        N is the number of points in the data set, and   \n"
    "        D is the dimension of the parameter space.       \n"
    "       Note: if M is an aligned array of doubles (not    \n"
    "        necessarily contiguous) then data will not be    \n"
    "        copied. Otherwise, an internal copy will be made.\n"
    "   leafsize : positive integer (default = 20)           \n"
    "               number of points at which to switch       \n"
    "               to brute-force                            \n"
    "                                                         \n",/* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    BallTree_methods,          /* tp_methods */
    BallTree_members,          /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)BallTree_init,   /* tp_init */
    0,                         /* tp_alloc */
    0,                         /* tp_new */
};


//Wrapper for Brute-Force neighbor search
static PyObject*
BallTree_knn_brute(PyObject *self, PyObject *args, PyObject *kwds){
    long int k = 1;
    std::vector<BallTree_Point*> Points;

    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arr1 = NULL;
    PyObject *arr2 = NULL;
    PyObject *nbrs = NULL;

    long int* nbrs_data;
    PyArrayIterObject *arr2_iter = NULL;
    PyArrayIterObject *nbrs_iter = NULL;
    static char *kwlist[] = {"x", "pt", "k", NULL};

    npy_intp* dim;
    int nd, pt_size, pt_inc;
    long int N;
    long int D;

  //parse arguments.  If k is not provided, the default is 1
    if(!PyArg_ParseTupleAndKeywords(args,kwds,"OO|l",kwlist,
        &arg1,&arg2,&k))
        goto fail;

  //First array should be a 2D array of doubles
    arr1 = PyArray_FROM_OTF(arg1,NPY_DOUBLE,0);
    if(arr1==NULL)
        goto fail;
    if( PyArray_NDIM(arr1) != 2){
        PyErr_SetString(PyExc_ValueError,
            "x must be two dimensions");
        goto fail;
    }

  //Second array should be a 1D array of doubles
    arr2 = PyArray_FROM_OTF(arg2,NPY_DOUBLE,0);
    if(arr2==NULL)
        goto fail;

    nd = PyArray_NDIM(arr2);

    if(nd == 0){
        PyErr_SetString(PyExc_ValueError,
            "pt cannot be zero-sized array");
        goto fail;
    }
    pt_size = PyArray_DIM(arr2,nd-1);

  //Check that dimensions match
    N = PyArray_DIMS(arr1)[0];
    D = PyArray_DIMS(arr1)[1];
    if( pt_size != D ){
        PyErr_SetString(PyExc_ValueError,
            "pt must be same dimension as x");
        goto fail;
    }

  //check the value of k
    if(k<1){
        PyErr_SetString(PyExc_ValueError,
            "k must be a positive integer");
        goto fail;
    }

    if(k>N){
        PyErr_SetString(PyExc_ValueError,
            "k must be less than the number of points");
        goto fail;
    }

  //create a neighbors array and distance array
    dim = new npy_intp[nd];
    for(int i=0; i<nd-1;i++)
        dim[i] = PyArray_DIM(arr2,i);
    dim[nd-1] = k;
    nbrs = (PyObject*)PyArray_SimpleNew(nd,dim,PyArray_LONG);
    delete[] dim;

    if(nbrs==NULL)
        goto fail;

  //create iterators to cycle through points
    nd-=1;
    arr2_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(arr2,&nd);
    nbrs_iter = (PyArrayIterObject*)PyArray_IterAllButAxis(nbrs,&nd);
    nd+=1;

    if( arr2_iter==NULL ||
        nbrs_iter==NULL ||
        (arr2_iter->size != nbrs_iter->size) ){
        PyErr_SetString(PyExc_ValueError,
            "failure constructing iterators");
        goto fail;
    }

    pt_inc = PyArray_STRIDES(arr2)[nd-1] / PyArray_DESCR(arr2)->elsize;
    if(PyArray_STRIDES(nbrs)[nd-1] != PyArray_DESCR(nbrs)->elsize ){
        PyErr_SetString(PyExc_ValueError,
            "nbrs not allocated as a C-array");
        goto fail;
    }

  //create the list of points
    pt_inc = PyArray_STRIDES(arr1)[1]/PyArray_DESCR(arr1)->elsize;
    Points.resize(N);
    for(int i=0;i<N;i++)
        Points[i] = new BallTree_Point(arr1,
        (double*)PyArray_GETPTR2(arr1,i,0),
        pt_inc, PyArray_DIM(arr1,1));

  //iterate through points and determine neighbors
  //warning: if nbrs is not a C-array, or if we're not iterating
  // over the last dimension, this may cause a seg fault.
    while(arr2_iter->index < arr2_iter->size){
        BallTree_Point Query_Point(arr2,
            (double*)PyArray_ITER_DATA(arr2_iter),
            pt_inc,pt_size);
        nbrs_data = (long int*)(PyArray_ITER_DATA(nbrs_iter));
        BruteForceNeighbors(Points, Query_Point,
            k, nbrs_data );
        PyArray_ITER_NEXT(arr2_iter);
        PyArray_ITER_NEXT(nbrs_iter);
    }

    for(int i=0;i<N;i++)
        delete Points[i];

    return nbrs;

    fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
    Py_XDECREF(nbrs);
    Py_XDECREF(arr2_iter);
    Py_XDECREF(nbrs_iter);
    return NULL;
}


//module methods: none thus far.
static PyMethodDef BTmodule_methods[] = {
  {"knn_brute", (PyCFunction)BallTree_knn_brute,
   METH_VARARGS|METH_KEYWORDS,
   "knn_brute(x, pt, k=1)                                     \n"
   "                                                          \n"
   "  Brute-Force k-nearest neighbor search.                  \n"
   "                                                          \n"
   "  Parameters                                              \n"
   "  ----------                                              \n"
   "  x  : array of shape [N,D]                               \n"
   "       representing N points in D dimensions              \n"
   "  pt : array-like, last dimension D                       \n"
   "       An array of points to query                        \n"
   "  k  : a positive integer, giving the number of nearest   \n"
   "       neighbors to query                                 \n"
   "                                                          \n"
   "  Returns                                                 \n"
   "  -------                                                 \n"
   "  nbrs : array of integers - shape: pt.shape[:-1] + (k,)  \n"
   "           each entry gives the list of indices of        \n"
   "           neighbors of the corresponding point           \n"},
  {NULL, NULL, 0, NULL}        /* End of Methods */
};

//initialize the module

#define MODULE_DOC              \
"Ball Tree package                                   \n"  \
" Written by Jake VanderPlas, January 2010           \n"  \
"   vanderplas@astro.washington.edu                  \n"  \
"   http://www.astro.washington.edu/users/vanderplas \n"  \
"                                                    \n"  \
" A Ball Tree is a data structure which can be used  \n"  \
"  to perform fast neighbor searches in data sets of \n"  \
"  very high dimensionality.  For low dimensional    \n"  \
"  problems (dimension less than 5-10) a KD tree is  \n"  \
"  a better choice (see, e.g. scipy.spatial.cKDTree) \n"  \
"                                                    \n"  \
" This package also provides an optimized brute-force\n"  \
"  neighbor search (knn_brute) which has better      \n"  \
"  performance than either tree algorithm for smaller\n"  \
"  data-sets (number of points less than ~1000),     \n"  \
"  especially when querying for more than one nearest\n"  \
"  neighbor.                                         \n"  \

struct module_state {
    PyObject *error;
};

#ifdef IS_PY3K

/* Some python3 specific methods */

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int ball_tree_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int ball_tree_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef BTmodule_def = {
    PyModuleDef_HEAD_INIT,
    "ball_tree",
    MODULE_DOC,
    sizeof (struct module_state),
    BTmodule_methods,
    NULL,
    ball_tree_traverse,
    ball_tree_clear,
    NULL
};

#define INITERROR return NULL

PyObject * PyInit_ball_tree (void)
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#define INITERROR return

PyMODINIT_FUNC initball_tree(void)
#endif
{
    PyObject* m;
    BallTreeType.tp_new = PyType_GenericNew;

#if defined (IS_PY3K)
    if (PyType_Ready(&BallTreeType) < 0) return NULL;
    m = PyModule_Create(&BTmodule_def);
    if (m==NULL) return NULL;
#else
    if (PyType_Ready(&BallTreeType) < 0) return;
    m = Py_InitModule3("ball_tree", BTmodule_methods, MODULE_DOC);
    if (m==NULL) return;
#endif

    Py_INCREF(&BallTreeType);
    PyModule_AddObject(m,"BallTree", (PyObject*)&BallTreeType);
    import_array();

#if defined (IS_PY3K)
    return m;
#endif
}
