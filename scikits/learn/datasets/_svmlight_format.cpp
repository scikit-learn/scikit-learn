/*
 * Author: Mathieu Blondel.
 *
 * This module implements _load_svmlight_format, a fast and memory efficient
 * function to load the file format originally created for svmlight and now used
 * by many other libraries, including libsvm.
 *
 * The function loads the file directly in a CSR sparse matrix without memory
 * copying.  The approach taken is to use 4 C++ vectors (data, indices, indptr
 * and labels) and to incrementally feed them with elements. Ndarrays are then
 * instanciated by PyArray_SimpleNewFromData, i.e., no memory is
 * copied.
 *
 * Since the memory is not allocated by the ndarray, the ndarray does own the
 * memory and thus cannot deallocate it. To automatically deallocate memory, the
 * technique described at http://blog.enthought.com/?p=62 is used. The main idea
 * is to use an additional object that the ndarray does own and that will be
 * responsible for deallocating the memory.
 */


#include <Python.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numpy/arrayobject.h>

// An object responsible for deallocating the memory
typedef struct {
  PyObject_HEAD
  void *memory;
  int typenum; // NPY_DOUBLE or NPY_INT
} _DeallocObject;


static void
_dealloc(PyObject *self)
{
  if (((_DeallocObject *)self)->typenum == NPY_DOUBLE) {
    std::vector<double> *v;
    v = (std::vector<double>*) ((_DeallocObject *)self)->memory;
    v->clear();
  }
  if (((_DeallocObject *)self)->typenum == NPY_INT) {
    std::vector<int> *v;
    v = (std::vector<int>*) ((_DeallocObject *)self)->memory;
    v->clear();
  }
  self->ob_type->tp_free(self);
}

static PyTypeObject _DeallocType = {
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "deallocator", /*tp_name*/
  sizeof(_DeallocObject), /*tp_basicsize*/
  0, /*tp_itemsize*/
  _dealloc, /*tp_dealloc*/
  0, /*tp_print*/
  0, /*tp_getattr*/
  0, /*tp_setattr*/
  0, /*tp_compare*/
  0, /*tp_repr*/
  0, /*tp_as_number*/
  0, /*tp_as_sequence*/
  0, /*tp_as_mapping*/
  0, /*tp_hash */
  0, /*tp_call*/
  0, /*tp_str*/
  0, /*tp_getattro*/
  0, /*tp_setattro*/
  0, /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT, /*tp_flags*/
  "Internal deallocator object", /* tp_doc */
};

// convert a C++ vector to a 1d-ndarray WITHOUT memory copying
template <class T>
static PyObject*
_to_1d_array(std::vector<T> *data, int typenum)
{
  npy_intp dims[1] = {data->size()};

  PyObject *arr = NULL;

  // A C++ vector's first element is guaranteed to point to the internally used
  // array of memory (memory is contiguous).
  arr = PyArray_SimpleNewFromData(1, dims, typenum, (void *) &(*data)[0]);

  if (arr == NULL)
    goto fail;

  _DeallocObject *newobj;
  newobj = PyObject_New(_DeallocObject, &_DeallocType);

  if (newobj == NULL)
    goto fail;

  ((_DeallocObject *)newobj)->memory = (void *)data;
  ((_DeallocObject *)newobj)->typenum = typenum;

  PyArray_BASE(arr) = (PyObject *)newobj;

  return arr;

fail:
  data->clear();
  Py_XDECREF(arr);

  return NULL;
}

static bool
_parse_line(const std::string& line,
            std::vector<double> *data,
            std::vector<int> *indices,
            std::vector<int> *indptr,
            std::vector<double> *labels)
{
  const char *in_string = line.c_str();

  int length = strlen(in_string);
  if (length == 0)
    return false;

  // Parse label
  double y;

  if (!sscanf(in_string, "%lf", &y)) {
    return false;
  }

  labels->push_back(y);

  const char* position;
  position = strchr(in_string, ' ') + 1;

  indptr->push_back(data->size());

  // Parse feature-value pairs
  for ( ;
       (position < in_string + length
      && position - 1 != NULL
      && position[0] != '#');
       position = strchr(position, ' ') + 1) {

    // Consume multiple spaces, if needed.
    if (position[0] == ' ' || position[0] == '\n' ||
        position[0] == '\v' || position[0] == '\r') {
      continue;
    };

    // Parse the feature-value pair.
    int id = atoi(position);
    position = strchr(position, ':') + 1;
    double value = atof(position);
    indices->push_back(id);
    data->push_back(value);
  }

  return true;
}

static bool
_parse_file(char *file_path,
            long int buffer_size,
            std::vector<double> *data,
            std::vector<int> *indices,
            std::vector<int> *indptr,
            std::vector<double> *labels)
{
  char* buffer = new char[buffer_size];

  std::ifstream file_stream(file_path, std::ifstream::in);
  file_stream.rdbuf()->pubsetbuf(buffer, buffer_size);

  if (file_stream) {
    std::string line;
    while (getline(file_stream, line)) {
      bool succeeded = _parse_line(line, data, indices, indptr, labels);
      if (!succeeded)
        return false;
    }
    indptr->push_back(data->size());
  }

  delete[] buffer;

  return true;
}

static char _load_svmlight_format_doc[] =
"Load file in svmlight format and return a CSR.";

static PyObject*
_load_svmlight_format(PyObject *self, PyObject *args)
{

  // initialization
  _import_array();

  _DeallocType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&_DeallocType) < 0)
    return NULL;

  std::vector<double> *data = new std::vector<double>;
  std::vector<int> *indices = new std::vector<int>;
  std::vector<int> *indptr = new std::vector<int>; // of size n_samples + 1
  std::vector<double> *labels = new std::vector<double>;

  // read function arguments
  char *file_path;
  int buffer_mb;

  if (!PyArg_ParseTuple(args, "si", &file_path, &buffer_mb)) {
    return NULL;
  }

  long int buffer_size = buffer_mb * 1024 * 1024;

  // parse file
  bool succeeded = _parse_file(file_path, buffer_size, data, indices, indptr, labels);

  if (succeeded) {
    // return a tuple of ndarrays
    return Py_BuildValue("OOOO",
                         _to_1d_array(data, NPY_DOUBLE),
                         _to_1d_array(indices, NPY_INT),
                         _to_1d_array(indptr, NPY_INT),
                         _to_1d_array(labels, NPY_DOUBLE));
  }
  else {
    return Py_BuildValue("()");
  }

}

static PyMethodDef _svmlight_format_methods[] = {
  {"_load_svmlight_format", _load_svmlight_format,
    METH_VARARGS, _load_svmlight_format_doc},
  {NULL, NULL, 0, NULL}
};

static char _svmlight_format_doc[] =
"Module _svmlight_format.";

PyMODINIT_FUNC
init_svmlight_format(void)
{

  Py_InitModule3("_svmlight_format",
                 _svmlight_format_methods,
                 _svmlight_format_doc);
}


