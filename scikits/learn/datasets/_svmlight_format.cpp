/*
 * Authors: Mathieu Blondel <mathieu@mblondel.org>
 *          Lars Buitinck <L.J.Buitinck@uva.nl>
 *
 * License: Simple BSD
 *
 * This module implements _load_svmlight_format, a fast and memory efficient
 * function to load the file format originally created for svmlight and now used
 * by many other libraries, including libsvm.
 *
 * The function loads the file directly in a CSR sparse matrix without memory
 * copying.  The approach taken is to use 4 C++ vectors (data, indices, indptr
 * and labels) and to incrementally feed them with elements. Ndarrays are then
 * instantiated by PyArray_SimpleNewFromData, i.e., no memory is
 * copied.
 *
 * Since the memory is not allocated by the ndarray, the ndarray doesn't own the
 * memory and thus cannot deallocate it. To automatically deallocate memory, the
 * technique described at http://blog.enthought.com/?p=62 is used. The main idea
 * is to use an additional object that the ndarray does own and that will be
 * responsible for deallocating the memory.
 */


#include <Python.h>
#include <numpy/arrayobject.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>


/*
 * A Python object responsible for memory management of our vectors.
 */
template <typename T>
struct VectorOwner {
  PyObject_HEAD
  std::vector<T> v;
};

/*
 * Deallocators (tp_dealloc). Since a template function can't have C linkage,
 * we define three functions.
 */
template <typename T>
static void destroy_vector_owner(PyObject *self)
{
  // Note: explicit call to destructor because of placement new;
  // memory management for VectorOwner is performed by the interpreter.
  // Compiler-generated dtor will release memory from vector member.
  VectorOwner<T> &obj = *reinterpret_cast<VectorOwner<T> *>(self);
  obj.~VectorOwner();
}

extern "C" {
static void destroy_int_vector(PyObject *self)
{
  destroy_vector_owner<int>(self);
}

static void destroy_double_vector(PyObject *self)
{
  destroy_vector_owner<double>(self);
}
}


/*
 * Type objects for above.
 */
static PyTypeObject IntVOwnerType    = { PyObject_HEAD_INIT(NULL) },
                    DoubleVOwnerType = { PyObject_HEAD_INIT(NULL) };

/*
 * Set the fields of the owner type objects.
 */
static void init_type_objs()
{
  IntVOwnerType.tp_flags = DoubleVOwnerType.tp_flags = Py_TPFLAGS_DEFAULT;
  IntVOwnerType.tp_name  = DoubleVOwnerType.tp_name  = "deallocator";
  IntVOwnerType.tp_doc   = DoubleVOwnerType.tp_doc   = "deallocator object";
  IntVOwnerType.tp_new   = DoubleVOwnerType.tp_new   = PyType_GenericNew;

  IntVOwnerType.tp_basicsize     = sizeof(VectorOwner<int>);
  DoubleVOwnerType.tp_basicsize  = sizeof(VectorOwner<double>);
  IntVOwnerType.tp_dealloc       = destroy_int_vector;
  DoubleVOwnerType.tp_dealloc    = destroy_double_vector;
}

PyTypeObject &vector_owner_type(int typenum)
{
  switch (typenum) {
    case NPY_INT: return IntVOwnerType;
    case NPY_DOUBLE: return DoubleVOwnerType;
  }
  throw std::logic_error("invalid argument to vector_owner_type");
}


/*
 * Convert a C++ vector to a 1d-ndarray WITHOUT memory copying.
 * Steals v's contents and leaves it empty.
 */
template <typename T>
static PyObject *to_1d_array(std::vector<T> &v, int typenum)
{
  npy_intp dims[1] = {v.size()};

  // A C++ vector's contents are guaranteed to be in a contiguous array.
  PyObject *arr = PyArray_SimpleNewFromData(1, dims, typenum, &v[0]);

  try {
    if (!arr)
      throw std::bad_alloc();

    VectorOwner<T> *owner = PyObject_New(VectorOwner<T>,
                                         &vector_owner_type(typenum));
    if (!owner)
      throw std::bad_alloc();

    // transfer ownership of v's contents to the VectorOwner
    new (&owner->v) std::vector<T>();
    owner->v.swap(v);

    PyArray_BASE(arr) = (PyObject *)owner;

    return arr;
  } catch (std::runtime_error const &e) {
    // let's assume the exception is set correctly
    Py_XDECREF(arr);
    return 0;
  }
}


/*
 * Parsing.
 */

class SyntaxError : public std::runtime_error {
public:
  SyntaxError(char const *msg)
   : std::runtime_error(std::string(msg) + " in SVMlight/libSVM file")
  {
  }
};

/*
 * Parse single line. Throws exception on failure.
 */
void parse_line(const std::string& line,
                std::vector<double> &data,
                std::vector<int> &indices,
                std::vector<int> &indptr,
                std::vector<double> &labels)
{
  if (line.length() == 0)
    throw SyntaxError("empty line");

  // Parse label
  // FIXME: this should be done using standard C++ IOstream facilities,
  // so we don't need to read the lines into strings first and get better
  // error handling.
  const char *in_string = line.c_str();
  double y;

  if (!std::sscanf(in_string, "%lf", &y))
    throw SyntaxError("non-numeric or missing label");

  labels.push_back(y);

  const char* position;
  position = std::strchr(in_string, ' ') + 1;

  indptr.push_back(data.size());

  // Parse feature-value pairs
  for ( ;
       (position
      && position < in_string + line.length()
      && position[0] != '#');
       position = std::strchr(position, ' ')) {

    // Consume multiple spaces, if needed.
    while (std::isspace(*position))
      position++;

    // Parse the feature-value pair.
    int id = std::atoi(position);
    position = std::strchr(position, ':') + 1;
    double value = std::atof(position);
    indices.push_back(id);
    data.push_back(value);
  }
}

/*
 * Parse entire file. Throws exception on failure.
 */
void parse_file(char const *file_path,
                size_t buffer_size,
                std::vector<double> &data,
                std::vector<int> &indices,
                std::vector<int> &indptr,
                std::vector<double> &labels)
{
  std::vector<char> buffer(buffer_size);

  std::ifstream file_stream;
  file_stream.exceptions(std::ios::badbit);
  file_stream.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
  file_stream.open(file_path);

  std::string line;
  while (std::getline(file_stream, line))
    parse_line(line, data, indices, indptr, labels);
  indptr.push_back(data.size());
}


static const char load_svmlight_format_doc[] =
  "Load file in svmlight format and return a CSR.";

extern "C" {
static PyObject *load_svmlight_format(PyObject *self, PyObject *args)
{
  // initialization
  if (PyType_Ready(&DoubleVOwnerType) < 0
   || PyType_Ready(&IntVOwnerType)    < 0)
    return 0;

  try {
    std::vector<double> data, labels;
    std::vector<int> indices, indptr;

    // read function arguments
    char const *file_path;
    int buffer_mb;

    if (!PyArg_ParseTuple(args, "si", &file_path, &buffer_mb))
      return 0;

    buffer_mb = std::max(buffer_mb, 1);
    size_t buffer_size = buffer_mb * 1024 * 1024;

    parse_file(file_path, buffer_size, data, indices, indptr, labels);
    return Py_BuildValue("OOOO",
                         to_1d_array(data, NPY_DOUBLE),
                         to_1d_array(indices, NPY_INT),
                         to_1d_array(indptr, NPY_INT),
                         to_1d_array(labels, NPY_DOUBLE));
  } catch (std::exception const &e) {
    // FIXME: we're losing information about the error here.
    return Py_BuildValue("()");
  }
}
}


/*
 * Python module setup.
 */

static PyMethodDef svmlight_format_methods[] = {
  {"_load_svmlight_format", load_svmlight_format,
    METH_VARARGS, load_svmlight_format_doc},
  {NULL, NULL, 0, NULL}
};

static const char svmlight_format_doc[] =
  "Loader for svmlight / libsvm datasets - C++ helper routines";

extern "C" {
PyMODINIT_FUNC init_svmlight_format(void)
{
  _import_array();
  init_type_objs();
  Py_InitModule3("_svmlight_format",
                 svmlight_format_methods,
                 svmlight_format_doc);
}
}
