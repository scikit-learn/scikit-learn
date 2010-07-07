/* -*- C -*-  (not really, but good for syntax highlighting) */
#ifdef SWIGPYTHON

%{
#ifndef SWIG_FILE_WITH_INIT
#  define NO_IMPORT_ARRAY
#endif
#include "stdio.h"
#include <numpy/arrayobject.h>

/* The following code originally appeared in
 * enthought/kiva/agg/src/numeric.i, author unknown.  It was
 * translated from C++ to C by John Hunter.  Bill Spotz has modified
 * it slightly to fix some minor bugs, upgrade to numpy (all
 * versions), add some comments and some functionality.
 */

/* Macros to extract array attributes.
 */
#define is_array(a)            ((a) && PyArray_Check((PyArrayObject *)a))
#define array_type(a)          (int)(PyArray_TYPE(a))
#define array_dimensions(a)    (((PyArrayObject *)a)->nd)
#define array_size(a,i)        (((PyArrayObject *)a)->dimensions[i])
#define array_is_contiguous(a) (PyArray_ISCONTIGUOUS(a))

/* Support older NumPy data type names
*/
#if NDARRAY_VERSION < 0x01000000
#define NPY_BOOL        PyArray_BOOL
#define NPY_BYTE        PyArray_BYTE
#define NPY_UBYTE       PyArray_UBYTE
#define NPY_SHORT       PyArray_SHORT
#define NPY_USHORT      PyArray_USHORT
#define NPY_INT         PyArray_INT
#define NPY_UINT        PyArray_UINT
#define NPY_LONG        PyArray_LONG
#define NPY_ULONG       PyArray_ULONG
#define NPY_LONGLONG    PyArray_LONGLONG
#define NPY_ULONGLONG   PyArray_ULONGLONG
#define NPY_FLOAT       PyArray_FLOAT
#define NPY_DOUBLE      PyArray_DOUBLE
#define NPY_LONGDOUBLE  PyArray_LONGDOUBLE
#define NPY_CFLOAT      PyArray_CFLOAT
#define NPY_CDOUBLE     PyArray_CDOUBLE
#define NPY_CLONGDOUBLE PyArray_CLONGDOUBLE
#define NPY_OBJECT      PyArray_OBJECT
#define NPY_STRING      PyArray_STRING
#define NPY_UNICODE     PyArray_UNICODE
#define NPY_VOID        PyArray_VOID
#define NPY_NTYPES      PyArray_NTYPES
#define NPY_NOTYPE      PyArray_NOTYPE
#define NPY_CHAR        PyArray_CHAR
#define NPY_USERDEF     PyArray_USERDEF
#define npy_intp        intp
#endif

/* Given a PyObject, return a string describing its type.
 */
char* pytype_string(PyObject* py_obj) {
  if (py_obj == NULL          ) return "C NULL value";
  if (PyCallable_Check(py_obj)) return "callable"    ;
  if (PyString_Check(  py_obj)) return "string"      ;
  if (PyInt_Check(     py_obj)) return "int"         ;
  if (PyFloat_Check(   py_obj)) return "float"       ;
  if (PyDict_Check(    py_obj)) return "dict"        ;
  if (PyList_Check(    py_obj)) return "list"        ;
  if (PyTuple_Check(   py_obj)) return "tuple"       ;
  if (PyFile_Check(    py_obj)) return "file"        ;
  if (PyModule_Check(  py_obj)) return "module"      ;
  if (PyInstance_Check(py_obj)) return "instance"    ;

  return "unkown type";
}

/* Given a NumPy typecode, return a string describing the type.
 */
char* typecode_string(int typecode) {
  static char* type_names[25] = {"bool", "byte", "unsigned byte",
				 "short", "unsigned short", "int",
				 "unsigned int", "long", "unsigned long",
				 "long long", "unsigned long long",
				 "float", "double", "long double",
				 "complex float", "complex double",
				 "complex long double", "object",
				 "string", "unicode", "void", "ntypes",
				 "notype", "char", "unknown"};
  return typecode < 24 ? type_names[typecode] : type_names[24];
}

/* Make sure input has correct numpy type.  Allow character and byte
 * to match.  Also allow int and long to match.  This is deprecated.
 * You should use PyArray_EquivTypenums() instead.
 */
int type_match(int actual_type, int desired_type) {
  return PyArray_EquivTypenums(actual_type, desired_type);
}

/* Given a PyObject pointer, cast it to a PyArrayObject pointer if
 * legal.  If not, set the python error string appropriately and
 * return NULL.
 */
PyArrayObject* obj_to_array_no_conversion(PyObject* input, int typecode) {
  PyArrayObject* ary = NULL;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input), typecode))) {
    ary = (PyArrayObject*) input;
  }
  else if is_array(input) {
    char* desired_type = typecode_string(typecode);
    char* actual_type  = typecode_string(array_type(input));
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  Array of type '%s' given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  else {
    char * desired_type = typecode_string(typecode);
    char * actual_type  = pytype_string(input);
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  A '%s' was given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  return ary;
}

/* Convert the given PyObject to a NumPy array with the given
 * typecode.  On success, return a valid PyArrayObject* with the
 * correct type.  On failure, the python error string will be set and
 * the routine returns NULL.
 */
PyArrayObject* obj_to_array_allow_conversion(PyObject* input, int typecode,
                                             int* is_new_object) {
  PyArrayObject* ary = NULL;
  PyObject* py_obj;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input),typecode))) {
    ary = (PyArrayObject*) input;
    *is_new_object = 0;
  }
  else {
    py_obj = PyArray_FromObject(input, typecode, 0, 0);
    /* If NULL, PyArray_FromObject will have set python error value.*/
    ary = (PyArrayObject*) py_obj;
    *is_new_object = 1;
  }
  return ary;
}

/* Given a PyArrayObject, check to see if it is contiguous.  If so,
 * return the input pointer and flag it as not a new object.  If it is
 * not contiguous, create a new PyArrayObject using the original data,
 * flag it as a new object and return the pointer.
 */
PyArrayObject* make_contiguous(PyArrayObject* ary, int* is_new_object,
                               int min_dims, int max_dims) {
  PyArrayObject* result;
  if (array_is_contiguous(ary)) {
    result = ary;
    *is_new_object = 0;
  }
  else {
    result = (PyArrayObject*) PyArray_ContiguousFromObject((PyObject*)ary, 
							   array_type(ary), 
							   min_dims,
							   max_dims);
    *is_new_object = 1;
  }
  return result;
}

/* Convert a given PyObject to a contiguous PyArrayObject of the
 * specified type.  If the input object is not a contiguous
 * PyArrayObject, a new one will be created and the new object flag
 * will be set.
 */
PyArrayObject* obj_to_array_contiguous_allow_conversion(PyObject* input,
                                                        int typecode,
                                                        int* is_new_object) {
  int is_new1 = 0;
  int is_new2 = 0;
  PyArrayObject* ary2;
  PyArrayObject* ary1 = obj_to_array_allow_conversion(input, typecode, 
						      &is_new1);
  if (ary1) {
    ary2 = make_contiguous(ary1, &is_new2, 0, 0);
    if ( is_new1 && is_new2) {
      Py_DECREF(ary1);
    }
    ary1 = ary2;    
  }
  *is_new_object = is_new1 || is_new2;
  return ary1;
}

/* Test whether a python object is contiguous.  If array is
 * contiguous, return 1.  Otherwise, set the python error string and
 * return 0.
 */
int require_contiguous(PyArrayObject* ary) {
  int contiguous = 1;
  if (!array_is_contiguous(ary)) {
    PyErr_SetString(PyExc_TypeError,
		    "Array must be contiguous.  A non-contiguous array was given");
    contiguous = 0;
  }
  return contiguous;
}

/* Require the given PyArrayObject to have a specified number of
 * dimensions.  If the array has the specified number of dimensions,
 * return 1.  Otherwise, set the python error string and return 0.
 */
int require_dimensions(PyArrayObject* ary, int exact_dimensions) {
  int success = 1;
  if (array_dimensions(ary) != exact_dimensions) {
    PyErr_Format(PyExc_TypeError, 
		 "Array must have %d dimensions.  Given array has %d dimensions", 
		 exact_dimensions, array_dimensions(ary));
    success = 0;
  }
  return success;
}

/* Require the given PyArrayObject to have one of a list of specified
 * number of dimensions.  If the array has one of the specified number
 * of dimensions, return 1.  Otherwise, set the python error string
 * and return 0.
 */
int require_dimensions_n(PyArrayObject* ary, int* exact_dimensions, int n) {
  int success = 0;
  int i;
  char dims_str[255] = "";
  char s[255];
  for (i = 0; i < n && !success; i++) {
    if (array_dimensions(ary) == exact_dimensions[i]) {
      success = 1;
    }
  }
  if (!success) {
    for (i = 0; i < n-1; i++) {
      sprintf(s, "%d, ", exact_dimensions[i]);                
      strcat(dims_str,s);
    }
    sprintf(s, " or %d", exact_dimensions[n-1]);            
    strcat(dims_str,s);
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have %s dimensions.  Given array has %d dimensions",
		 dims_str, array_dimensions(ary));
  }
  return success;
}    

/* Require the given PyArrayObject to have a specified shape.  If the
 * array has the specified shape, return 1.  Otherwise, set the python
 * error string and return 0.
 */
int require_size(PyArrayObject* ary, npy_intp* size, int n) {
  int i;
  int success = 1;
  int len;
  char desired_dims[255] = "[";
  char s[255];
  char actual_dims[255] = "[";
  for(i=0; i < n;i++) {
    if (size[i] != -1 &&  size[i] != array_size(ary,i)) {
      success = 0;    
    }
  }
  if (!success) {
    for (i = 0; i < n; i++) {
      if (size[i] == -1) {
	sprintf(s, "*,");                
      }
      else
      {
	sprintf(s, "%d,", size[i]);                
      }    
      strcat(desired_dims,s);
    }
    len = strlen(desired_dims);
    desired_dims[len-1] = ']';
    for (i = 0; i < n; i++) {
      sprintf(s, "%d,", array_size(ary,i));                            
      strcat(actual_dims,s);
    }
    len = strlen(actual_dims);
    actual_dims[len-1] = ']';
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have shape of %s.  Given array has shape of %s",
		 desired_dims, actual_dims);
  }
  return success;
}

/* End John Hunter translation (with modifications by Bill Spotz)
 */

%}

/* %numpy_typemaps() macro
 *
 * This macro defines a family of typemaps that allow pure input C
 * arguments of the form
 *
 *     (DATA_TYPE IN_ARRAY1[ANY])
 *     (DATA_TYPE* IN_ARRAY1, DIM_TYPE DIM1)
 *     (DIM_TYPE DIM1, DATA_TYPE* IN_ARRAY1)
 *
 *     (DATA_TYPE IN_ARRAY2[ANY][ANY])
 *     (DATA_TYPE* IN_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 *     (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_ARRAY2)
 *
 *     (DATA_TYPE INPLACE_ARRAY1[ANY])
 *     (DATA_TYPE* INPLACE_ARRAY1, DIM_TYPE DIM1)
 *     (DIM_TYPE DIM1, DATA_TYPE* INPLACE_ARRAY1)
 *
 *     (DATA_TYPE INPLACE_ARRAY2[ANY][ANY])
 *     (DATA_TYPE* INPLACE_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 *     (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
 *
 *     (DATA_TYPE ARGOUT_ARRAY1[ANY])
 *     (DATA_TYPE* ARGOUT_ARRAY1, DIM_TYPE DIM1)
 *     (DIM_TYPE DIM1, DATA_TYPE* ARGOUT_ARRAY1)
 *
 *     (DATA_TYPE ARGOUT_ARRAY2[ANY][ANY])
 *
 * where "DATA_TYPE" is any type supported by the NumPy module, and
 * "DIM_TYPE" is any int-like type suitable for specifying dimensions.
 * In python, the dimensions will not need to be specified (except for
 * the "DATA_TYPE* ARGOUT_ARRAY1" typemaps).  The IN_ARRAYs can be a
 * numpy array or any sequence that can be converted to a numpy array
 * of the specified type.  The INPLACE_ARRAYs must be numpy arrays of
 * the appropriate type.  The ARGOUT_ARRAYs will be returned as numpy
 * arrays of the appropriate type.
 *
 * These typemaps can be applied to existing functions using the
 * %apply directive:
 *
 *     %apply (double IN_ARRAY1[ANY]) {(double vector[ANY])};
 *     double length(double vector[3]);
 *
 *     %apply (double* IN_ARRAY1, int DIM1) {(double* series, int length)};
 *     double prod(double* series, int length);
 *
 *     %apply (int DIM1, double* IN_ARRAY1) {(int length, double* series)}
 *     double sum(int length, double* series)
 *
 *     %apply (double IN_ARRAY2[ANY][ANY]) {(double matrix[2][2])};
 *     double det(double matrix[2][2]);
 *
 *     %apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};
 *     double max(double* matrix, int rows, int cols);
 *
 *     %apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int rows, int cols, double* matrix)}
 *     double min(int length, double* series)
 *
 *     %apply (double INPLACE_ARRAY1[ANY]) {(double vector[3])};
 *     void reverse(double vector[3]);
 *
 *     %apply (double* INPLACE_ARRAY1, int DIM1) {(double* series, int length)};
 *     void ones(double* series, int length);
 *
 *     %apply (int DIM1, double* INPLACE_ARRAY1) {(int length, double* series)}
 *     double zeros(int length, double* series)
 *
 *     %apply (double INPLACE_ARRAY2[ANY][ANY]) {(double matrix[3][3])};
 *     void scale(double matrix[3][3]);
 *
 *     %apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};
 *     void floor(double* matrix, int rows, int cols);
 *
 *     %apply (int DIM1, int DIM2, double* INPLACE_ARRAY2) {(int rows, int cols, double* matrix)};
 *     void ceil(int rows, int cols, double* matrix);
 *
 *     %apply (double IN_ARRAY1[ANY]    ) {(double vector[ANY])};
 *     %apply (double ARGOUT_ARRAY1[ANY]) {(double even[    3])};
 *     %apply (double ARGOUT_ARRAY1[ANY]) {(double odd[     3])};
 *     void eoSplit(double vector[3], double even[3], double odd[3]);
 *
 *     %apply (double* ARGOUT_ARRAY1, int DIM1) {(double* twoVec, int size)};
 *     void twos(double* twoVec, int size);
 *
 *     %apply (int DIM1, double* ARGOUT_ARRAY1) {(int size, double* threeVec)};
 *     void threes(int size, double* threeVec);
 *
 *     %apply (double IN_ARRAY2[ANY][ANY])     {(double matrix[2][2])};
 *     %apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double upper[ 3][3])};
 *     %apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double lower[ 3][3])};
 *     void luSplit(double matrix[3][3], double upper[3][3], double lower[3][3]);
 *
 * or directly with
 *
 *     double length(double IN_ARRAY1[ANY]);
 *     double prod(double* IN_ARRAY1, int DIM1);
 *     double sum( int DIM1, double* IN_ARRAY1)
 *
 *     double det(double IN_ARRAY2[ANY][ANY]);
 *     double max(double* IN_ARRAY2, int DIM1, int DIM2);
 *     double min(int DIM1, int DIM2, double* IN_ARRAY2)
 *
 *     void reverse(double INPLACE_ARRAY1[ANY]);
 *     void ones( double* INPLACE_ARRAY1, int DIM1);
 *     void zeros(int DIM1, double* INPLACE_ARRAY1)
 *
 *     void scale(double INPLACE_ARRAY2[ANY][ANY]);
 *     void floor(double* INPLACE_ARRAY2, int DIM1, int DIM2, double floor);
 *     void ceil( int DIM1, int DIM2, double* INPLACE_ARRAY2, double ceil );
 *
 *     void eoSplit(double IN_ARRAY1[ANY], double ARGOUT_ARRAY1[ANY],
 *                  double ARGOUT_ARRAY1[ANY]);
 *     void twos(double* ARGOUT_ARRAY1, int DIM1)
 *     void threes(int DIM1, double* ARGOUT_ARRAY1)
 *
 *     void luSplit(double IN_ARRAY2[ANY][ANY], double ARGOUT_ARRAY2[ANY][ANY],
 *                  double ARGOUT_ARRAY2[ANY][ANY]);
 */

%define %numpy_typemaps(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

/* Typemap suite for (DATA_TYPE IN_ARRAY1[ANY])
 */
%typemap(in)
  (DATA_TYPE IN_ARRAY1[ANY])
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[1] = { $1_dim0 };
  if (!array || !require_dimensions(array, 1) || !require_size(array, size, 1)) SWIG_fail;
  $1 = ($1_ltype) array->data;
}
%typemap(freearg)
  (DATA_TYPE IN_ARRAY1[ANY])
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DATA_TYPE* IN_ARRAY1, DIM_TYPE DIM1)
 */
%typemap(in)
  (DATA_TYPE* IN_ARRAY1, DIM_TYPE DIM1)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[1] = { -1 };
  if (!array || !require_dimensions(array, 1) || !require_size(array, size, 1)) SWIG_fail;
  $1 = (DATA_TYPE*) array->data;
  $2 = (DIM_TYPE) array->dimensions[0];
}
%typemap(freearg)
  (DATA_TYPE* IN_ARRAY1, DIM_TYPE DIM1)
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DATA_TYPE* IN_ARRAY1)
 */
%typemap(in)
  (DIM_TYPE DIM1, DATA_TYPE* IN_ARRAY1)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[1] = {-1};
  if (!array || !require_dimensions(array, 1) || !require_size(array, size, 1)) SWIG_fail;
  $1 = (DIM_TYPE) array->dimensions[0];
  $2 = (DATA_TYPE*) array->data;
}
%typemap(freearg)
  (DIM_TYPE DIM1, DATA_TYPE* IN_ARRAY1)
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DATA_TYPE IN_ARRAY2[ANY][ANY])
 */
%typemap(in)
  (DATA_TYPE IN_ARRAY2[ANY][ANY])
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[2] = { $1_dim0, $1_dim1 };
  if (!array || !require_dimensions(array, 2) || !require_size(array, size, 2)) SWIG_fail;
  $1 = ($1_ltype) array->data;
}
%typemap(freearg)
  (DATA_TYPE IN_ARRAY2[ANY][ANY])
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DATA_TYPE* IN_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 */
%typemap(in)
  (DATA_TYPE* IN_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[2] = {-1,-1};
  if (!array || !require_dimensions(array, 2) || !require_size(array, size, 1)) SWIG_fail;
  $1 = (DATA_TYPE*) array->data;
  $2 = (DIM_TYPE) array->dimensions[0];
  $3 = (DIM_TYPE) array->dimensions[1];
}
%typemap(freearg)
  (DATA_TYPE* IN_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_ARRAY2)
 */
%typemap(in)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_ARRAY2)
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE, &is_new_object);
  npy_intp size[2] = {-1,-1};
  if (!array || !require_dimensions(array, 2) || !require_size(array, size, 1)) SWIG_fail;
  $1 = (DIM_TYPE) array->dimensions[0];
  $2 = (DIM_TYPE) array->dimensions[1];
  $3 = (DATA_TYPE*) array->data;
}
%typemap(freearg)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_ARRAY2)
{
  if (is_new_object$argnum && array$argnum) Py_DECREF(array$argnum);
}

/* Typemap suite for (DATA_TYPE INPLACE_ARRAY1[ANY])
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DATA_TYPE INPLACE_ARRAY1[ANY])
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DATA_TYPE INPLACE_ARRAY1[ANY])
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  npy_intp size[1] = { $1_dim0 };
  if (!temp  || !require_dimensions(temp,1) || !require_size(temp, size, 1)) SWIG_fail;
  $1 = ($1_ltype) temp->data;
}

/* Typemap suite for (DATA_TYPE* INPLACE_ARRAY1, DIM_TYPE DIM1)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DATA_TYPE* INPLACE_ARRAY1, DIM_TYPE DIM1)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DATA_TYPE* INPLACE_ARRAY1, DIM_TYPE DIM1)
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!temp  || !require_contiguous(temp)) SWIG_fail;
  $1 = (DATA_TYPE*) temp->data;
  $2 = 1;
  for (int i=0; i<temp->nd; ++i) $2 *= temp->dimensions[i];
}

/* Typemap suite for (DIM_TYPE DIM1, DATA_TYPE* INPLACE_ARRAY1)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DIM_TYPE DIM1, DATA_TYPE* INPLACE_ARRAY1)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DIM_TYPE DIM1, DATA_TYPE* INPLACE_ARRAY1)
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!temp  || !require_contiguous(temp)) SWIG_fail;
  $1 = 1;
  for (int i=0; i<temp->nd; ++i) $1 *= temp->dimensions[i];
  $2 = (DATA_TYPE*) temp->data;
}

/* Typemap suite for (DATA_TYPE INPLACE_ARRAY2[ANY][ANY])
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DATA_TYPE INPLACE_ARRAY2[ANY][ANY])
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DATA_TYPE INPLACE_ARRAY2[ANY][ANY])
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  npy_intp size[2] = { $1_dim0, $1_dim1 };
  if (!temp  || !require_dimensions(temp,2) || !require_size(temp, size, 2)) SWIG_fail;
  $1 = ($1_ltype) temp->data;
}

/* Typemap suite for (DATA_TYPE* INPLACE_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DATA_TYPE* INPLACE_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DATA_TYPE* INPLACE_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!temp || !require_contiguous(temp)) SWIG_fail;
  $1 = (DATA_TYPE*) temp->data;
  $2 = (DIM_TYPE) temp->dimensions[0];
  $3 = (DIM_TYPE) temp->dimensions[1];
}

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),DATA_TYPECODE);
}
%typemap(in)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
  (PyArrayObject* temp=NULL)
{
  temp = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!temp || !require_contiguous(temp)) SWIG_fail;
  $1 = (DIM_TYPE) temp->dimensions[0];
  $2 = (DIM_TYPE) temp->dimensions[1];
  $3 = (DATA_TYPE*) temp->data;
}

/* Typemap suite for (DATA_TYPE ARGOUT_ARRAY1[ANY])
 */
%typemap(in,numinputs=0)
  (DATA_TYPE ARGOUT_ARRAY1[ANY])
  (PyObject * array = NULL)
{
  npy_intp dims[1] = { $1_dim0 };
  array = PyArray_SimpleNew(1, dims, DATA_TYPECODE);
  $1 = ($1_ltype)((PyArrayObject*)array)->data;
}
%typemap(argout)
  (DATA_TYPE ARGOUT_ARRAY1[ANY])
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DATA_TYPE* ARGOUT_ARRAY1, DIM_TYPE DIM1)
 */
%typemap(in,numinputs=1)
  (DATA_TYPE* ARGOUT_ARRAY1, DIM_TYPE DIM1)
  (PyObject * array = NULL)
{
  if (!PyInt_Check($input)) {
    char* typestring = pytype_string($input);
    PyErr_Format(PyExc_TypeError, 
		 "Int dimension expected.  '%s' given.", 
		 typestring);
    SWIG_fail;
  }
  $2 = (DIM_TYPE) PyInt_AsLong($input);
  npy_intp dims[1] = { (npy_intp) $2 };
  array = PyArray_SimpleNew(1, dims, DATA_TYPECODE);
  $1 = (DATA_TYPE*)((PyArrayObject*)array)->data;
}
%typemap(argout)
  (DATA_TYPE* ARGOUT_ARRAY1, DIM_TYPE DIM1)
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DATA_TYPE* ARGOUT_ARRAY1)
 */
%typemap(in,numinputs=1)
  (DIM_TYPE DIM1, DATA_TYPE* ARGOUT_ARRAY1)
  (PyObject * array = NULL)
{
  if (!PyInt_Check($input)) {
    char* typestring = pytype_string($input);
    PyErr_Format(PyExc_TypeError, 
		 "Int dimension expected.  '%s' given.", 
		 typestring);
    SWIG_fail;
  }
  $1 = (DIM_TYPE) PyInt_AsLong($input);
  npy_intp dims[1] = { (npy_intp) $1 };
  array = PyArray_SimpleNew(1, dims, DATA_TYPECODE);
  $2 = (DATA_TYPE*)((PyArrayObject*)array)->data;
}
%typemap(argout)
  (DIM_TYPE DIM1, DATA_TYPE* ARGOUT_ARRAY1)
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

/* Typemap suite for (DATA_TYPE ARGOUT_ARRAY2[ANY][ANY])
 */
%typemap(in,numinputs=0)
  (DATA_TYPE ARGOUT_ARRAY2[ANY][ANY])
  (PyObject * array = NULL)
{
  npy_intp dims[2] = { $1_dim0, $1_dim1 };
  array = PyArray_SimpleNew(2, dims, DATA_TYPECODE);
  $1 = ($1_ltype)((PyArrayObject*)array)->data;
}
%typemap(argout)
  (DATA_TYPE ARGOUT_ARRAY2[ANY][ANY])
{
  $result = SWIG_Python_AppendOutput($result,array$argnum);
}

%enddef    /* %numpy_typemaps() macro */


/* Concrete instances of the %numpy_typemaps() macro: Each invocation
 * below applies all of the typemaps above to the specified data type.
 */
%numpy_typemaps(signed char       , NPY_BYTE     , int)
%numpy_typemaps(unsigned char     , NPY_UBYTE    , int)
%numpy_typemaps(short             , NPY_SHORT    , int)
%numpy_typemaps(unsigned short    , NPY_USHORT   , int)
%numpy_typemaps(int               , NPY_INT      , int)
%numpy_typemaps(unsigned int      , NPY_UINT     , int)
%numpy_typemaps(long              , NPY_LONG     , int)
%numpy_typemaps(unsigned long     , NPY_ULONG    , int)
%numpy_typemaps(long long         , NPY_LONGLONG , int)
%numpy_typemaps(unsigned long long, NPY_ULONGLONG, int)
%numpy_typemaps(float             , NPY_FLOAT    , int)
%numpy_typemaps(double            , NPY_DOUBLE   , int)

/* ***************************************************************
 * The follow macro expansion does not work, because C++ bool is 4
 * bytes and NPY_BOOL is 1 byte
 */
/*%numpy_typemaps(bool, NPY_BOOL)
 */

/* ***************************************************************
 * On my Mac, I get the following warning for this macro expansion:
 * 'swig/python detected a memory leak of type 'long double *', no destructor found.'
 */
/*%numpy_typemaps(long double, NPY_LONGDOUBLE)
 */

/* ***************************************************************
 * Swig complains about a syntax error for the following macros
 * expansions:
 */
/*%numpy_typemaps(complex float,       NPY_CFLOAT     )
 */
/*%numpy_typemaps(complex double,      NPY_CDOUBLE    )
 */
/*%numpy_typemaps(complex long double, NPY_CLONGDOUBLE)
 */

#endif /* SWIGPYTHON */
