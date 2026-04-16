//
//  From NumpySrc/numpy/core/src/multiarray/compiled_base.c

/** @brief find index of a sorted array such that arr[i] <= key < arr[i + 1].
 *
 * If an starting index guess is in-range, the array values around this
 * index are first checked.  This allows for repeated calls for well-ordered
 * keys (a very common case) to use the previous index as a very good guess.
 *
 * If the guess value is not useful, bisection of the array is used to
 * find the index.  If there is no such index, the return values are:
 *     key < arr[0] -- -1
 *     key == arr[len - 1] -- len - 1
 *     key > arr[len - 1] -- len
 * The array is assumed contiguous and sorted in ascending order.
 *
 * @param key key value.
 * @param arr contiguous sorted array to be searched.
 * @param len length of the array.
 * @param guess initial guess of index
 * @return index
 */
#include "pythonic/numpy/isnan.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"

#define LIKELY_IN_CACHE_SIZE 8
PYTHONIC_NS_BEGIN

template <typename npy_intp, typename npy_double, class T>
static npy_intp binary_search_with_guess(const npy_double key, const T &arr, npy_intp len,
                                         npy_intp guess)
{
  npy_intp imin = 0;
  npy_intp imax = len;

  /* Handle keys outside of the arr range first */
  if (key > arr[len - 1]) {
    return len;
  } else if (key < arr[0]) {
    return -1;
  }

  /*
   * If len <= 4 use linear search.
   * From above we know key >= arr[0] when we start.
   */
  if (len <= 4) {
    npy_intp i;

    for (i = 1; i < len && key >= arr[i]; ++i)
      ;
    return i - 1;
  }

  if (guess > len - 3) {
    guess = len - 3;
  }
  if (guess < 1) {
    guess = 1;
  }

  /* check most likely values: guess - 1, guess, guess + 1 */
  if (key < arr[guess]) {
    if (key < arr[guess - 1]) {
      imax = guess - 1;
      /* last attempt to restrict search to items in cache */
      if (guess > LIKELY_IN_CACHE_SIZE && key >= arr[guess - LIKELY_IN_CACHE_SIZE]) {
        imin = guess - LIKELY_IN_CACHE_SIZE;
      }
    } else {
      /* key >= arr[guess - 1] */
      return guess - 1;
    }
  } else {
    /* key >= arr[guess] */
    if (key < arr[guess + 1]) {
      return guess;
    } else {
      /* key >= arr[guess + 1] */
      if (key < arr[guess + 2]) {
        return guess + 1;
      } else {
        /* key >= arr[guess + 2] */
        imin = guess + 2;
        /* last attempt to restrict search to items in cache */
        if (guess < len - LIKELY_IN_CACHE_SIZE - 1 && key < arr[guess + LIKELY_IN_CACHE_SIZE]) {
          imax = guess + LIKELY_IN_CACHE_SIZE;
        }
      }
    }
  }

  /* finally, find index by bisection */
  while (imin < imax) {
    const npy_intp imid = imin + ((imax - imin) >> 1);
    if (key >= arr[imid]) {
      imin = imid + 1;
    } else {
      imax = imid;
    }
  }
  return imin - 1;
}
//
// #undef LIKELY_IN_CACHE_SIZE
//
// NPY_NO_EXPORT PyObject *
// arr_interp(PyObject *NPY_UNUSED(self), PyObject *args, PyObject *kwdict)
//{
//
//    PyObject *fp, *xp, *x;
//    PyObject *left = NULL, *right = NULL;
//    PyArrayObject *afp = NULL, *axp = NULL, *ax = NULL, *af = NULL;
//    npy_intp i, lenx, lenxp;
//    double lval, rval;
//    const double *dy, *dx, *dz;
//    double *dres, *slopes = NULL;
//
//    static char *kwlist[] = {"x", "xp", "fp", "left", "right", NULL};
//
//    NPY_BEGIN_THREADS_DEF;
//
//    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOO|OO:interp", kwlist,
//                                     &x, &xp, &fp, &left, &right)) {
//        return NULL;
//    }
//
//    afp = (PyArrayObject *)PyArray_ContiguousFromAny(fp, double, 1, 1);
//    if (afp == NULL) {
//        return NULL;
//    }
//    axp = (PyArrayObject *)PyArray_ContiguousFromAny(xp, double, 1, 1);
//    if (axp == NULL) {
//        goto fail;
//    }
//    ax = (PyArrayObject *)PyArray_ContiguousFromAny(x, double, 0, 0);
//    if (ax == NULL) {
//        goto fail;
//    }
//    lenxp = PyArray_SIZE(axp);
//    if (lenxp == 0) {
//        PyErr_SetString(PyExc_ValueError,
//                        "array of sample points is empty");
//        goto fail;
//    }
//    if (PyArray_SIZE(afp) != lenxp) {
//        PyErr_SetString(PyExc_ValueError,
//                        "fp and xp are not of the same length.");
//        goto fail;
//    }
//
//    af = (PyArrayObject *)PyArray_SimpleNew(PyArray_NDIM(ax),
//                                            PyArray_DIMS(ax), double);
//    if (af == NULL) {
//        goto fail;
//    }
//    lenx = PyArray_SIZE(ax);
//
//    dy = (const double *)PyArray_DATA(afp);
//    dx = (const double *)PyArray_DATA(axp);
//    dz = (const double *)PyArray_DATA(ax);
//    dres = (double *)PyArray_DATA(af);
//    /* Get left and right fill values. */
//    if ((left == NULL) || (left == Py_None)) {
//        lval = dy[0];
//    }
//    else {
//        lval = PyFloat_AsDouble(left);
//        if (error_converting(lval)) {
//            goto fail;
//        }
//    }
//    if ((right == NULL) || (right == Py_None)) {
//        rval = dy[lenxp - 1];
//    }
//    else {
//        rval = PyFloat_AsDouble(right);
//        if (error_converting(rval)) {
//            goto fail;
//        }
//    }

// xp->dx   fp->dy  x -> dz
// This is the output type, based on the type of T, which can be complex.
template <class T>
using out_type =
    std::conditional_t<types::is_complex<typename T::dtype>::value, std::complex<double>, double>;

template <typename npy_intp, typename T5, class T1, class T2, class T3, class T4>
void do_interp(const T1 &dz, const T2 &dx, const T3 &dy, T4 &dres, npy_intp lenxp, npy_intp lenx,
               T5 lval, T5 rval)
{
  npy_intp i;
  out_type<T3> *slopes = NULL;
  std::vector<out_type<T3>, utils::allocator<out_type<T3>>> slope_vect;
  /* binary_search_with_guess needs at least a 3 item long array */
  if (lenxp == 1) {
    const double xp_val = dx[0];
    const out_type<T3> fp_val = dy[0];

    //        NPY_BEGIN_THREADS_THRESHOLDED(lenx);
    for (i = 0; i < lenx; ++i) {
      const double x_val = dz[i];
      dres[i] = (x_val < xp_val) ? lval : ((x_val > xp_val) ? rval : fp_val);
    }
    //        NPY_END_THREADS;
  } else {
    npy_intp j = 0;

    /* only pre-calculate slopes if there are relatively few of them. */
    if (lenxp <= lenx) {
      slope_vect.resize(lenxp - 1);
      slopes = slope_vect.data();
    }

    //        NPY_BEGIN_THREADS;

    if (slopes != NULL) {
      for (i = 0; i < lenxp - 1; ++i) {
        slopes[i] = (dy[i + 1] - dy[i]) / (dx[i + 1] - dx[i]);
      }
    }

    for (i = 0; i < lenx; ++i) {
      const double x_val = dz[i];

      if (pythonic::numpy::functor::isnan()(x_val)) {
        dres[i] = x_val;
        continue;
      }

      j = binary_search_with_guess(x_val, dx, lenxp, j);
      if (j == -1) {
        dres[i] = lval;
      } else if (j == lenxp) {
        dres[i] = rval;
      } else if (j == lenxp - 1) {
        dres[i] = dy[j];
      } else if (dx[j] == x_val) {
        /* Avoid potential non-finite interpolation */
        dres[i] = dy[j];
      } else {
        const out_type<T3> slope =
            (slopes != NULL) ? slopes[j] : (dy[j + 1] - dy[j]) / (dx[j + 1] - dx[j]);
        dres[i] = slope * (x_val - dx[j]) + dy[j];
      }
    }

    //        NPY_END_THREADS;
  }
}
PYTHONIC_NS_END
