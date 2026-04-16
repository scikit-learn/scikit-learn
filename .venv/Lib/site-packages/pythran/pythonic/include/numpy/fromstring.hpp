#ifndef PYTHONIC_INCLUDE_NUMPY_FROMSTRING_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMSTRING_HPP

#include "pythonic/include/builtins/pythran/kwonly.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <limits>
#include <sstream>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromstring(types::str const &string, dtype d = dtype(), long count = -1, types::kwonly = {},
             types::str const &sep = {});

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromstring(types::str const &string, dtype d = dtype(), long count = -1,
             types::str const &sep = {})
  {
    return fromstring(string, d, count, types::kwonly{}, sep);
  }

  types::ndarray<typename functor::float64::type, types::pshape<long>> inline fromstring(
      types::str const &string, types::kwonly, types::str const &sep = {})
  {
    return fromstring(string, functor::float64{}, -1, types::kwonly{}, sep);
  }

  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromstring(types::str const &string, dtype d, types::kwonly, types::str const &sep)
  {
    return fromstring(string, d, -1, types::kwonly{}, sep);
  }

  DEFINE_FUNCTOR(pythonic::numpy, fromstring);
} // namespace numpy
PYTHONIC_NS_END

#endif
