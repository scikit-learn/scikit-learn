#ifndef PYTHONIC_NUMPY_ARRAYEQUAL_HPP
#define PYTHONIC_NUMPY_ARRAYEQUAL_HPP

#include "pythonic/include/numpy/array_equal.hpp"

#include "pythonic/numpy/all.hpp"
#include "pythonic/numpy/equal.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class U, class V>
  bool array_equal(U const &u, V const &v)
  {
    if (sutils::getshape(u) == sutils::getshape(v))
      return all(functor::equal{}(u, v));
    return false;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
