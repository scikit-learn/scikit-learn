#ifndef PYTHONIC_NUMPY_ONESLIKE_HPP
#define PYTHONIC_NUMPY_ONESLIKE_HPP

#include "pythonic/include/numpy/ones_like.hpp"

#include "pythonic/numpy/ones.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  auto ones_like(E const &expr, dtype d) -> decltype(ones(sutils::getshape(expr), d))
  {
    return ones(sutils::getshape(expr), d);
  }

  template <class E>
  auto ones_like(E const &expr, types::none_type)
      -> decltype(ones(sutils::getshape(expr), types::dtype_t<typename E::dtype>()))
  {
    return ones(sutils::getshape(expr), types::dtype_t<typename E::dtype>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
