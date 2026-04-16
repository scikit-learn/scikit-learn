#ifndef PYTHONIC_NUMPY_FULLLIKE_HPP
#define PYTHONIC_NUMPY_FULLLIKE_HPP

#include "pythonic/include/numpy/full_like.hpp"

#include "pythonic/numpy/full.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class F, class dtype>
  auto full_like(E const &expr, F fill_value, dtype d)
      -> decltype(full(sutils::getshape(expr), fill_value, d))
  {
    return full(sutils::getshape(expr), fill_value, d);
  }

  template <class E, class F>
  auto full_like(E const &expr, F fill_value, types::none_type)
      -> decltype(full(sutils::getshape(expr), fill_value, types::dtype_t<typename E::dtype>()))
  {
    return full(sutils::getshape(expr), fill_value, types::dtype_t<typename E::dtype>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
