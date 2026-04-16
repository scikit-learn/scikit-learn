#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_GAMMA_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_GAMMA_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS>
    types::ndarray<double, pS> gamma(double shape, double scale, pS const &array_shape);

    auto gamma(double shape, double scale, long size)
        -> decltype(gamma(shape, scale, types::array_tuple<long, 1>{{size}}));

    double gamma(double shape = 0.0, double scale = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, gamma);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
