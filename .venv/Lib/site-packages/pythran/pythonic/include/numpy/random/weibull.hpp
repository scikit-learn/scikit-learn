#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_WEIBULL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_WEIBULL_HPP

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
    types::ndarray<double, pS> weibull(double a, pS const &shape);

    auto weibull(double a, long size) -> decltype(weibull(a, types::array_tuple<long, 1>{{size}}));

    double weibull(double a, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, weibull);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
