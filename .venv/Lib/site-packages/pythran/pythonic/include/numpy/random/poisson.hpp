#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_POISSON_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_POISSON_HPP

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
    types::ndarray<double, pS> poisson(double lam, pS const &shape);

    auto poisson(double lam, long size)
        -> decltype(poisson(lam, types::array_tuple<long, 1>{{size}}));

    double poisson(double lam = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, poisson);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
