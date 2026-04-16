#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_BINOMIAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_BINOMIAL_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/types/numpy_expr.hpp"

#include "pythonic/include/numpy/random/generator.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS>
    types::ndarray<long, pS> binomial(double n, double p, pS const &shape);

    auto binomial(double n, double p, long size)
        -> decltype(binomial(n, p, types::array_tuple<long, 1>{{size}}));

    long binomial(double n, double p, types::none_type d = types::none_type());

    DEFINE_FUNCTOR(pythonic::numpy::random, binomial);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
