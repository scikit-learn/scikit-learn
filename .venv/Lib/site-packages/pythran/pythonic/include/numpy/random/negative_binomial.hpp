#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_NEGATIVE_BINOMIAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_NEGATIVE_BINOMIAL_HPP

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
    types::ndarray<long, pS> negative_binomial(long n, double p, pS const &shape);

    auto negative_binomial(long n, double p, long size)
        -> decltype(negative_binomial(n, p, types::array_tuple<long, 1>{{size}}));

    long negative_binomial(long n, double p, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, negative_binomial);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
