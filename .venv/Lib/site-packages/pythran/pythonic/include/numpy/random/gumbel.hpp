#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_GUMBEL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_GUMBEL_HPP

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
    types::ndarray<double, pS> gumbel(double loc, double scale, pS const &shape);

    auto gumbel(double loc, double scale, long size)
        -> decltype(gumbel(loc, scale, types::array_tuple<long, 1>{{size}}));

    double gumbel(double loc = 0.0, double scale = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, gumbel);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
