#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_EXPONENTIAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_EXPONENTIAL_HPP

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
    types::ndarray<double, pS> exponential(double scale, pS const &shape);

    auto exponential(double scale, long size)
        -> decltype(exponential(scale, types::array_tuple<long, 1>{{size}}));

    double exponential(double scale = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, exponential);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
