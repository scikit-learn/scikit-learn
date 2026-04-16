#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_GEOMETRIC_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_GEOMETRIC_HPP

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
    types::ndarray<double, pS> geometric(double p, pS const &shape);

    auto geometric(double p, long size)
        -> decltype(geometric(p, types::array_tuple<long, 1>{{size}}));

    double geometric(double, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, geometric);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
