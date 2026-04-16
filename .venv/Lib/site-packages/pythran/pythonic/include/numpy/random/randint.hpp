#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDINT_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDINT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS>
    std::enable_if_t<!std::is_integral<pS>::value, types::ndarray<long, pS>>
    randint(long min, long max, pS const &shape);

    template <class pS>
    std::enable_if_t<std::is_integral<pS>::value, types::ndarray<long, types::pshape<long>>>
    randint(long min, long max, pS const &shape);

    template <class pS>
    auto randint(long max, types::none_type, pS const &shape) -> decltype(randint(0, max, shape));

    long randint(long min, long max);

    long randint(long max, types::none_type = {});

    auto randint(long min, long max, long size)
        -> decltype(randint(min, max, types::array_tuple<long, 1>{{size}}));

    DEFINE_FUNCTOR(pythonic::numpy::random, randint);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
