#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_CHISQUARE_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_CHISQUARE_HPP

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
    types::ndarray<double, pS> chisquare(double df, pS const &shape);

    auto chisquare(double df, long size)
        -> decltype(chisquare(df, types::array_tuple<long, 1>{{size}}));

    double chisquare(double df, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, chisquare);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
