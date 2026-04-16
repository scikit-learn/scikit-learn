#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_F_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_F_HPP

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
    types::ndarray<double, pS> f(double dfnum, double dfden, pS const &shape);

    auto f(double dfnum, double dfden, long size)
        -> decltype(f(dfnum, dfden, types::array_tuple<long, 1>{{size}}));

    double f(double dfnum, double dfden, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, f);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
