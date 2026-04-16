#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_DIRICHLET_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_DIRICHLET_HPP

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
    types::ndarray<double, pS> dirichlet(double alpha, pS const &shape);

    auto dirichlet(double alpha, long size)
        -> decltype(dirichlet(alpha, types::array_tuple<long, 1>{{size}}));

    double dirichlet(double alpha, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, dirichlet);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
