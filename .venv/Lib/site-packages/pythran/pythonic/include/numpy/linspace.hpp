#ifndef PYTHONIC_INCLUDE_NUMPY_LINSPACE_HPP
#define PYTHONIC_INCLUDE_NUMPY_LINSPACE_HPP

#include "pythonic/include/numpy/arange.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype = types::dtype_t<double>>
  types::ndarray<typename dtype::type, types::pshape<long>>
  linspace(double start, double stop, long num = 50, bool endpoint = true, bool retstep = false,
           dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, linspace);
} // namespace numpy
PYTHONIC_NS_END

#endif
