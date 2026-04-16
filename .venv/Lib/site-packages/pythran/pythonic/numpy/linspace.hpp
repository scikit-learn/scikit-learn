#ifndef PYTHONIC_NUMPY_LINSPACE_HPP
#define PYTHONIC_NUMPY_LINSPACE_HPP

#include "pythonic/include/numpy/linspace.hpp"

#include "pythonic/numpy/arange.hpp"
#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>>
  linspace(double start, double stop, long num, bool endpoint, bool retstep, dtype d)
  {
    assert(!retstep && "retstep not supported");
    if (num <= 1)
      endpoint = 0;
    double step = 1.;
    if (stop == start || num == 0) // Special case, return [start] if num>0 and [] if num=0
      stop = start + ((num > 0) ? 1 : 0);
    else
      step = (stop - start) / (num - (endpoint ? 1 : 0));
    if (std::is_integral<typename dtype::type>::value)
      return asarray(arange(start, stop + (endpoint ? step * .5 : 0), step), d);
    else
      return arange(start, stop + (endpoint ? step * .5 : 0), step, d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
