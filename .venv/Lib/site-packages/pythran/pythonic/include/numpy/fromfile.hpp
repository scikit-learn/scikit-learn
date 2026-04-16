#ifndef PYTHONIC_INCLUDE_NUMPY_FROMFILE_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMFILE_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromfile(types::str const &file_name, dtype d = dtype(), long count = -1,
           types::str const &sep = {}, long offset = 0);

  DEFINE_FUNCTOR(pythonic::numpy, fromfile);
} // namespace numpy
PYTHONIC_NS_END

#endif
