#ifndef PYTHONIC_NUMPY_FLIPLR_HPP
#define PYTHONIC_NUMPY_FLIPLR_HPP

#include "pythonic/include/numpy/fliplr.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto fliplr(E &&expr)
      -> decltype(std::forward<E>(expr)(types::cstride_slice<1>{builtins::None, builtins::None},
                                        types::slice{builtins::None, builtins::None, -1}))
  {
    return std::forward<E>(expr)(types::cstride_slice<1>{builtins::None, builtins::None},
                                 types::slice{builtins::None, builtins::None, -1});
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
