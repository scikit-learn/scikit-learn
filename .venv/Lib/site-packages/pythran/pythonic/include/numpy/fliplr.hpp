#ifndef PYTHONIC_INCLUDE_NUMPY_FLIPLR_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLIPLR_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto fliplr(E &&expr)
      -> decltype(std::forward<E>(expr)(types::cstride_slice<1>{builtins::None, builtins::None},
                                        types::slice{builtins::None, builtins::None, -1}));

  DEFINE_FUNCTOR(pythonic::numpy, fliplr);
} // namespace numpy
PYTHONIC_NS_END

#endif
