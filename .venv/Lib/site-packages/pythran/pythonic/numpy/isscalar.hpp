#ifndef PYTHONIC_NUMPY_ISSCALAR_HPP
#define PYTHONIC_NUMPY_ISSCALAR_HPP

#include "pythonic/include/numpy/isscalar.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"

#include <type_traits>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  constexpr bool isscalar(E const &)
  {
    return types::is_dtype<E>::value || std::is_same<E, types::str>::value;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
