#ifndef PYTHONIC_NUMPY_ISSCTYPE_HPP
#define PYTHONIC_NUMPY_ISSCTYPE_HPP

#include "pythonic/include/numpy/issctype.hpp"

#include "pythonic/numpy/isscalar.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  constexpr auto issctype(E const &expr)
      -> std::enable_if_t<!types::is_dtype<E>::value && !std::is_same<E, types::str>::value, bool>
  {
    return isscalar(typename E::type());
  }

  template <class E>
  constexpr auto issctype(E const &expr)
      -> std::enable_if_t<types::is_dtype<E>::value || std::is_same<E, types::str>::value, bool>
  {
    return false;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
