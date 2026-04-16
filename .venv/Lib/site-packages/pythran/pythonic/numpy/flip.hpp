#ifndef PYTHONIC_NUMPY_FLIP_HPP
#define PYTHONIC_NUMPY_FLIP_HPP

#include "pythonic/include/numpy/flip.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class E, class S, size_t... I>
    auto flip(E const &expr, S const &slices, std::index_sequence<I...>)
        -> decltype(expr(slices[I]...))
    {
      return expr(slices[I]...);
    }
  } // namespace details

  template <class E>
  auto flip(E const &expr, long axis)
      -> decltype(details::flip(expr, std::array<types::slice, E::value>{},
                                std::make_index_sequence<E::value>{}))
  {
    std::array<types::slice, E::value> slices;
    slices[axis].step = -1;
    return details::flip(expr, slices, std::make_index_sequence<E::value>{});
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
