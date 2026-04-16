#ifndef PYTHONIC_INCLUDE_NUMPY_FLIP_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLIP_HPP

#include "pythonic/include/types/numpy_gexpr.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class E, class S, size_t... I>
    auto flip(E const &expr, S const &slices, std::index_sequence<I...>)
        -> decltype(expr(slices[I]...));
  }

  template <class E>
  auto flip(E const &expr, long axis)
      -> decltype(details::flip(expr, std::array<types::slice, E::value>{},
                                std::make_index_sequence<E::value>{}));

  DEFINE_FUNCTOR(pythonic::numpy, flip);
} // namespace numpy
PYTHONIC_NS_END

#endif
