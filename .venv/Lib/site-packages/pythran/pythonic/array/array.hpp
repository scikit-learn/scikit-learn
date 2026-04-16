#ifndef PYTHONIC_ARRAY_ARRAY_HPP
#define PYTHONIC_ARRAY_ARRAY_HPP

#include "pythonic/include/array/array.hpp"

#include "pythonic/types/array.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{
  namespace details
  {

    template <char c>
    types::array<typename details::typecodes<c>::type> array(std::integral_constant<char, c>)
    {
      return {0};
    }

    template <char c, class E>
    types::array<typename details::typecodes<c>::type> array(std::integral_constant<char, c>,
                                                             E &&elts)
    {
      return {std::forward<E>(elts)};
    }
  } // namespace details

} // namespace array
PYTHONIC_NS_END

#endif
