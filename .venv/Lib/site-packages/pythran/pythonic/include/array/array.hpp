#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{
  namespace details
  {
    template <char c>
    struct typecodes;
    template <>
    struct typecodes<'b'> {
      using type = signed char;
    };
    template <>
    struct typecodes<'B'> {
      using type = unsigned char;
    };
    template <>
    struct typecodes<'u'> {
      using type = wchar_t;
    };
    template <>
    struct typecodes<'h'> {
      using type = signed short;
    };
    template <>
    struct typecodes<'H'> {
      using type = unsigned short;
    };
    template <>
    struct typecodes<'i'> {
      using type = signed int;
    };
    template <>
    struct typecodes<'I'> {
      using type = unsigned int;
    };
    template <>
    struct typecodes<'l'> {
      using type = signed long;
    };
    template <>
    struct typecodes<'L'> {
      using type = unsigned long;
    };
    template <>
    struct typecodes<'q'> {
      using type = signed long long;
    };
    template <>
    struct typecodes<'Q'> {
      using type = unsigned long long;
    };
    template <>
    struct typecodes<'f'> {
      using type = float;
    };
    template <>
    struct typecodes<'d'> {
      using type = double;
    };

    template <char c>
    types::array<typename details::typecodes<c>::type> array(std::integral_constant<char, c>);

    template <char c, class E>
    types::array<typename details::typecodes<c>::type> array(std::integral_constant<char, c>,
                                                             E &&elts);
  } // namespace details

  DEFINE_FUNCTOR(pythonic::array::details, array);
} // namespace array
PYTHONIC_NS_END

#endif
