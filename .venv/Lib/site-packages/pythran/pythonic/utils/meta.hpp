#ifndef PYTHONIC_UTILS_META_HPP
#define PYTHONIC_UTILS_META_HPP

#include "pythonic/include/utils/meta.hpp"

template <bool C, class... Types>
struct static_assert_check {
  static_assert(C, "Assertion failed <see below for more information>");
  static constexpr bool value = C;
};

#define pythran_static_assert(value, str, ...)                                                     \
  static_assert(static_assert_check<value, __VA_ARGS__>::value, str)

#endif
