#ifndef PYTHONIC_INCLUDE_BUILTIN_PRINT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PRINT_HPP

#include "pythonic/include/utils/functor.hpp"
#include <ostream>

PYTHONIC_NS_BEGIN

namespace builtins
{

  void print_nonl();

  template <typename T, typename... Types>
  void print_nonl(T const &value, Types const &...values);

  void print();

  template <typename T, typename... Types>
  void print(T const &value, Types const &...values);
  DEFINE_FUNCTOR(pythonic::builtins, print);
} // namespace builtins
PYTHONIC_NS_END

#endif
