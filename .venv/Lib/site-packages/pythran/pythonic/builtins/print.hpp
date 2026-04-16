#ifndef PYTHONIC_BUILTIN_PRINT_HPP
#define PYTHONIC_BUILTIN_PRINT_HPP

#include "pythonic/include/builtins/print.hpp"
#include "pythonic/utils/functor.hpp"

#include <iostream>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class T>
    std::ostream &print(std::ostream &os, T const &t)
    {
      return os << t;
    }

    inline std::ostream &print(std::ostream &os, bool const &t)
    {
      static char const repr[2][6] = {"False", "True\0"};
      return os << repr[t];
    }
  } // namespace details

  inline void print_nonl()
  {
  }

  template <typename T, typename... Types>
  void print_nonl(T const &value, Types const &...values)
  {
    details::print(std::cout, value);
    if (sizeof...(Types) > 0)
      std::cout << ' ';
    print_nonl(values...);
  }

  inline void print()
  {
    std::cout << std::endl;
  }

  template <typename T, typename... Types>
  void print(T const &value, Types const &...values)
  {
    details::print(std::cout, value);
    if (sizeof...(values) > 0)
      std::cout << ' ';
    print(values...);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
