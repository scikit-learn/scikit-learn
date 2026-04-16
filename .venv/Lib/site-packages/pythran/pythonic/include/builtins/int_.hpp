#ifndef PYTHONIC_INCLUDE_BUILTIN_INT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_INT_HPP

#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {

    struct int_ {
      using callable = void;
      using type = long;

      type operator()(char const t[], long base) const;
      type operator()(types::str const &t, long base) const;
      type operator()(types::chr const &t, long base) const;
      template <class T>
      type operator()(T &&t) const;
      type operator()() const;
      friend std::ostream &operator<<(std::ostream &os, int_)
      {
        return os << "int";
      }
    };
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#endif
