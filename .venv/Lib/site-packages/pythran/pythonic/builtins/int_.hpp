#ifndef PYTHONIC_BUILTIN_INT_HPP
#define PYTHONIC_BUILTIN_INT_HPP

#include "pythonic/include/builtins/int_.hpp"

#include "pythonic/types/str.hpp"

#include <cassert>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    inline int_::type int_::operator()(char const t[], long base) const
    {
      return std::strtol(t, nullptr, base);
    }
    inline int_::type int_::operator()(types::str const &t, long base) const
    {
      return (*this)(t.c_str(), base);
    }
    inline int_::type int_::operator()(types::chr const &t, long base) const
    {
      char tmp[2] = {t.c, 0};
      return (*this)(&tmp[0], base);
    }

    template <class T>
    int_::type int_::operator()(T &&t) const
    {
      return static_cast<int_::type>(t);
    }

    inline int_::type int_::operator()() const
    {
      return 0L;
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#endif
