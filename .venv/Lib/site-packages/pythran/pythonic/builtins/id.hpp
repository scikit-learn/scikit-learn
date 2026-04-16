#ifndef PYTHONIC_BUILTIN_ID_HPP
#define PYTHONIC_BUILTIN_ID_HPP

#include "pythonic/include/builtins/id.hpp"

#include "pythonic/utils/functor.hpp"

/*
 * We use uintptr_t conversion because on windows 64 bits, sizeof(void*) == 8
 * && sizeof(long) == 4. Because of this, void* to long is forbidden but
 * void* -> uintptr_t -> long is allowed
 * Accuracy is lost this way...
 */

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  long id(T const &t)
  {
    return t.id();
  }

  inline long id(long const &t)
  {
    return reinterpret_cast<uintptr_t>(&t);
  }

  inline long id(double const &t)
  {
    return reinterpret_cast<uintptr_t>(&t);
  }

  inline long id(bool const &t)
  {
    return reinterpret_cast<uintptr_t>(&t);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
