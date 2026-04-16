#ifndef PYTHONIC_INCLUDE_STRING_FIND_HPP
#define PYTHONIC_INCLUDE_STRING_FIND_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace string
{

  template <class T>
  long find(types::str const &s, T &&val);

  DEFINE_FUNCTOR(pythonic::string, find);
} // namespace string
PYTHONIC_NS_END

#endif
