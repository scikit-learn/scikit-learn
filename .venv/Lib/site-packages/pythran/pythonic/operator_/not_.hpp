#ifndef PYTHONIC_OPERATOR_NOT_HPP
#define PYTHONIC_OPERATOR_NOT_HPP

#include "pythonic/include/operator_/not_.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class T>
  auto not_(T &&a) -> decltype(!std::forward<T>(a))
  {
    return !std::forward<T>(a);
  }
  template <class T>
  bool not_(std::complex<T> const &a)
  {
    return !a.real() && !a.imag();
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
