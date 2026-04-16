#ifndef PYTHONIC_INCLUDE_OPERATOR_NOT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_NOT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class T>
  auto not_(T &&a) -> decltype(!std::forward<T>(a));

  template <class T>
  bool not_(std::complex<T> const &a);

  DEFINE_FUNCTOR(pythonic::operator_, not_);
} // namespace operator_
PYTHONIC_NS_END

#endif
