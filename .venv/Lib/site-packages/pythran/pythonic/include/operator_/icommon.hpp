#ifndef OPERATOR_NAME
#error OPERATOR_NAME ! defined
#endif

#ifndef OPERATOR_SYMBOL
#error OPERATOR_SYMBOL ! defined
#endif

#ifndef OPERATOR_ISYMBOL
#error OPERATOR_ISYMBOL ! defined
#endif

#include "pythonic/utils/functor.hpp"

#ifdef USE_XSIMD
#include <xsimd/xsimd.hpp>
#endif

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  auto OPERATOR_NAME(bool, A &&a, B &&b, ...)
      -> decltype(std::forward<A>(a) OPERATOR_SYMBOL std::forward<B>(b));

  template <class A, class B>
  auto OPERATOR_NAME(bool, A &&a, B &&b, std::nullptr_t)
      -> decltype(std::forward<A>(a) OPERATOR_ISYMBOL std::forward<B>(b));

  template <class A, class B>
  auto OPERATOR_NAME(A &&a, B &&b)
      -> decltype(OPERATOR_NAME(true, std::forward<A>(a), std::forward<B>(b), nullptr))
  {
    return OPERATOR_NAME(true, std::forward<A>(a), std::forward<B>(b), nullptr);
  }

  DEFINE_FUNCTOR(pythonic::operator_, OPERATOR_NAME);
} // namespace operator_
PYTHONIC_NS_END

#undef OPERATOR_NAME
#undef OPERATOR_SYMBOL
#undef OPERATOR_ISYMBOL
