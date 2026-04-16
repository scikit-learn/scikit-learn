#ifndef PYTHONIC_INCLUDE_OPERATOR_ICONCAT_HPP
#define PYTHONIC_INCLUDE_OPERATOR_ICONCAT_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  template <class A, class B>
  A iconcat(A a, B const &b);

  template <class A>
  auto iconcat(types::empty_list a, types::list<A> b) -> decltype(b);

  template <class K, class V>
  auto iconcat(types::empty_dict a, types::dict<K, V> b) -> decltype(b);

  template <class A>
  auto iconcat(types::empty_set a, types::set<A> b) -> decltype(b);

  DEFINE_FUNCTOR(pythonic::operator_, iconcat);
} // namespace operator_
PYTHONIC_NS_END

#endif
