#ifndef PYTHONIC_OPERATOR_IADD_HPP
#define PYTHONIC_OPERATOR_IADD_HPP

#include "pythonic/include/operator_/iadd.hpp"

#define OPERATOR_NAME iadd
#define OPERATOR_SYMBOL +
#define OPERATOR_ISYMBOL +=

#include "pythonic/operator_/icommon.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/types/set.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A>
  auto iadd(types::empty_list, types::list<A> const &b) -> decltype(b)
  {
    return b;
  }

  template <class K, class V>
  auto iadd(types::empty_dict, types::dict<K, V> const &b) -> decltype(b)
  {
    return b;
  }

  template <class A>
  auto iadd(types::empty_set, types::set<A> const &b) -> decltype(b)
  {
    return b;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
