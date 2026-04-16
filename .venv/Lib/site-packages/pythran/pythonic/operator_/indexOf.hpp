#ifndef PYTHONIC_OPERATOR_INDEXOF_HPP
#define PYTHONIC_OPERATOR_INDEXOF_HPP

#include "pythonic/include/operator_/indexOf.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/builtins/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace operator_
{

  template <class A, class B>
  long indexOf(A &&a, B &&b)
  {
    auto where = std::find(a.begin(), a.end(), b);
    if (where == a.end())
      throw types::ValueError(builtins::anonymous::str(b) + " is not in this sequence");
    return where - a.begin();
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
