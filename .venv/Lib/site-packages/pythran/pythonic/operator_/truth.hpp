#ifndef PYTHONIC_OPERATOR_TRUTH_HPP
#define PYTHONIC_OPERATOR_TRUTH_HPP

#include "pythonic/include/operator_/truth.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  inline bool truth(bool const &a)
  {
    return a;
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
