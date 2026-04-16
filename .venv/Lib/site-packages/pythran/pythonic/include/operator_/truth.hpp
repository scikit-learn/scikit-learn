#ifndef PYTHONIC_INCLUDE_OPERATOR_TRUTH_HPP
#define PYTHONIC_INCLUDE_OPERATOR_TRUTH_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  bool truth(bool const &a);

  DEFINE_FUNCTOR(pythonic::operator_, truth);
} // namespace operator_
PYTHONIC_NS_END

#endif
