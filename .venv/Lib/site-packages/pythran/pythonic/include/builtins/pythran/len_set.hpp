#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_LEN_SET_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_LEN_SET_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {

    template <class Iterable>
    long len_set(Iterable const &s);

    DEFINE_FUNCTOR(pythonic::builtins::pythran, len_set);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
