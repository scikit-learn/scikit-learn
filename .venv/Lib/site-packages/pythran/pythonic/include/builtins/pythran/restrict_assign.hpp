#ifndef PYTHONIC_INCLUDE_BUILTIN_RESTRICT_ASSIGN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_RESTRICT_ASSIGN_HPP

#include "pythonic/include/types/numpy_gexpr.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace pythran
  {
    template <class Arg, class... S, class E>
    void restrict_assign(types::numpy_gexpr<Arg, S...> &&target, E &&value)
    {
      std::move(target)._copy_restrict(std::forward<E>(value));
    }

    template <class T, class E>
    void restrict_assign(T &&target, E &&value)
    {
      target = std::forward<E>(value);
    }

    DEFINE_FUNCTOR(pythonic::builtins::pythran, restrict_assign);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
