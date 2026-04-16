#ifndef PYTHONIC_INCLUDE_BUILTIN_REVERSED_HPP
#define PYTHONIC_INCLUDE_BUILTIN_REVERSED_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class Iterable>
    struct reversed {

      using value_type = typename Iterable::value_type;
      using iterator = typename Iterable::reverse_iterator;
      using const_iterator = typename Iterable::const_reverse_iterator;

      Iterable iterable;

      reversed();
      reversed(Iterable const &iterable);
      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;
    };
  } // namespace details

  template <class Iterable>
  details::reversed<Iterable> reversed(Iterable const &iterable);

  DEFINE_FUNCTOR(pythonic::builtins, reversed);
} // namespace builtins
PYTHONIC_NS_END

#endif
