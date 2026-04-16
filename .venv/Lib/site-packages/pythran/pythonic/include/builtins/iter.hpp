#ifndef PYTHONIC_INCLUDE_BUILTIN_ITER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ITER_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class T>
    struct iter : T::iterator {
      using iterator = typename T::iterator;

      iterator _end;
      T data;

      iter();
      iter(T data);
      iterator &begin();
      iterator const &begin() const;
      iterator const &end() const;
    };
  } // namespace details

  template <class T>
  details::iter<std::remove_cv_t<std::remove_reference_t<T>>> iter(T &&t);

  DEFINE_FUNCTOR(pythonic::builtins, iter);
} // namespace builtins
PYTHONIC_NS_END

#endif
