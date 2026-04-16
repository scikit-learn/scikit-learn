#ifndef PYTHONIC_INCLUDE_TYPES_TYPE_HPP
#define PYTHONIC_INCLUDE_TYPES_TYPE_HPP

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct type_functor;

  template <class T>
  using type_t = typename type_functor<T>::type;
} // namespace types

PYTHONIC_NS_END

#endif
