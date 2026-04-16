#ifndef PYTHONIC_TYPES_VARIANT_FUNCTOR_HPP
#define PYTHONIC_TYPES_VARIANT_FUNCTOR_HPP

#include "pythonic/include/types/variant_functor.hpp"
#include "pythonic/utils/meta.hpp"

#include <cassert>
#include <utility>

PYTHONIC_NS_BEGIN

namespace types
{

  namespace details
  {

    template <class Type>
    variant_functor_impl<Type>::variant_functor_impl(char mem[], Type const &t)
        : fun(new(mem) Type(t))
    {
    }

    template <class Type>
    template <class OtherType>
    variant_functor_impl<Type>::variant_functor_impl(char mem[], OtherType const &t) : fun(nullptr)
    {
    }

    template <class Type>
    variant_functor_impl<Type>::variant_functor_impl(char mem[],
                                                     variant_functor_impl<Type> const &t)
        : fun(t.fun ? new(mem) Type(*t.fun) : nullptr)
    {
    }

    template <class Type>
    template <class... OtherTypes>
    variant_functor_impl<Type>::variant_functor_impl(
        char mem[], variant_functor_impl<Type, OtherTypes...> const &t)
        : variant_functor_impl(mem, t.head)
    {
    }

    template <class Type>
    template <class OtherType, class... OtherTypes>
    variant_functor_impl<Type>::variant_functor_impl(
        char mem[], variant_functor_impl<OtherType, OtherTypes...> const &t)
        : variant_functor_impl(mem, t.tail)
    {
    }

    template <class Type>
    void variant_functor_impl<Type>::assign(char mem[], variant_functor_impl<Type> const &other)
    {
      if (fun != nullptr)
        fun->~Type();
      if (other.fun)
        fun = new (mem) Type(*other.fun);
    }
    template <class Type>
    void variant_functor_impl<Type>::assign(char mem[], variant_functor<Type> const &other)
    {
      assign(mem, static_cast<variant_functor_impl<Type> const &>(other));
    }

    template <class Type>
    template <class OT0, class OT1, class... OtherTypes>
    void
    variant_functor_impl<Type>::assign(char mem[],
                                       variant_functor_impl<OT0, OT1, OtherTypes...> const &other)
    {
      assign(mem, other.head);
      assign(mem, other.tail);
    }

    template <class Type>
    template <class OT0, class OT1, class... OtherTypes>
    void variant_functor_impl<Type>::assign(char mem[],
                                            variant_functor<OT0, OT1, OtherTypes...> const &other)
    {
      assign(mem, static_cast<variant_functor_impl<OT0, OT1, OtherTypes...> const &>(other));
    }

    template <class Type>
    template <class OtherType>
    void variant_functor_impl<Type>::assign(char mem[],
                                            variant_functor_impl<OtherType> const &other)
    {
    }

    template <class Type>
    template <class OtherType>
    void variant_functor_impl<Type>::assign(char mem[], variant_functor<OtherType> const &other)
    {
    }

    template <class Type>
    void variant_functor_impl<Type>::assign(char mem[], Type const &other)
    {
      if (fun != nullptr)
        fun->~Type();
      fun = new (mem) Type(other);
    }

    template <class Type>
    variant_functor_impl<Type>::~variant_functor_impl()
    {
      if (fun != nullptr)
        fun->~Type();
    }

    template <class Type>
    template <class OtherType>
    void variant_functor_impl<Type>::assign(char mem[], OtherType const &other)
    {
    }

    template <class Type>
    template <class... Args>
    auto variant_functor_impl<Type>::operator()(Args &&...args)
        -> decltype(std::declval<Type>()(std::forward<Args>(args)...))
    {
      assert(fun && "handler defined");
      return (*fun)(std::forward<Args>(args)...);
    }

    template <class Type>
    template <class... Args>
    auto variant_functor_impl<Type>::operator()(Args &&...args) const
        -> decltype(std::declval<Type>()(std::forward<Args>(args)...))
    {
      assert(fun && "handler defined");
      return (*fun)(std::forward<Args>(args)...);
    }

    template <class Type, class... Types>
    template <class... OtherTypes>
    variant_functor_impl<Type, Types...>::variant_functor_impl(char mem[], OtherTypes const &...t)
        : head(mem, t...), tail(mem, t...)
    {
    }

    template <class Type, class... Types>
    template <class... OtherTypes>
    variant_functor_impl<Type, Types...>::variant_functor_impl(
        char mem[], variant_functor_impl<OtherTypes...> const &t)
        : head(mem, t), tail(mem, t)
    {
    }

    template <class Type, class... Types>
    void
    variant_functor_impl<Type, Types...>::assign(char mem[],
                                                 variant_functor_impl<Type, Types...> const &other)
    {
      head.assign(mem, other);
      tail.assign(mem, other);
    }

    template <class Type, class... Types>
    template <class OtherType>
    void variant_functor_impl<Type, Types...>::assign(char mem[], OtherType const &other)
    {
      head.assign(mem, other);
      tail.assign(mem, other);
    }

    template <class Type, class... Types>
    template <class... Args>
    auto variant_functor_impl<Type, Types...>::operator()(Args &&...args) ->
        typename __combined<decltype(std::declval<Type>()(std::forward<Args>(args)...)),
                            decltype(std::declval<Types>()(std::forward<Args>(args)...))...>::type
    {
      if (head.fun)
        return head(std::forward<Args>(args)...);
      else
        return tail(std::forward<Args>(args)...);
    }

    template <class Type, class... Types>
    template <class... Args>
    auto variant_functor_impl<Type, Types...>::operator()(Args &&...args) const ->
        typename __combined<decltype(std::declval<Type>()(std::forward<Args>(args)...)),
                            decltype(std::declval<Types>()(std::forward<Args>(args)...))...>::type
    {
      if (head.fun)
        return head(std::forward<Args>(args)...);
      else
        return tail(std::forward<Args>(args)...);
    }
  } // namespace details

  template <class... Types>
  variant_functor<Types...>::variant_functor(variant_functor const &other)
      : details::variant_functor_impl<Types...>(
            mem, static_cast<details::variant_functor_impl<Types...> const &>(other))
  {
  }

  template <class... Types>
  variant_functor<Types...> &
  variant_functor<Types...>::operator=(variant_functor<Types...> const &other)
  {
    details::variant_functor_impl<Types...>::assign(mem, other);
    return *this;
  }

  template <class... Types>
  template <class... OtherTypes>
  variant_functor<Types...> &
  variant_functor<Types...>::operator=(variant_functor<OtherTypes...> const &other)
  {
    details::variant_functor_impl<Types...>::assign(mem, other);
    return *this;
  }

  template <class... Types>
  template <class OtherType>
  variant_functor<Types...> &variant_functor<Types...>::operator=(OtherType const &other)
  {
    static_assert(utils::any_of<std::is_same<OtherType, Types>::value...>::value,
                  "consistent assign");
    details::variant_functor_impl<Types...>::assign(mem, other);
    return *this;
  }

  template <class... Types>
  template <class... OtherTypes>
  variant_functor<Types...>::variant_functor(OtherTypes const &...t)
      : details::variant_functor_impl<Types...>(mem, t...)
  {
  }

  template <class... Types>
  template <class... OtherTypes>
  variant_functor<Types...>::variant_functor(variant_functor<OtherTypes...> const &t)
      : details::variant_functor_impl<Types...>(
            mem, static_cast<details::variant_functor_impl<OtherTypes...> const &>(t))
  {
  }
} // namespace types
PYTHONIC_NS_END
#endif
