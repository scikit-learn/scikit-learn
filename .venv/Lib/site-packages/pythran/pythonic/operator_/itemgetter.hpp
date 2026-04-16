#ifndef PYTHONIC_OPERATOR_ITEMGETTER_HPP
#define PYTHONIC_OPERATOR_ITEMGETTER_HPP

#include "pythonic/include/operator_/itemgetter.hpp"

#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{

  inline itemgetter_return::itemgetter_return(long const &item) : i(item)
  {
  }

  template <class A>
  auto itemgetter_return::operator()(A const &a) const -> decltype(a[i])
  {
    return a[i];
  }

  inline itemgetter_return itemgetter(long item)
  {
    return itemgetter_return(item);
  }

  template <typename... Types>
  itemgetter_tuple_return<Types...>::itemgetter_tuple_return(Types... items) : items(items...)
  {
  }

  template <typename... Types>
  itemgetter_tuple_return<Types...>::itemgetter_tuple_return()
  {
  }

  template <typename... Types>
  template <class T, class A, size_t I>
  void itemgetter_tuple_return<Types...>::helper(T &t, A const &a, utils::int_<I>) const
  {
    std::get<I>(t) = a[std::get<I>(items)];
    helper(t, a, utils::int_<I - 1>());
  }

  template <typename... Types>
  template <class T, class A>
  void itemgetter_tuple_return<Types...>::helper(T &t, A const &a, utils::int_<0>) const
  {
    std::get<0>(t) = a[std::get<0>(items)];
  }

  template <typename... Types>
  template <class A>
  auto itemgetter_tuple_return<Types...>::operator()(A const &a) const -> std::tuple<
      std::remove_cv_t<std::remove_reference_t<decltype(a[std::declval<Types>()])>>...>
  {
    std::tuple<std::remove_cv_t<std::remove_reference_t<decltype(a[std::declval<Types>()])>>...> t;
    helper(t, a, utils::int_<sizeof...(Types) - 1>());
    return t;
  }

  template <class... L>
  itemgetter_tuple_return<long, long, L...> itemgetter(long const &item1, long const &item2,
                                                       L... items)
  {
    return {item1, item2, items...};
  }
} // namespace operator_
PYTHONIC_NS_END

#endif
