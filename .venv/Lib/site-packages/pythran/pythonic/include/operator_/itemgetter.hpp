#ifndef PYTHONIC_INCLUDE_OPERATOR_ITEMGETTER_HPP
#define PYTHONIC_INCLUDE_OPERATOR_ITEMGETTER_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"

PYTHONIC_NS_BEGIN

namespace operator_
{
  struct itemgetter_return {
    long i;
    itemgetter_return(long const &item = -1);
    template <class A>
    auto operator()(A const &a) const -> decltype(a[i]);
  };

  itemgetter_return itemgetter(long item);

  template <typename... Types>
  struct itemgetter_tuple_return {

    std::tuple<Types...> items;

    itemgetter_tuple_return(Types... items);

    itemgetter_tuple_return();

    template <class T, class A, size_t I>
    void helper(T &t, A const &a, utils::int_<I>) const;

    template <class T, class A>
    void helper(T &t, A const &a, utils::int_<0>) const;

    template <class A>
    auto operator()(A const &a) const -> std::tuple<
        std::remove_cv_t<std::remove_reference_t<decltype(a[std::declval<Types>()])>>...>;
  };

  template <class... L>
  itemgetter_tuple_return<long, long, L...> itemgetter(long const &item1, long const &item2,
                                                       L... items);

  DEFINE_FUNCTOR(pythonic::operator_, itemgetter);
} // namespace operator_
PYTHONIC_NS_END

#endif
