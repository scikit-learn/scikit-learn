#ifndef PYTHONIC_BUILTIN_MAP_HPP
#define PYTHONIC_BUILTIN_MAP_HPP

#include "pythonic/include/builtins/map.hpp"

#include "pythonic/itertools/common.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/fwd.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/iterator.hpp"
#include "pythonic/utils/seq.hpp"

#include <iterator>
#include <tuple>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace details
  {

    template <typename Operator, typename... Iters>
    template <size_t... I>
    map_iterator<Operator, Iters...>::map_iterator(Operator const &op, std::tuple<Iters...> &_iters,
                                                   std::index_sequence<I...>)
        : it(std::get<I>(_iters).begin()...), _op(op)
    {
    }

    template <typename Operator, typename... Iters>
    template <size_t... I>
    map_iterator<Operator, Iters...>::map_iterator(itertools::npos, Operator const &op,
                                                   std::tuple<Iters...> &_iters,
                                                   std::index_sequence<I...>)
        : it(std::get<I>(_iters).end()...), _op(op)
    {
    }

    template <typename Operator, typename... Iters>
    template <size_t... I>
    typename map_res<Operator, Iters...>::type
    map_iterator<Operator, Iters...>::get_value(std::index_sequence<I...>, std::false_type) const
    {
      return _op(*std::get<I>(it)...);
    }

    template <typename Operator, typename... Iters>
    template <size_t... I>
    typename map_res<Operator, Iters...>::type
    map_iterator<Operator, Iters...>::get_value(std::index_sequence<I...>, std::true_type) const
    {
      return types::make_tuple(*std::get<I>(it)...);
    }

    template <typename Operator, typename... Iters>
    typename map_res<Operator, Iters...>::type map_iterator<Operator, Iters...>::operator*() const
    {
      return get_value(std::make_index_sequence<sizeof...(Iters)>{},
                       std::is_same<Operator, types::none_type>());
    }

    template <typename Operator, typename... Iters>
    template <size_t... I>
    void map_iterator<Operator, Iters...>::next(std::index_sequence<I...>)
    {
      utils::fwd(++std::get<I>(it)...);
    }

    template <typename Operator, typename... Iters>
    map_iterator<Operator, Iters...> &map_iterator<Operator, Iters...>::operator++()
    {
      next(std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    template <typename Operator, typename... Iters>
    template <size_t I>
    void map_iterator<Operator, Iters...>::advance(long i, utils::int_<I>)
    {
      std::get<I>(it) += i;
      advance(i, utils::int_<I - 1>());
    }

    template <typename Operator, typename... Iters>
    void map_iterator<Operator, Iters...>::advance(long i, utils::int_<0>)
    {
      std::get<0>(it) += i;
    }

    template <typename Operator, typename... Iters>
    map_iterator<Operator, Iters...> &map_iterator<Operator, Iters...>::operator+=(long i)
    {
      advance(i, utils::int_<sizeof...(Iters) - 1>());
      return *this;
    }

    template <typename Operator, typename... Iters>
    map_iterator<Operator, Iters...> map_iterator<Operator, Iters...>::operator+(long i) const
    {
      map_iterator<Operator, Iters...> other(*this);
      other += i;
      return other;
    }

    template <typename Operator, typename... Iters>
    template <size_t N>
    bool map_iterator<Operator, Iters...>::equal(map_iterator<Operator, Iters...> const &other,
                                                 utils::int_<N>) const
    {
      return std::get<N>(other.it) == std::get<N>(it) || equal(other, utils::int_<N - 1>());
    }

    template <typename Operator, typename... Iters>
    bool map_iterator<Operator, Iters...>::equal(map_iterator<Operator, Iters...> const &other,
                                                 utils::int_<0>) const
    {
      return std::get<0>(other.it) == std::get<0>(it);
    }

    template <typename Operator, typename... Iters>
    bool map_iterator<Operator, Iters...>::operator==(
        map_iterator<Operator, Iters...> const &other) const
    {
      return equal(other, utils::int_<sizeof...(Iters) - 1>());
    }

    template <typename Operator, typename... Iters>
    bool map_iterator<Operator, Iters...>::operator!=(
        map_iterator<Operator, Iters...> const &other) const
    {
      return !(*this == other);
    }

    template <typename Operator, typename... Iters>
    template <size_t N>
    bool map_iterator<Operator, Iters...>::lt(map_iterator<Operator, Iters...> const &other,
                                              utils::int_<N>) const
    {
      return std::get<N>(it) < std::get<N>(other.it) ||
             ((std::get<N>(it) == std::get<N>(other.it)) && lt(other, utils::int_<N - 1>()));
    }

    template <typename Operator, typename... Iters>
    bool map_iterator<Operator, Iters...>::lt(map_iterator<Operator, Iters...> const &other,
                                              utils::int_<0>) const
    {
      return std::get<0>(it) < std::get<0>(other.it);
    }

    template <typename Operator, typename... Iters>
    bool
    map_iterator<Operator, Iters...>::operator<(map_iterator<Operator, Iters...> const &other) const
    {
      return lt(other, utils::int_<sizeof...(Iters) - 1>());
    }

    template <typename Operator, typename... Iters>
    bool map_iterator<Operator, Iters...>::operator<=(
        map_iterator<Operator, Iters...> const &other) const
    {
      return (*this == other) || (*this < other);
    }

    template <typename Operator, typename... Iters>
    template <size_t N>
    long map_iterator<Operator, Iters...>::min_len(map_iterator<Operator, Iters...> const &other,
                                                   utils::int_<N>) const
    {
      return std::min((long)(std::get<N>(it) - std::get<N>(other.it)),
                      min_len(other, utils::int_<N - 1>()));
    }

    template <typename Operator, typename... Iters>
    long map_iterator<Operator, Iters...>::min_len(map_iterator<Operator, Iters...> const &other,
                                                   utils::int_<0>) const
    {
      return std::get<0>(it) - std::get<0>(other.it);
    }

    template <typename Operator, typename... Iters>
    long
    map_iterator<Operator, Iters...>::operator-(map_iterator<Operator, Iters...> const &other) const
    {
      return min_len(other, utils::int_<sizeof...(Iters) - 1>());
    }

    template <typename Operator, typename... Iters>
    template <class... Types>
    map<Operator, Iters...>::map(Operator const &_op, Types &&..._iters)
        : utils::iterator_reminder<true, Iters...>(std::forward<Types>(_iters)...),
          map_iterator<Operator, Iters...>(_op, this->values,
                                           std::make_index_sequence<sizeof...(Iters)>{}),
          end_iter(itertools::npos(), _op, this->values,
                   std::make_index_sequence<sizeof...(Iters)>{})
    {
    }

    template <typename Operator, typename... Iters>
    typename map<Operator, Iters...>::iterator &map<Operator, Iters...>::begin()
    {
      return *this;
    }

    template <typename Operator, typename... Iters>
    typename map<Operator, Iters...>::iterator const &map<Operator, Iters...>::begin() const
    {
      return *this;
    }

    template <typename Operator, typename... Iters>
    typename map<Operator, Iters...>::iterator const &map<Operator, Iters...>::end() const
    {
      return end_iter;
    }
  } // namespace details

  template <typename Operator, typename... Iter>
  details::map<std::remove_cv_t<std::remove_reference_t<Operator>>,
               typename types::iterator<std::remove_cv_t<std::remove_reference_t<Iter>>>::type...>
  map(Operator &&_op, Iter &&...iters)
  {
    return {std::forward<Operator>(_op), std::forward<Iter>(iters)...};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
