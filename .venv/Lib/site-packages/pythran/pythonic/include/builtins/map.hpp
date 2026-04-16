#ifndef PYTHONIC_INCLUDE_BUILTIN_MAP_HPP
#define PYTHONIC_INCLUDE_BUILTIN_MAP_HPP

#include "pythonic/include/itertools/common.hpp"
#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/iterator.hpp"
#include "pythonic/include/utils/seq.hpp"

#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {

    template <class Operator, class... Iters>
    struct map_res {
      using type = decltype(std::declval<Operator>()(
          std::declval<typename std::iterator_traits<typename Iters::iterator>::value_type>()...));
    };

    template <class... Iters>
    struct map_res<types::none_type, Iters...> {
      using type = decltype(types::make_tuple(
          std::declval<typename std::iterator_traits<typename Iters::iterator>::value_type>()...));
    };

    template <typename Operator, typename... Iters>
    struct map_iterator
        : std::iterator<typename utils::iterator_min<typename Iters::iterator...>::type,
                        typename map_res<Operator, Iters...>::type> {
      std::tuple<typename Iters::iterator...> it;
      Operator _op;

      map_iterator() = default;
      template <size_t... I>
      map_iterator(Operator const &_op, std::tuple<Iters...> &_iters, std::index_sequence<I...>);
      template <size_t... I>
      map_iterator(itertools::npos, Operator const &_op, std::tuple<Iters...> &_iters,
                   std::index_sequence<I...>);

      typename map_res<Operator, Iters...>::type operator*() const;
      map_iterator &operator++();
      map_iterator &operator+=(long i);
      map_iterator operator+(long i) const;
      bool operator==(map_iterator const &other) const;
      bool operator!=(map_iterator const &other) const;
      bool operator<(map_iterator const &other) const;
      bool operator<=(map_iterator const &other) const;
      long operator-(map_iterator const &other) const;

    private:
      template <size_t N>
      long min_len(map_iterator<Operator, Iters...> const &other, utils::int_<N>) const;
      long min_len(map_iterator<Operator, Iters...> const &other, utils::int_<0>) const;

      template <size_t N>
      bool equal(map_iterator const &other, utils::int_<N>) const;
      bool equal(map_iterator const &other, utils::int_<0>) const;

      template <size_t N>
      bool lt(map_iterator const &other, utils::int_<N>) const;
      bool lt(map_iterator const &other, utils::int_<0>) const;

      template <size_t I>
      void advance(long i, utils::int_<I>);
      void advance(long i, utils::int_<0>);

      template <size_t... I>
      void next(std::index_sequence<I...>);

      template <size_t... I>
      typename map_res<Operator, Iters...>::type get_value(std::index_sequence<I...>,
                                                           std::true_type) const;
      template <size_t... I>
      typename map_res<Operator, Iters...>::type get_value(std::index_sequence<I...>,
                                                           std::false_type) const;
    };

    template <typename Operator, typename... Iters>
    struct map : utils::iterator_reminder<true, Iters...>, map_iterator<Operator, Iters...> {
      using iterator = map_iterator<Operator, Iters...>;
      using value_type = typename iterator::value_type;
      using dtype = typename types::dtype_of<value_type>::type;
      static constexpr long value = 1 + utils::nested_container_depth<value_type>::value;

      iterator end_iter;

      map() = default;
      // Use an extra template to enable forwarding
      template <class... Types>
      map(Operator const &_op, Types &&..._iters);

      iterator &begin();
      iterator const &begin() const;
      iterator const &end() const;
    };
  } // namespace details

  template <typename Operator, typename... Iter>
  details::map<std::remove_cv_t<std::remove_reference_t<Operator>>,
               typename types::iterator<std::remove_cv_t<std::remove_reference_t<Iter>>>::type...>
  map(Operator &&_op, Iter &&...iters);

  DEFINE_FUNCTOR(pythonic::builtins, map);
} // namespace builtins

namespace types
{

  template <class Op, class Iter>
  struct len_of<pythonic::builtins::details::map<Op, Iter>> {
    static constexpr long value = len_of<std::remove_cv_t<std::remove_reference_t<Iter>>>::value;
  };

  template <class Op, class I0, class I1, class... Iter>
  struct len_of<pythonic::builtins::details::map<Op, I0, I1, Iter...>> {
    static constexpr long _head = len_of<std::remove_cv_t<std::remove_reference_t<I0>>>::value;
    static constexpr long _tail = len_of<pythonic::builtins::details::map<Op, I1, Iter...>>::value;
    // take the minimal value. If one is negative, it will be automatically
    // selected
    static constexpr long value = (_head < _tail ? _head : _tail);
  };
} // namespace types
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class Op, class... Iter>
struct __combined<E, pythonic::builtins::details::map<Op, Iter...>> {
  using type = typename __combined<
      E, container<typename pythonic::builtins::details::map<Op, Iter...>::value_type>>::type;
};

#endif
