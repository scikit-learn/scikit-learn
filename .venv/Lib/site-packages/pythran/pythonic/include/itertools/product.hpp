#ifndef PYTHONIC_INCLUDE_ITERTOOLS_PRODUCT_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_PRODUCT_HPP

#include "pythonic/include/itertools/common.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/iterator.hpp"
#include "pythonic/include/utils/seq.hpp"

#include <iterator>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {

    // FIXME : should be a combined_iterator_tag
    template <typename... Iters>
    struct product_iterator : std::iterator<std::forward_iterator_tag,
                                            types::make_tuple_t<typename Iters::value_type...>> {

      std::tuple<typename Iters::iterator...> const it_begin;
      std::tuple<typename Iters::iterator...> const it_end;
      std::tuple<typename Iters::iterator...> it;
      bool end;

      product_iterator() = default;
      template <size_t... I>
      product_iterator(std::tuple<Iters...> &_iters, std::index_sequence<I...> const &);
      template <size_t... I>
      product_iterator(npos, std::tuple<Iters...> &_iters, std::index_sequence<I...> const &);
      types::make_tuple_t<typename Iters::value_type...> operator*() const;
      product_iterator &operator++();
      bool operator==(product_iterator const &other) const;
      bool operator!=(product_iterator const &other) const;
      bool operator<(product_iterator const &other) const;

    private:
      template <size_t N>
      void advance(utils::int_<N>);
      void advance(utils::int_<0>);
      template <size_t... I>
      types::make_tuple_t<typename Iters::value_type...>
      get_value(std::index_sequence<I...> const &) const;
    };

    template <typename... Iters>
    struct product : utils::iterator_reminder<true, Iters...>, product_iterator<Iters...> {

      using value_type = types::make_tuple_t<typename Iters::value_type...>;
      using iterator = product_iterator<Iters...>;

      iterator end_iter;

      product() = default;
      product(Iters const &..._iters);

      iterator &begin();
      iterator const &begin() const;
      iterator const &end() const;
    };
  } // namespace details

  template <typename... Iter>
  details::product<std::remove_cv_t<std::remove_reference_t<Iter>>...> product(Iter &&...iters);

  DEFINE_FUNCTOR(pythonic::itertools, product);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class... Iter>
struct __combined<E, pythonic::itertools::details::product<Iter...>> {
  using type = typename __combined<
      E, container<typename pythonic::itertools::details::product<Iter...>::value_type>>::type;
};

/* } */

#endif
