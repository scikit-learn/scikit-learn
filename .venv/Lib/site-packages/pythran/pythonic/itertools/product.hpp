#ifndef PYTHONIC_ITERTOOLS_PRODUCT_HPP
#define PYTHONIC_ITERTOOLS_PRODUCT_HPP

#include "pythonic/include/itertools/product.hpp"
#include "pythonic/itertools/common.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/iterator.hpp"
#include "pythonic/utils/seq.hpp"

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {

    /// product iterator implementation

    template <typename... Iters>
    template <size_t... I>
    product_iterator<Iters...>::product_iterator(std::tuple<Iters...> &_iters,
                                                 std::index_sequence<I...> const &)
        : it_begin(std::get<I>(_iters).begin()...), it_end(std::get<I>(_iters).end()...),
          it(std::get<I>(_iters).begin()...), end(it_begin == it_end)
    {
    }

    template <typename... Iters>
    template <size_t... I>
    product_iterator<Iters...>::product_iterator(npos, std::tuple<Iters...> &_iters,
                                                 std::index_sequence<I...> const &)
        : it_begin(std::get<I>(_iters).end()...), it_end(std::get<I>(_iters).end()...),
          it(std::get<I>(_iters).end()...), end(true)
    {
    }

    template <typename... Iters>
    template <size_t... I>
    types::make_tuple_t<typename Iters::value_type...>
    product_iterator<Iters...>::get_value(std::index_sequence<I...> const &) const
    {
      return types::make_tuple(*std::get<I>(it)...);
    }

    template <typename... Iters>
    types::make_tuple_t<typename Iters::value_type...> product_iterator<Iters...>::operator*() const
    {
      return get_value(std::make_index_sequence<sizeof...(Iters)>{});
    }

    template <typename... Iters>
    template <size_t N>
    void product_iterator<Iters...>::advance(utils::int_<N>)
    {
      if (++std::get<N>(it) == std::get<N>(it_end)) {
        std::get<N>(it) = std::get<N>(it_begin);
        advance(utils::int_<N - 1>());
      }
    }

    template <typename... Iters>
    void product_iterator<Iters...>::advance(utils::int_<0>)
    {
      if (++std::get<0>(it) == std::get<0>(it_end))
        end = true;
    }

    template <typename... Iters>
    product_iterator<Iters...> &product_iterator<Iters...>::operator++()
    {
      advance(utils::int_<sizeof...(Iters) - 1>{});
      return *this;
    }

    template <typename... Iters>
    bool product_iterator<Iters...>::operator==(product_iterator<Iters...> const &other) const
    {
      return end == other.end;
    }

    template <typename... Iters>
    bool product_iterator<Iters...>::operator!=(product_iterator<Iters...> const &other) const
    {
      return end != other.end;
    }

    template <typename... Iters>
    bool product_iterator<Iters...>::operator<(product_iterator<Iters...> const &other) const
    {
      return end != other.end;
    }

    /// details product implementation

    // FIXME: Iterators need to be evaluated as they may be used multiple
    // times
    template <typename... Iters>
    product<Iters...>::product(Iters const &..._iters)
        : utils::iterator_reminder<true, Iters...>(_iters...),
          iterator(this->values, std::make_index_sequence<sizeof...(Iters)>{}),
          end_iter(npos(), this->values, std::make_index_sequence<sizeof...(Iters)>{})
    {
    }

    template <typename... Iters>
    typename product<Iters...>::iterator &product<Iters...>::begin()
    {
      return *this;
    }

    template <typename... Iters>
    typename product<Iters...>::iterator const &product<Iters...>::begin() const
    {
      return *this;
    }

    template <typename... Iters>
    typename product<Iters...>::iterator const &product<Iters...>::end() const
    {
      return end_iter;
    }
  } // namespace details

  template <typename... Iter>
  details::product<std::remove_cv_t<std::remove_reference_t<Iter>>...> product(Iter &&...iters)
  {
    return {std::forward<Iter>(iters)...};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
