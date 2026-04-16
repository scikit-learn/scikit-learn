#ifndef PYTHONIC_TYPES_NUMPY_EXPR_HPP
#define PYTHONIC_TYPES_NUMPY_EXPR_HPP

#include "pythonic/include/types/numpy_expr.hpp"

#include "pythonic/types/nditerator.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/meta.hpp"

#include "pythonic/builtins/ValueError.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  namespace details
  {

    inline long best_of()
    {
      return 1;
    }
    template <class V0, class... Vs>
    long best_of(V0 v0, Vs... vs)
    {
      long vtail = best_of(vs...);
      return ((long)v0 == vtail) ? v0 : (v0 * vtail);
    }

    template <size_t I, class Args, size_t... Is>
    long init_shape_element(Args const &args, std::index_sequence<Is...>)
    {
      return best_of(std::get<Is>(args).template shape<I>()...);
    }
  } // namespace details

  template <class Op, class... Args>
  numpy_expr<Op, Args...>::numpy_expr(Args const &...args) : args(args...)
  {
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::const_iterator
  numpy_expr<Op, Args...>::_begin(std::index_sequence<I...>) const
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            const_cast<std::decay_t<Args> const &>(std::get<I>(args)).begin()...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::const_iterator numpy_expr<Op, Args...>::begin() const
  {
    return _begin(std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::const_iterator
  numpy_expr<Op, Args...>::_end(std::index_sequence<I...>) const
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            const_cast<std::decay_t<Args> const &>(std::get<I>(args)).end()...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::const_iterator numpy_expr<Op, Args...>::end() const
  {
    return _end(std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::const_fast_iterator
  numpy_expr<Op, Args...>::begin(types::fast) const
  {
    return {*this, 0};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::const_fast_iterator
  numpy_expr<Op, Args...>::end(types::fast) const
  {
    return {*this, size()};
  }

  template <class Op, class... Args>
  template <size_t... I>
  bool numpy_expr<Op, Args...>::_no_broadcast(std::index_sequence<I...>) const
  {
    bool child_broadcast = false;
    (void)std::initializer_list<bool>{
        (child_broadcast |= !utils::no_broadcast(std::get<I>(args)))...};
    if (child_broadcast)
      return false;

    bool same_shape = true;
    (void)std::initializer_list<bool>{
        (same_shape &= (is_trivial_broadcast<I, decltype(args)>() ||
                        std::get<I>(args).template shape<0>() == size()))...};
    return same_shape;
  }

  template <class Op, class... Args>
  template <size_t... I>
  bool numpy_expr<Op, Args...>::_no_broadcast_ex(std::index_sequence<I...>) const
  {
    bool child_broadcast = false;
    (void)std::initializer_list<bool>{
        (child_broadcast |= !utils::no_broadcast_ex(std::get<I>(args)))...};
    if (child_broadcast)
      return false;

    bool same_shape = true;
    auto shp = sutils::getshape(*this);
    (void)std::initializer_list<bool>{
        (same_shape &= (is_trivial_broadcast<I, decltype(args)>() ||
                        sutils::getshape(std::get<I>(args)) == shp))...};
    return same_shape;
  }

  template <class Op, class... Args>
  template <size_t... I>
  bool numpy_expr<Op, Args...>::_no_broadcast_vectorize(std::index_sequence<I...>) const
  {
    bool child_broadcast = false;
    (void)std::initializer_list<bool>{
        (child_broadcast |= !utils::no_broadcast_vectorize(std::get<I>(args)))...};
    if (child_broadcast)
      return false;

    bool same_shape = true;
    (void)std::initializer_list<bool>{
        (same_shape &= ((long)std::get<I>(args).template shape<0>() == size()))...};
    return same_shape;
  }

  template <class Op, class... Args>
  bool numpy_expr<Op, Args...>::no_broadcast() const
  {
    return _no_broadcast(std::make_index_sequence<sizeof...(Args)>{});
  }
  template <class Op, class... Args>
  bool numpy_expr<Op, Args...>::no_broadcast_ex() const
  {
    return _no_broadcast_ex(std::make_index_sequence<sizeof...(Args)>{});
  }
  template <class Op, class... Args>
  bool numpy_expr<Op, Args...>::no_broadcast_vectorize() const
  {
    return _no_broadcast_vectorize(std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::iterator
  numpy_expr<Op, Args...>::_begin(std::index_sequence<I...>)
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            const_cast<std::decay_t<Args> &>(std::get<I>(args)).begin()...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::iterator numpy_expr<Op, Args...>::begin()
  {
    return _begin(std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::iterator
  numpy_expr<Op, Args...>::_end(std::index_sequence<I...>)
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            const_cast<std::decay_t<Args> &>(std::get<I>(args)).end()...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::iterator numpy_expr<Op, Args...>::end()
  {
    return _end(std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  auto numpy_expr<Op, Args...>::fast(long i) const
      -> decltype(this->_fast(i, std::make_index_sequence<sizeof...(Args)>{}))
  {
    return _fast(i, std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <class... Indices>
  auto numpy_expr<Op, Args...>::map_fast(Indices... indices) const
      -> decltype(this->_map_fast(array_tuple<long, sizeof...(Indices)>{{indices...}},
                                  std::make_index_sequence<sizeof...(Args)>{}))
  {
    static_assert(sizeof...(Indices) == sizeof...(Args), "compatible call");
    return _map_fast(array_tuple<long, sizeof...(Indices)>{{indices...}},
                     std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  auto numpy_expr<Op, Args...>::operator[](long i) const -> decltype(this->fast(i))
  {
    if (i < 0)
      i += size();
    return fast(i);
  }

#ifdef USE_XSIMD
  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::simd_iterator
  numpy_expr<Op, Args...>::_vbegin(vectorize, std::index_sequence<I...>) const
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            std::make_tuple(xsimd::batch<typename std::remove_reference_t<Args>::value_type>(
                *std::get<I>(args).begin())...),
            std::get<I>(args).vbegin(vectorize{})...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::simd_iterator numpy_expr<Op, Args...>::vbegin(vectorize) const
  {
    return _vbegin(vectorize{}, std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::simd_iterator
  numpy_expr<Op, Args...>::_vend(vectorize, std::index_sequence<I...>) const
  {
    return {{make_step(size(), std::get<I>(args).template shape<0>())...},
            {},
            std::get<I>(args).vend(vectorize{})...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::simd_iterator numpy_expr<Op, Args...>::vend(vectorize) const
  {
    return _vend(vectorize{}, std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::simd_iterator_nobroadcast
  numpy_expr<Op, Args...>::_vbegin(vectorize_nobroadcast, std::index_sequence<I...>) const
  {
    return {std::get<I>(args).vbegin(vectorize_nobroadcast{})...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::simd_iterator_nobroadcast
  numpy_expr<Op, Args...>::vbegin(vectorize_nobroadcast) const
  {
    return _vbegin(vectorize_nobroadcast{}, std::make_index_sequence<sizeof...(Args)>{});
  }

  template <class Op, class... Args>
  template <size_t... I>
  typename numpy_expr<Op, Args...>::simd_iterator_nobroadcast
  numpy_expr<Op, Args...>::_vend(vectorize_nobroadcast, std::index_sequence<I...>) const
  {
    return {std::get<I>(args).vend(vectorize_nobroadcast{})...};
  }

  template <class Op, class... Args>
  typename numpy_expr<Op, Args...>::simd_iterator_nobroadcast
  numpy_expr<Op, Args...>::vend(vectorize_nobroadcast) const
  {
    return _vend(vectorize_nobroadcast{}, std::make_index_sequence<sizeof...(Args)>{});
  }

#endif

  template <class Op, class... Args>
  template <class... S>
  auto numpy_expr<Op, Args...>::operator()(S const &...s) const
      -> decltype(this->_get(std::make_index_sequence<sizeof...(Args)>{}, s...))
  {
    return _get(std::make_index_sequence<sizeof...(Args)>{}, s...);
  }

  template <class Op, class... Args>
  template <class F>
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       !is_pod_array<F>::value,
                   numpy_vexpr<numpy_expr<Op, Args...>, ndarray<long, pshape<long>>>>
  numpy_expr<Op, Args...>::fast(F const &filter) const
  {
    long sz = filter.template shape<0>();
    long *raw = utils::allocate<long>(sz);
    long n = 0;
    for (long i = 0; i < sz; ++i)
      if (filter.fast(i))
        raw[n++] = i;
    // reallocate(raw, n);
    long shp[1] = {n};
    return this->fast(ndarray<long, pshape<long>>(raw, shp, types::ownership::owned));
  }

  template <class Op, class... Args>
  template <class F>
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       !is_pod_array<F>::value,
                   numpy_vexpr<numpy_expr<Op, Args...>, ndarray<long, pshape<long>>>>
  numpy_expr<Op, Args...>::operator[](F const &filter) const
  {
    return fast(filter);
  }
  template <class Op, class... Args>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                       !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_expr<Op, Args...>, F>>
  numpy_expr<Op, Args...>::operator[](F const &filter) const
  {
    return {*this, filter};
  }

  template <class Op, class... Args>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                       !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_expr<Op, Args...>, F>>
  numpy_expr<Op, Args...>::fast(F const &filter) const
  {
    return {*this, filter};
  }

  template <class Op, class... Args>
  numpy_expr<Op, Args...>::operator bool() const
  {
    if (sutils::any_of(*this, [](long n) { return n != 1; }))
      throw ValueError("The truth value of an array with more than one element "
                       "is ambiguous. Use a.any() or a.all()");
    array_tuple<long, value> first = {0};
    return operator[](first);
  }

  template <class Op, class... Args>
  long numpy_expr<Op, Args...>::flat_size() const
  {
    return prod_helper(*this, std::make_index_sequence<value>());
  }

  template <class Op, class... Args>
  long numpy_expr<Op, Args...>::size() const
  {
    return this->template shape<0>();
  }
} // namespace types
PYTHONIC_NS_END

#endif
