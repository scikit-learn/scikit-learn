#ifndef PYTHONIC_NUMPY_REDUCE_HPP
#define PYTHONIC_NUMPY_REDUCE_HPP

#include "pythonic/include/numpy/reduce.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/neutral.hpp"

#ifdef USE_XSIMD
#include <xsimd/xsimd.hpp>
#endif

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class Op, size_t N, class vector_form>
  struct _reduce {
    template <class E, class F>
    F operator()(E &&e, F acc)
    {
      for (auto &&value : std::forward<E>(e))
        acc = _reduce<Op, N - 1, vector_form>{}(std::forward<decltype(value)>(value), acc);

      return acc;
    }
  };

  template <class Op, class vector_form>
  struct _reduce<Op, 1, vector_form> {
    template <class E, class F>
    F operator()(E &&e, F acc)
    {
      for (auto value : std::forward<E>(e)) {
        Op{}(acc, value);
      }

      return acc;
    }
  };

  template <class Op, size_t N>
  struct _reduce<Op, N, types::novectorize_nobroadcast> {
    template <class E, class F, class... Indices>
    F operator()(E &&e, F acc, Indices... indices)
    {
      for (long i = 0, n = e.template shape<std::decay_t<E>::value - N>(); i < n; ++i) {
        acc = _reduce<Op, N - 1, types::novectorize_nobroadcast>{}(e, acc, indices..., i);
      }
      return acc;
    }
  };

  template <class Op>
  struct _reduce<Op, 1, types::novectorize_nobroadcast> {
    template <class E, class F, class... Indices>
    F operator()(E &&e, F acc, Indices... indices)
    {
      for (long i = 0, n = e.template shape<std::decay_t<E>::value - 1>(); i < n; ++i) {
        Op{}(acc, e.load(indices..., i));
      }
      return acc;
    }
  };

#ifdef USE_XSIMD
  template <class vectorizer, class Op, class E, class F>
  F vreduce(E e, F acc)
  {
    using T = typename E::dtype;
    using vT = xsimd::batch<T>;
    static const size_t vN = vT::size;
    const long n = e.size();
    auto viter = vectorizer::vbegin(e), vend = vectorizer::vend(e);
    const long bound = std::distance(viter, vend);
    if (bound > 0) {
      auto vacc = *viter;
      for (++viter; viter != vend; ++viter)
        Op{}(vacc, *viter);
      alignas(sizeof(vT)) T stored[vN];
      vacc.store_aligned(&stored[0]);
      for (size_t j = 0; j < vN; ++j)
        Op{}(acc, stored[j]);
    }
    auto iter = e.begin() + bound * vN;

    for (long i = bound * vN; i < n; ++i, ++iter) {
      Op{}(acc, *iter);
    }
    return acc;
  }

  template <class Op>
  struct _reduce<Op, 1, types::vectorizer> {
    template <class E, class F>
    F operator()(E &&e, F acc)
    {
      return vreduce<types::vectorizer, Op>(std::forward<E>(e), acc);
    }
  };
  template <class Op>
  struct _reduce<Op, 1, types::vectorizer_nobroadcast> {
    template <class E, class F>
    F operator()(E &&e, F acc)
    {
      return vreduce<types::vectorizer_nobroadcast, Op>(std::forward<E>(e), acc);
    }
  };
#else
  template <class Op, size_t N>
  struct _reduce<Op, N, types::vectorizer_nobroadcast>
      : _reduce<Op, N, types::novectorize_nobroadcast> {
  };
  template <class Op>
  struct _reduce<Op, 1, types::vectorizer_nobroadcast>
      : _reduce<Op, 1, types::novectorize_nobroadcast> {
  };
#endif
  template <class Op, class E, bool vector_form>
  struct reduce_helper;

  template <class Op, class E>
  struct reduce_helper<Op, E, false> {
    template <class T>
    reduce_result_type<Op, E> operator()(E const &expr, T p) const
    {
      if (utils::no_broadcast_ex(expr))
        return _reduce<Op, E::value, types::novectorize_nobroadcast>{}(expr, p);
      else
        return _reduce<Op, E::value, types::novectorize>{}(expr, p);
    }
  };
  template <class Op, class E>
  struct reduce_helper<Op, E, true> {
    template <class T>
    reduce_result_type<Op, E> operator()(E const &expr, T p) const
    {
      if (utils::no_broadcast_vectorize(expr))
        return _reduce<Op, E::value, types::vectorizer_nobroadcast>{}(expr, p);
      else
        return _reduce<Op, E::value, types::vectorizer>{}(expr, p);
    }
  };
  template <class Op, class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, E>
  reduce(E const &expr, types::none_type)
  {
    return expr;
  }

  template <class Op, class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, E>
  reduce(E const &array, long axis)
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return reduce<Op>(array);
  }

  template <class Op, class E, class dtype>
  std::enable_if_t<types::is_numexpr_arg<E>::value, reduce_result_type<Op, E, dtype>>
  reduce(E const &expr, types::none_type axis, dtype)
  {
    using rrt = reduce_result_type<Op, E, dtype>;
    bool constexpr is_vectorizable = E::is_vectorizable &&
                                     !std::is_same<typename E::dtype, bool>::value &&
                                     std::is_same<rrt, typename E::dtype>::value;
    rrt p = utils::neutral<Op, rrt>::value;
    return reduce_helper<Op, E, is_vectorizable>{}(expr, p);
  }

  template <class Op, class E, class dtype>
  std::enable_if_t<E::value == 1, reduce_result_type<Op, E, dtype>>
  reduce(E const &array, long axis, dtype d, types::none_type)
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return reduce<Op>(array, types::none_type{}, d);
  }

  template <class Op, class E, class Out>
  std::enable_if_t<E::value == 1, reduce_result_type<Op, E>> reduce(E const &array, long axis,
                                                                    types::none_type, Out &&out)
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return std::forward<Out>(out) = reduce<Op>(array);
  }

  template <class Op, size_t N>
  struct _reduce_axisb {
    template <class E, class F, class EIndices, class FIndices>
    void operator()(E &&e, F &&f, long axis, EIndices &&e_indices, FIndices &&f_indices)
    {
      for (long i = 0, n = e.template shape<std::decay_t<E>::value - N>(); i < n; ++i) {
        _reduce_axisb<Op, N - 1>{}(e, f, axis, std::tuple_cat(e_indices, std::make_tuple(i)),
                                   std::tuple_cat(f_indices, std::make_tuple(i)));
      }
    }
  };

  template <class Op>
  struct _reduce_axisb<Op, 0> {
    template <class E, class F, class EIndices, class FIndices, size_t... Es, size_t... Fs>
    void helper(E &&e, F &&f, EIndices &&e_indices, FIndices &&f_indices,
                std::index_sequence<Es...>, std::index_sequence<Fs...>)
    {
      f.template update<Op>(e.load(std::get<Es>(e_indices)...), (long)std::get<Fs>(f_indices)...);
    }
    template <class E, class F, class EIndices, class FIndices>
    void operator()(E &&e, F &&f, long axis, EIndices &&e_indices, FIndices &&f_indices)
    {
      helper(std::forward<E>(e), std::forward<F>(f), e_indices, f_indices,
             std::make_index_sequence<std::tuple_size<std::decay_t<EIndices>>::value>(),
             std::make_index_sequence<std::tuple_size<std::decay_t<FIndices>>::value>());
    }
  };

  template <class Op, size_t N>
  struct _reduce_axis {
    template <class E, class F, class EIndices, class FIndices>
    void operator()(E &&e, F &&f, long axis, EIndices &&e_indices, FIndices &&f_indices)
    {
      if (axis == std::decay_t<E>::value - N) {
        for (long i = 0, n = e.template shape<std::decay_t<E>::value - N>(); i < n; ++i) {
          _reduce_axisb<Op, N - 1>{}(e, f, axis, std::tuple_cat(e_indices, std::make_tuple(i)),
                                     std::forward<FIndices>(f_indices));
        }
      } else {
        for (long i = 0, n = e.template shape<std::decay_t<E>::value - N>(); i < n; ++i) {
          _reduce_axis<Op, N - 1>{}(e, f, axis, std::tuple_cat(e_indices, std::make_tuple(i)),
                                    std::tuple_cat(f_indices, std::make_tuple(i)));
        }
      }
    }
  };
  template <class Op>
  struct _reduce_axis<Op, 0> {
    template <class E, class F, class EIndices, class FIndices>
    void operator()(E &&e, F &&f, long axis, EIndices &&e_indices, FIndices &&f_indices)
    {
    }
  };

  template <class Op, class E, class dtype>
  std::enable_if_t<E::value != 1, reduced_type<E, Op, dtype>> reduce(E const &array, long axis,
                                                                     dtype, types::none_type)
  {
    if (axis < 0)
      axis += E::value;
    if (axis < 0 || size_t(axis) >= E::value)
      throw types::ValueError("axis out of bounds");
    types::array_tuple<long, E::value - 1> shp;
    auto tmp = sutils::getshape(array);
    auto next = std::copy(tmp.begin(), tmp.begin() + axis, shp.begin());
    std::copy(tmp.begin() + axis + 1, tmp.end(), next);
    reduced_type<E, Op, dtype> out{shp, builtins::None};
    std::fill(out.begin(), out.end(), utils::neutral<Op, typename E::dtype>::value);
    return reduce<Op>(array, axis, types::none_type{}, out);
  }
  template <class Op, class E, class Out>
  std::enable_if_t<E::value != 1, reduced_type<E, Op>> reduce(E const &array, long axis,
                                                              types::none_type, Out &&out)
  {
    if (axis < 0)
      axis += E::value;
    if (axis < 0 || size_t(axis) >= E::value)
      throw types::ValueError("axis out of bounds");
    if (utils::no_broadcast(array)) {
      std::fill(out.begin(), out.end(), utils::neutral<Op, typename E::dtype>::value);
      _reduce_axis<Op, E::value>{}(array, std::forward<Out>(out), axis, std::make_tuple(),
                                   std::make_tuple());
      return std::forward<Out>(out);
    } else {
      if (axis == 0) {
        std::fill(out.begin(), out.end(), utils::neutral<Op, typename E::dtype>::value);
        return _reduce<Op, 1, types::novectorize /* not on scalars*/>{}(array,
                                                                        std::forward<Out>(out));
      } else {
        std::transform(array.begin(), array.end(), out.begin(),
                       [axis](typename E::const_iterator::value_type other) {
                         return reduce<Op>(other, axis - 1);
                       });
        return std::forward<Out>(out);
      }
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
