#ifndef PYTHONIC_NUMPY_ARGMINMAX_HPP
#define PYTHONIC_NUMPY_ARGMINMAX_HPP

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class P, size_t... Is>
    P iota(std::index_sequence<Is...>)
    {
      return {static_cast<typename P::value_type>(Is)...};
    }

    template <class P>
    P iota()
    {
      return iota<P>(std::make_index_sequence<P::size>());
    }
  } // namespace details
  template <class Op, class E, class T>
  long _argminmax_seq(E const &elts, T &minmax_elts)
  {
    long index = 0;
    long res = -1;
    for (auto const &elt : elts) {
      if (Op::value(elt, minmax_elts)) {
        minmax_elts = elt;
        res = index;
      }
      ++index;
    }

    return res;
  }
  template <class Op, class E, class T>
#ifdef USE_XSIMD
  std::enable_if_t<!E::is_vectorizable || !types::is_vector_op<typename Op::op, T, T>::value ||
                       std::is_same<typename E::dtype, bool>::value,
                   long>
#else
  long
#endif
  _argminmax(E const &elts, T &minmax_elts, utils::int_<1>)
  {
    return _argminmax_seq<Op>(elts, minmax_elts);
  }

  template <class Op, class E, class T, class... Indices>
  std::tuple<long, long> _argminmax_fast(E const &elts, T &minmax_elts, long current_pos,
                                         utils::int_<1>, Indices... indices)
  {
    long res = -1;
    long n = elts.template shape<std::decay_t<E>::value - 1>();
    for (long i = 0; i < n; ++i) {
      auto elt = elts.load(indices..., i);
      if (Op::value(elt, minmax_elts)) {
        minmax_elts = elt;
        res = current_pos + i;
      }
    }

    return std::make_tuple(res, current_pos + n);
  }

#ifdef USE_XSIMD
  template <class Op, class E, class T>
  std::enable_if_t<E::is_vectorizable && types::is_vector_op<typename Op::op, T, T>::value &&
                       !std::is_same<typename E::dtype, bool>::value,
                   long>
  _argminmax(E const &elts, T &minmax_elts, utils::int_<1>)
  {
    using vT = xsimd::batch<T>;
    using iT = xsimd::as_integer_t<T>;
    static const size_t vN = vT::size;
    const long n = elts.size();
    if (n >= std::numeric_limits<iT>::max()) {
      return _argminmax_seq<Op>(elts, minmax_elts);
    }

    auto viter = types::vectorizer_nobroadcast::vbegin(elts),
         vend = types::vectorizer_nobroadcast::vend(elts);

    const long bound = std::distance(viter, vend);
    long minmax_index = -1;
    if (bound > 0) {
      auto vacc = *viter;
      alignas(sizeof(vT)) iT iota[vN] = {0};
      for (long i = 0; i < (long)vN; ++i)
        iota[i] = i;
      xsimd::batch<iT> curr = xsimd::load_aligned(iota);
      xsimd::batch<iT> indices = curr;
      xsimd::batch<iT> step(vN);

      for (++viter, curr += step; viter != vend; ++viter, curr += step) {
        auto c = *viter;
        // In order to keep the first element that matches the condition,
        // we compare to the previous vacc and select indices based on that
        // mask.
        auto prev_vacc = vacc;
        vacc = typename Op::op{}(vacc, c);
        auto mask = xsimd::batch_bool_cast<iT>(prev_vacc != vacc);
        indices = xsimd::select(mask, curr, indices);
      }

      alignas(sizeof(vT)) T stored[vN];
      vacc.store_aligned(&stored[0]);
      alignas(sizeof(vT)) long indexed[vN];
      indices.store_aligned(&indexed[0]);

      for (size_t j = 0; j < vN; ++j) {
        if (Op::value(stored[j], minmax_elts)) {
          minmax_elts = stored[j];
          minmax_index = indexed[j];
        }
      }
    }
    auto iter = elts.begin() + bound * vN;

    for (long i = bound * vN; i < n; ++i, ++iter) {
      if (Op::value(*iter, minmax_elts)) {
        minmax_elts = *iter;
        minmax_index = i;
      }
    }
    return minmax_index;
  }
#endif

  template <class Op, class E, size_t N, class T>
  long _argminmax(E const &elts, T &minmax_elts, utils::int_<N>)
  {
    long current_pos = 0;
    long current_minmaxarg = 0;
    for (auto &&elt : elts) {
      long v = _argminmax<Op>(elt, minmax_elts, utils::int_<N - 1>());
      if (v >= 0)
        current_minmaxarg = current_pos + v;
      current_pos += elt.flat_size();
    }
    return current_minmaxarg;
  }
  template <class Op, class E, size_t N, class T, class... Indices>
  std::enable_if_t<N != 1, std::tuple<long, long>> _argminmax_fast(E const &elts, T &minmax_elts,
                                                                   long current_pos, utils::int_<N>,
                                                                   Indices... indices)
  {
    long current_minmaxarg = 0;
    for (long i = 0, n = elts.template shape<std::decay_t<E>::value - N>(); i < n; ++i) {
      long v;
      std::tie(v, current_pos) =
          _argminmax_fast<Op>(elts, minmax_elts, current_pos, utils::int_<N - 1>(), indices..., i);
      if (v >= 0)
        current_minmaxarg = v;
    }
    return std::make_tuple(current_minmaxarg, current_pos);
  }

  template <class Op, class E>
  long argminmax(E const &expr)
  {
    if (!expr.flat_size())
      throw types::ValueError("empty sequence");
    using elt_type = typename E::dtype;
    elt_type argminmax_value = Op::limit();
#ifndef USE_XSIMD
    if (utils::no_broadcast_ex(expr)) {
      return std::get<0>(_argminmax_fast<Op>(expr, argminmax_value, 0, utils::int_<E::value>()));
    } else
#endif
      return _argminmax<Op>(expr, argminmax_value, utils::int_<E::value>());
  }

  template <class Op, size_t Dim, size_t Axis, class T, class E, class V>
  void _argminmax_tail(T &&out, E const &expr, long curr, V &&curr_minmax,
                       std::integral_constant<size_t, 0>)
  {
    if (Op::value(expr, curr_minmax)) {
      out = curr;
      curr_minmax = expr;
    }
  }

  template <class Op, size_t Dim, size_t Axis, class T, class E, class V, size_t N>
  std::enable_if_t<Axis != (Dim - N), void> _argminmax_tail(T &&out, E const &expr, long curr,
                                                            V &&curr_minmax,
                                                            std::integral_constant<size_t, N>)
  {
    static_assert(N >= 1, "specialization ok");
    long i = 0;
    for (auto &&elt : expr) {
      _argminmax_tail<Op, Dim, Axis>(out.fast(i), elt, curr, curr_minmax.fast(i),
                                     std::integral_constant<size_t, N - 1>());
      ++i;
    }
  }

  template <class Op, size_t Dim, size_t Axis, class T, class E>
  std::enable_if_t<Axis == (Dim - 1), void> _argminmax_head(T &&out, E const &expr,
                                                            std::integral_constant<size_t, 1>)
  {
    typename E::dtype val = Op::limit();
    long i = 0;
    for (auto &&elt : expr)
      _argminmax_tail<Op, Dim, Axis>(out, elt, i++, val, std::integral_constant<size_t, 0>());
  }

  template <class Op, size_t Dim, size_t Axis, class T, class E, size_t N>
  std::enable_if_t<Axis == (Dim - N), void> _argminmax_head(T &&out, E const &expr,
                                                            std::integral_constant<size_t, N>)
  {
    static_assert(N > 1, "specialization ok");
    types::ndarray<typename E::dtype, types::array_tuple<long, N - 1>> val{sutils::getshape(out),
                                                                           Op::limit()};
    long i = 0;
    for (auto &&elt : expr) {
      _argminmax_tail<Op, Dim, Axis>(out, elt, i++, val, std::integral_constant<size_t, N - 1>());
    }
  }

  template <class Op, size_t Dim, size_t Axis, class T, class E, size_t N>
  std::enable_if_t<Axis != (Dim - N), void> _argminmax_head(T &&out, E const &expr,
                                                            std::integral_constant<size_t, N>)
  {
    static_assert(N >= 1, "specialization ok");
    auto out_iter = out.begin();
    for (auto &&elt : expr) {
      _argminmax_head<Op, Dim, Axis>(*out_iter, elt, std::integral_constant<size_t, N - 1>());
      ++out_iter;
    }
  }

  template <class Op, size_t N, class T, class E, size_t... Axis>
  void _argminmax_pick_axis(long axis, T &&out, E const &expr, std::index_sequence<Axis...>)
  {
    (void)std::initializer_list<bool>{
        ((Axis == axis) &&
         (_argminmax_head<Op, N, Axis>(out, expr, std::integral_constant<size_t, N>()), true))...};
  }

  template <class Op, class E>
  types::ndarray<long, types::array_tuple<long, E::value - 1>> argminmax(E const &array, long axis)
  {
    if (axis < 0)
      axis += E::value;
    if (axis < 0 || size_t(axis) >= E::value)
      throw types::ValueError("axis out of bounds");
    auto shape = sutils::getshape(array);
    types::array_tuple<long, E::value - 1> shp;
    auto next = std::copy(shape.begin(), shape.begin() + axis, shp.begin());
    std::copy(shape.begin() + axis + 1, shape.end(), next);
    types::ndarray<long, types::array_tuple<long, E::value - 1>> out{shp, builtins::None};
    _argminmax_pick_axis<Op, E::value>(axis, out, array, std::make_index_sequence<E::value>());
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
