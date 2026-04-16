#ifndef PYTHONIC_TYPES_NUMPY_GEXPR_HPP
#define PYTHONIC_TYPES_NUMPY_GEXPR_HPP

#include "pythonic/include/types/numpy_gexpr.hpp"

#include "pythonic/builtins/ValueError.hpp"

#include "pythonic/operator_/iadd.hpp"
#include "pythonic/operator_/idiv.hpp"
#include "pythonic/operator_/imul.hpp"
#include "pythonic/operator_/ior.hpp"
#include "pythonic/operator_/isub.hpp"
#include "pythonic/operator_/ixor.hpp"
#include "pythonic/types/numpy_iexpr.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/meta.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class T>
    struct arange_index;
  }
} // namespace numpy

namespace types
{
  template <class S0, class S1>
  bool slices_may_overlap(S0 const &s0, S1 const &s1)
  {
    if (s0.step >= 0 && s1.step >= 0)
      return s0.lower > s1.lower;
    else
      return s0.lower < s1.lower;
  }
  template <class S>
  bool slices_may_overlap(S const &s, long const &i)
  {
    return s.lower <= i && i < s.upper;
  }
  template <class S>
  bool slices_may_overlap(long const &i, S const &s)
  {
    return s.lower <= i && i < s.upper;
  }

  template <class T>
  constexpr intptr_t baseid_helper(T const &, std::true_type)
  {
    return 0;
  }

  template <class T>
  intptr_t baseid_helper(T const &e, std::false_type)
  {
    return e.baseid();
  }

  template <class T>
  intptr_t baseid(T const &e)
  {
    return baseid_helper(e, std::integral_constant<bool, types::is_dtype<T>::value>{});
  }

  template <class E0, class E1>
  bool may_overlap(E0 const &e0, E1 const &e1)
  {
    return baseid(e0) == baseid(e1);
  }

  template <class E0, class T1>
  bool may_overlap(E0 const &e0, types::list<T1> const &e1)
  {
    return false;
  }

  template <class E0, class T1>
  bool may_overlap(E0 const &e0, types::broadcasted<T1> const &e1)
  {
    return may_overlap(e0, e1.ref);
  }

  template <class E0, class T1, size_t N, class S>
  bool may_overlap(E0 const &e0, types::array_base<T1, N, S> const &e1)
  {
    return false;
  }

  template <class E0, class E1, class E2>
  bool may_overlap(E0 const &e0, types::numpy_vexpr<E1, E2> const &e1)
  {
    return may_overlap(e0, e1.data_);
  }

  template <class E0, class Op, class... Args, size_t... Is>
  bool may_overlap(E0 const &e0, types::numpy_expr<Op, Args...> const &e1,
                   std::index_sequence<Is...>)
  {
    bool overlaps[] = {may_overlap(e0, std::get<Is>(e1.args))...};
    return std::any_of(std::begin(overlaps), std::end(overlaps), [](bool b) { return b; });
  }

  template <class E0, class Op, class... Args>
  bool may_overlap(E0 const &e0, types::numpy_expr<Op, Args...> const &e1)
  {
    return may_overlap(e0, e1, std::make_index_sequence<sizeof...(Args)>());
  }

  template <class E0, class T1>
  bool may_overlap(E0 const &e0, pythonic::numpy::details::arange_index<T1> const &e1)
  {
    return false;
  }

  template <class E0, class T0, class T1>
  bool may_overlap(E0 const &, broadcast<T0, T1> const &)
  {
    return false;
  }
  template <class Arg, class E1, class... S>
  std::enable_if_t<std::is_scalar<E1>::value, bool> may_overlap(numpy_gexpr<Arg, S...> const &gexpr,
                                                                E1 const &)
  {
    return false;
  }

  template <class T, class pS, class Tp, class pSp, class E0, class E1>
  bool may_gexpr_overlap(E0 const &gexpr, E1 const &expr)
  {
    if (!std::is_same<T, Tp>::value) {
      return false;
    }
    if (std::tuple_size<pS>::value != std::tuple_size<pSp>::value) {
      return false;
    }
    if (gexpr.arg.id() != expr.arg.id()) {
      return false;
    }
    if (!slices_may_overlap(std::get<0>(gexpr.slices), std::get<0>(expr.slices)))
      return false;
    return true;
  }
  template <class T, class pS, class Tp, class pSp, class... S, class... Sp>
  bool may_overlap(numpy_gexpr<ndarray<T, pS> const &, S...> const &gexpr,
                   numpy_gexpr<ndarray<Tp, pSp> const &, Sp...> const &expr)
  {
    return may_gexpr_overlap<T, pS, Tp, pSp>(gexpr, expr);
  }
  template <class T, class pS, class Tp, class pSp, class... S, class... Sp>
  bool may_overlap(numpy_gexpr<ndarray<T, pS> &, S...> const &gexpr,
                   numpy_gexpr<ndarray<Tp, pSp> &, Sp...> const &expr)
  {
    return may_gexpr_overlap<T, pS, Tp, pSp>(gexpr, expr);
  }
  template <class T, class pS, class Tp, class pSp, class... S, class... Sp>
  bool may_overlap(numpy_gexpr<ndarray<T, pS> &, S...> const &gexpr,
                   numpy_gexpr<ndarray<Tp, pSp> const &, Sp...> const &expr)
  {
    return may_gexpr_overlap<T, pS, Tp, pSp>(gexpr, expr);
  }
  template <class T, class pS, class Tp, class pSp, class... S, class... Sp>
  bool may_overlap(numpy_gexpr<ndarray<T, pS> const &, S...> const &gexpr,
                   numpy_gexpr<ndarray<Tp, pSp> &, Sp...> const &expr)
  {
    return may_gexpr_overlap<T, pS, Tp, pSp>(gexpr, expr);
  }

  template <class T>
  T to_slice<T>::operator()(T value)
  {
    return value;
  }

  inline fast_contiguous_slice to_slice<none_type>::operator()(none_type)
  {
    return {0, 1};
  }
  template <class T>
  T to_normalized_slice<T>::operator()(T value)
  {
    return value;
  }

  inline cstride_normalized_slice<1> to_normalized_slice<none_type>::operator()(none_type)
  {
    return {0, 1};
  }

  /* helper to build a new shape out of a shape and a slice with new axis
   */
  template <size_t N, class pS, class IsNewAxis>
  auto make_reshape(pS const &shape, IsNewAxis is_new_axis)
      -> decltype(sutils::copy_new_axis<pS::value + N>(shape, is_new_axis))
  {
    return sutils::copy_new_axis<pS::value + N>(shape, is_new_axis);
  }

  /* helper to build an extended slice aka numpy_gexpr out of a subscript
   */

  namespace details
  {

    template <size_t I, class S>
    std::tuple<> merge_gexpr<std::tuple<>, std::tuple<>>::run(S const &, std::tuple<> const &t0,
                                                              std::tuple<> const &)
    {
      return t0;
    }

    template <class... T0>
    template <size_t I, class S>
    std::tuple<T0...> merge_gexpr<std::tuple<T0...>, std::tuple<>>::run(S const &,
                                                                        std::tuple<T0...> const &t0,
                                                                        std::tuple<>)
    {
      return t0;
    }

    template <class T, size_t... Is>
    constexpr long count_new_axis_helper(std::index_sequence<Is...>)
    {
      return count_new_axis<std::tuple_element_t<Is, T>...>::value;
    }

    template <size_t I, class S, class T, size_t... Is>
    auto normalize_all(S const &s, T const &t, std::index_sequence<Is...>)
        -> decltype(std::make_tuple(normalize(
            std::get<Is>(t),
            s.template shape<I + Is -
                             count_new_axis_helper<T>(std::make_index_sequence<1 + Is>())>())...))
    {
      return std::make_tuple(normalize(
          std::get<Is>(t),
          s.template shape<I + Is -
                           count_new_axis_helper<T>(std::make_index_sequence<1 + Is>())>())...);
    }

    template <class... T1>
    template <size_t I, class S>
    std::tuple<normalize_t<T1>...>
    merge_gexpr<std::tuple<>, std::tuple<T1...>>::run(S const &s, std::tuple<>,
                                                      std::tuple<T1...> const &t1)
    {
      return normalize_all<I>(s, t1, std::make_index_sequence<sizeof...(T1)>());
    }

    template <class Arg, class... Sp>
    std::enable_if_t<count_new_axis<Sp...>::value == 0, numpy_gexpr<Arg, Sp...>>
    _make_gexpr(Arg arg, std::tuple<Sp...> const &t)
    {
      return {arg, t};
    }

    template <class Arg, class S, size_t... Is>
    numpy_gexpr<Arg, typename to_normalized_slice<std::tuple_element_t<Is, S>>::type...>
    _make_gexpr_helper(Arg arg, S const &s, std::index_sequence<Is...>)
    {
      return {arg, to_normalized_slice<std::tuple_element_t<Is, S>>{}(std::get<Is>(s))...};
    }

    template <class Arg, class... Sp>
    auto _make_gexpr(Arg arg, std::tuple<Sp...> const &s) -> std::enable_if_t<
        count_new_axis<Sp...>::value != 0,
        decltype(_make_gexpr_helper(
            arg.reshape(make_reshape<count_new_axis<Sp...>::value>(
                arg, std::tuple<std::integral_constant<bool, to_slice<Sp>::is_new_axis>...>())),
            s, std::make_index_sequence<sizeof...(Sp)>()))>
    {
      return _make_gexpr_helper(
          arg.reshape(make_reshape<count_new_axis<Sp...>::value>(
              arg, std::tuple<std::integral_constant<bool, to_slice<Sp>::is_new_axis>...>())),
          s, std::make_index_sequence<sizeof...(Sp)>());
    }

    template <class Arg, class... S>
    template <size_t... Is>
    numpy_gexpr<Arg, normalize_t<S>...>
    make_gexpr<Arg, S...>::operator()(Arg arg, std::tuple<S...> s, std::index_sequence<Is...>)
    {
      return {arg, normalize(std::get<Is>(s), arg.template shape<Is>())...};
    }

    template <class Arg, class... S>
    numpy_gexpr<Arg, normalize_t<S>...> make_gexpr<Arg, S...>::operator()(Arg arg, S const &...s)
    {
      return operator()(arg, std::tuple<S...>(s...), std::make_index_sequence<sizeof...(S)>());
    }
  } // namespace details

  template <class Arg, class... S>
  auto make_gexpr(Arg &&arg, S const &...s)
      -> decltype(details::make_gexpr<Arg, S...>{}(std::forward<Arg>(arg), s...))
  {
    return details::make_gexpr<Arg, S...>{}(std::forward<Arg>(arg), s...);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...>::numpy_gexpr() : buffer(nullptr)
  {
  }

  template <class Arg, class... S>
  template <class Argp> // not using the default one, to make it possible to
  // accept reference && non reference version of Argp
  numpy_gexpr<Arg, S...>::numpy_gexpr(numpy_gexpr<Argp, S...> const &other)
      : arg(other.arg), slices(other.slices), _shape(other._shape), buffer(other.buffer),
        _strides(other._strides)
  {
    static_assert(
        std::is_same<typename returnable<Arg>::type, typename returnable<Argp>::type>::value,
        "this constructor is only here to adapt reference / non "
        "reference type, nothing else");
    assert(buffer);
  }

  template <class Arg, class... S>
  template <size_t J, class Slice>
  std::enable_if_t<is_normalized_slice<Slice>::value, void>
  numpy_gexpr<Arg, S...>::init_shape(Slice const &s, utils::int_<1>, utils::int_<J>)
  {
    buffer += s.lower * arg.template strides<sizeof...(S) - 1>();
    sutils::assign(std::get<J>(_strides), s.step * arg.template strides<sizeof...(S) - 1>());
    sutils::assign(std::get<J>(_shape), std::get<sizeof...(S) - 1>(slices).size());
  }

  template <class Arg, class... S>
  template <size_t I, size_t J, class Slice>
  std::enable_if_t<is_normalized_slice<Slice>::value, void>
  numpy_gexpr<Arg, S...>::init_shape(Slice const &s, utils::int_<I>, utils::int_<J>)
  {
    sutils::assign(std::get<J>(_shape), std::get<sizeof...(S) - I>(slices).size());
    buffer += s.lower * arg.template strides<sizeof...(S) - I>();
    sutils::assign(std::get<J>(_strides), s.step * arg.template strides<sizeof...(S) - I>());
    init_shape(std::get<sizeof...(S) - I + 1>(slices), utils::int_<I - 1>(), utils::int_<J + 1>());
  }

  template <class Arg, class... S>
  template <size_t J>
  void numpy_gexpr<Arg, S...>::init_shape(long cs, utils::int_<1>, utils::int_<J>)
  {
    assert(cs >= 0 && "normalized");
    buffer += cs * arg.template strides<sizeof...(S) - 1>();
  }

  template <class Arg, class... S>
  template <size_t I, size_t J>
  void numpy_gexpr<Arg, S...>::init_shape(long cs, utils::int_<I>, utils::int_<J>)
  {
    assert(cs >= 0 && "normalized");
    buffer += cs * arg.template strides<sizeof...(S) - I>();
    init_shape(std::get<sizeof...(S) - I + 1>(slices), utils::int_<I - 1>(), utils::int_<J>());
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...>::numpy_gexpr(Arg const &arg, std::tuple<S const &...> const &values)
      : arg(arg), slices(values), buffer(const_cast<dtype *>(this->arg.buffer))
  {
    assert(buffer);
    init_shape(std::get<0>(slices), utils::int_<sizeof...(S)>(), utils::int_<0>());

    sutils::copy_shape<sizeof...(S) - count_long<S...>::value, count_long<S...>::value>(
        _shape, arg, std::make_index_sequence<value - (sizeof...(S) - count_long<S...>::value)>());

    sutils::copy_strides<sizeof...(S) - count_long<S...>::value, count_long<S...>::value>(
        _strides, arg,
        std::make_index_sequence<value - (sizeof...(S) - count_long<S...>::value)>());
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...>::numpy_gexpr(Arg const &arg, S const &...s)
      : numpy_gexpr(arg, std::tuple<S const &...>(s...))
  {
  }
  template <class Arg, class... S>
  template <class Argp, class... Sp>
  numpy_gexpr<Arg, S...>::numpy_gexpr(numpy_gexpr<Argp, Sp...> const &expr, Arg arg)
      : arg(arg), slices(tuple_pop(expr.slices)), buffer(expr.buffer)
  {
    assert(buffer);
    sutils::copy_shape<0, 1>(_shape, expr, std::make_index_sequence<value>());
    buffer += arg.buffer - expr.arg.buffer;
    sutils::copy_strides<0, 1>(_strides, expr, std::make_index_sequence<value>());
  }

  template <class Arg, class... S>
  template <class G>
  numpy_gexpr<Arg, S...>::numpy_gexpr(G const &expr, Arg &&arg)
      : arg(std::forward<Arg>(arg)), slices(tuple_pop(expr.slices)), buffer(expr.buffer)
  {
    assert(buffer);
    sutils::copy_shape<0, 1>(_shape, expr, std::make_index_sequence<value>());
    buffer += (arg.buffer - expr.arg.buffer);
    sutils::copy_strides<0, 1>(_strides, expr, std::make_index_sequence<value>());
  }

  template <class Arg, class... S>
  template <class E>
  std::enable_if_t<may_overlap_gexpr<E>::value, numpy_gexpr<Arg, S...> &>
  numpy_gexpr<Arg, S...>::_copy(E const &expr)
  {
    static_assert(value >= utils::dim_of<E>::value, "dimensions match");
    /* at this point, we could not statically check that there is not an
     * aliasing issue that would require an extra copy because of the vector
     * assignment
     * perform a fuzzy alias check dynamically!
     */
    assert(buffer);
    constexpr bool vectorize = is_vectorizable &&
                               std::is_same<dtype, typename dtype_of<E>::type>::value &&
                               is_vectorizable_array<E>::value;
    if (may_overlap(*this, expr)) {
      return utils::broadcast_copy<numpy_gexpr &, ndarray<typename E::dtype, typename E::shape_t>,
                                   value, value - utils::dim_of<E>::value, vectorize>(
          *this, ndarray<typename E::dtype, typename E::shape_t>(expr));
    } else {
      // 100% sure there's no overlap
      return utils::broadcast_copy<numpy_gexpr &, E, value, value - utils::dim_of<E>::value,
                                   vectorize>(*this, expr);
    }
  }

  template <class Arg, class... S>
  template <class E>
  std::enable_if_t<!may_overlap_gexpr<E>::value, numpy_gexpr<Arg, S...> &>
  numpy_gexpr<Arg, S...>::_copy(E const &expr)
  {
    return _copy_restrict(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::_copy_restrict(E const &expr)
  {
    constexpr bool vectorize = is_vectorizable &&
                               std::is_same<dtype, typename dtype_of<E>::type>::value &&
                               is_vectorizable_array<E>::value;
    assert(buffer);
    return utils::broadcast_copy<numpy_gexpr &, E, value, (int)value - (int)utils::dim_of<E>::value,
                                 vectorize>(*this, expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator=(E const &expr)
  {
    return _copy(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator=(numpy_gexpr<Arg, S...> const &expr)
  {
    if (buffer == nullptr) {
      // arg = expr.arg;
      const_cast<std::decay_t<Arg> &>(arg) = expr.arg;
      slices = expr.slices;
      assert(expr.buffer);
      buffer = arg.buffer + (expr.buffer - expr.arg.buffer);
      _shape = expr._shape;
      _strides = expr._strides;
      assert(sutils::getshape(*this) == sutils::getshape(expr) && "compatible sizes");
      return *this;
    } else {
      return _copy(expr);
    }
  }

  template <class Arg, class... S>
  template <class Argp>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator=(numpy_gexpr<Argp, S...> const &expr)
  {
    if (buffer == nullptr) {
      // arg = expr.arg;
      const_cast<std::decay_t<Arg> &>(arg) = expr.arg;
      slices = expr.slices;
      assert(expr.buffer);
      buffer = arg.buffer + (expr.buffer - expr.arg.buffer);
      _shape = expr._shape;
      _strides = expr._strides;
      return *this;
    } else {
      return _copy(expr);
    }
  }

  template <class Arg, class... S>
  template <class Op, class E>
  std::enable_if_t<!may_overlap_gexpr<E>::value, numpy_gexpr<Arg, S...> &>
  numpy_gexpr<Arg, S...>::update_(E const &expr)
  {
    using BExpr = std::conditional_t<std::is_scalar<E>::value, broadcast<E, dtype>, E const &>;
    BExpr bexpr = expr;
    // 100% sure there's no overlap
    return utils::broadcast_update < Op, numpy_gexpr &, BExpr, value,
           value - (std::is_scalar<E>::value + utils::dim_of<E>::value),
           is_vectorizable &&
               types::is_vectorizable<std::remove_cv_t<std::remove_reference_t<BExpr>>>::value &&
               std::is_same<dtype, typename dtype_of<std::decay_t<BExpr>>::type>::value >
                   (*this, bexpr);
  }

  template <class Arg, class... S>
  template <class Op, class E>
  std::enable_if_t<may_overlap_gexpr<E>::value, numpy_gexpr<Arg, S...> &>
  numpy_gexpr<Arg, S...>::update_(E const &expr)
  {
    using BExpr = std::conditional_t<std::is_scalar<E>::value, broadcast<E, dtype>, E const &>;
    BExpr bexpr = expr;

    if (may_overlap(*this, expr)) {
      using NBExpr = ndarray<typename std::remove_reference_t<BExpr>::dtype,
                             typename std::remove_reference_t<BExpr>::shape_t>;
      return utils::broadcast_update < Op, numpy_gexpr &, NBExpr, value,
             value - (std::is_scalar<E>::value + utils::dim_of<E>::value),
             is_vectorizable && types::is_vectorizable<E>::value &&
                 std::is_same<dtype, typename std::decay_t<BExpr>::dtype>::value >
                     (*this, NBExpr(bexpr));
    } else {
      // 100% sure there's no overlap
      return utils::broadcast_update < Op, numpy_gexpr &, BExpr, value,
             value - (std::is_scalar<E>::value + utils::dim_of<E>::value),
             is_vectorizable && types::is_vectorizable<E>::value &&
                 std::is_same<dtype, typename std::decay_t<BExpr>::dtype>::value > (*this, bexpr);
    }
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator+=(E const &expr)
  {
    return update_<pythonic::operator_::functor::iadd>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator+=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::iadd>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator-=(E const &expr)
  {
    return update_<pythonic::operator_::functor::isub>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator-=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::isub>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator*=(E const &expr)
  {
    return update_<pythonic::operator_::functor::imul>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator*=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::imul>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator/=(E const &expr)
  {
    return update_<pythonic::operator_::functor::idiv>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator/=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::idiv>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator|=(E const &expr)
  {
    return update_<pythonic::operator_::functor::ior>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator|=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::ior>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator&=(E const &expr)
  {
    return update_<pythonic::operator_::functor::iand>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator&=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::iand>(expr);
  }

  template <class Arg, class... S>
  template <class E>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator^=(E const &expr)
  {
    return update_<pythonic::operator_::functor::ixor>(expr);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...> &numpy_gexpr<Arg, S...>::operator^=(numpy_gexpr<Arg, S...> const &expr)
  {
    return update_<pythonic::operator_::functor::ixor>(expr);
  }

  template <class Arg, class... S>
  typename numpy_gexpr<Arg, S...>::const_iterator numpy_gexpr<Arg, S...>::begin() const
  {
    return make_const_nditerator<is_strided || value != 1>()(*this, 0);
  }

  template <class Arg, class... S>
  typename numpy_gexpr<Arg, S...>::const_iterator numpy_gexpr<Arg, S...>::end() const
  {
    return make_const_nditerator<is_strided || value != 1>()(*this, size());
  }

  template <class Arg, class... S>
  typename numpy_gexpr<Arg, S...>::iterator numpy_gexpr<Arg, S...>::begin()
  {
    return make_nditerator<is_strided || value != 1>()(*this, 0);
  }

  template <class Arg, class... S>
  typename numpy_gexpr<Arg, S...>::iterator numpy_gexpr<Arg, S...>::end()
  {
    return make_nditerator<is_strided || value != 1>()(*this, size());
  }

#ifdef USE_XSIMD
  template <class Arg, class... S>
  template <class vectorizer>
  typename numpy_gexpr<Arg, S...>::simd_iterator numpy_gexpr<Arg, S...>::vbegin(vectorizer) const
  {
    return {buffer};
  }

  template <class Arg, class... S>
  template <class vectorizer>
  typename numpy_gexpr<Arg, S...>::simd_iterator numpy_gexpr<Arg, S...>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {buffer + long(size() / vector_size * vector_size)};
  }

#endif

  template <class Arg, class... S>
  auto numpy_gexpr<Arg, S...>::operator[](long i) const -> decltype(this->fast(i))
  {
    if (i < 0)
      i += std::get<0>(_shape);
    return fast(i);
  }

  template <class Arg, class... S>
  auto numpy_gexpr<Arg, S...>::operator[](long i) -> decltype(this->fast(i))
  {
    if (i < 0)
      i += std::get<0>(_shape);
    return fast(i);
  }

  template <class Arg, class... S>
  template <class... Sp>
  auto numpy_gexpr<Arg, S...>::operator()(Sp const &...s) const -> decltype(make_gexpr(*this, s...))
  {
    return make_gexpr(*this, s...);
  }

  template <class Arg, class... S>
  template <class Sp>
  auto numpy_gexpr<Arg, S...>::operator[](Sp const &s) const
      -> std::enable_if_t<is_slice<Sp>::value, decltype(make_gexpr(*this, (s.lower, s)))>
  {
    return make_gexpr(*this, s);
  }

  template <class Arg, class... S>
  template <size_t M>
  auto numpy_gexpr<Arg, S...>::fast(
      array_tuple<long, M> const &indices) const & -> decltype(nget<M - 1>().fast(*this, indices))
  {
    return nget<M - 1>().fast(*this, indices);
  }

  template <class Arg, class... S>
  template <size_t M>
  auto numpy_gexpr<Arg, S...>::fast(array_tuple<long, M> const &indices)
      && -> decltype(nget<M - 1>().fast(std::move(*this), indices))
  {
    return nget<M - 1>().fast(std::move(*this), indices);
  }

  template <class Arg, class... S>
  template <size_t M>
  auto numpy_gexpr<Arg, S...>::operator[](
      array_tuple<long, M> const &indices) const & -> decltype(nget<M - 1>()(*this, indices))
  {
    return nget<M - 1>()(*this, indices);
  }

  template <class Arg, class... S>
  template <size_t M>
  auto numpy_gexpr<Arg, S...>::operator[](
      array_tuple<long, M> const &indices) && -> decltype(nget<M - 1>()(std::move(*this), indices))
  {
    return nget<M - 1>()(std::move(*this), indices);
  }

  template <class Arg, class... S>
  template <class F>
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value,
                   numpy_vexpr<numpy_gexpr<Arg, S...>, ndarray<long, pshape<long>>>>
  numpy_gexpr<Arg, S...>::fast(F const &filter) const
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

  template <class Arg, class... S>
  template <class F>
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value,
                   numpy_vexpr<numpy_gexpr<Arg, S...>, ndarray<long, pshape<long>>>>
  numpy_gexpr<Arg, S...>::operator[](F const &filter) const
  {
    return fast(filter);
  }

  template <class Arg, class... S>
  numpy_gexpr<Arg, S...>::operator bool() const
  {
    if (sutils::any_of(*this, [](long n) { return n != 1; }))
      throw ValueError("The truth value of an array with more than one element "
                       "is ambiguous. Use a.any() or a.all()");
    return *buffer;
  }

  template <class Arg, class... S>
  long numpy_gexpr<Arg, S...>::flat_size() const
  {
    return sutils::prod(*this);
  }

  template <class Arg, class... S>
  long numpy_gexpr<Arg, S...>::size() const
  {
    return std::get<0>(_shape);
  }
} // namespace types
PYTHONIC_NS_END

#endif
