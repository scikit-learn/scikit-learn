#ifndef PYTHONIC_TYPES_NUMPY_TEXPR_HPP
#define PYTHONIC_TYPES_NUMPY_TEXPR_HPP

#include "pythonic/include/types/numpy_texpr.hpp"

#include "pythonic/numpy/array.hpp"
#include "pythonic/numpy/transpose.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/allocate.hpp"

#include "pythonic/operator_/iadd.hpp"
#include "pythonic/operator_/iand.hpp"
#include "pythonic/operator_/idiv.hpp"
#include "pythonic/operator_/imul.hpp"
#include "pythonic/operator_/ior.hpp"
#include "pythonic/operator_/isub.hpp"
#include "pythonic/operator_/ixor.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  template <class E>
  numpy_texpr_2<E>::numpy_texpr_2()
  {
  }

  template <class E>
  numpy_texpr_2<E>::numpy_texpr_2(Arg const &arg) : arg(arg)
  {
  }

  template <class E>
  typename numpy_texpr_2<E>::const_iterator numpy_texpr_2<E>::begin() const
  {
    return {*this, 0};
  }

  template <class E>
  typename numpy_texpr_2<E>::const_iterator numpy_texpr_2<E>::end() const
  {
    return {*this, size()};
  }

  template <class E>
  typename numpy_texpr_2<E>::iterator numpy_texpr_2<E>::begin()
  {
    return {*this, 0};
  }

  template <class E>
  typename numpy_texpr_2<E>::iterator numpy_texpr_2<E>::end()
  {
    return {*this, size()};
  }

  template <class E>
  auto numpy_texpr_2<E>::fast(long i) const -> decltype(this->arg(
      fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), i))
  {
    return arg(cstride_slice<1>(pythonic::builtins::None, pythonic::builtins::None), i);
  }

  template <class E>
  auto numpy_texpr_2<E>::fast(long i) -> decltype(this->arg(
      fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), i))
  {
    return arg(cstride_slice<1>(pythonic::builtins::None, pythonic::builtins::None), i);
  }

#ifdef USE_XSIMD
  template <class E>
  template <class vectorizer>
  typename numpy_texpr_2<E>::simd_iterator numpy_texpr_2<E>::vbegin(vectorizer) const
  {
    return {*this};
  }

  template <class E>
  template <class vectorizer>
  typename numpy_texpr_2<E>::simd_iterator numpy_texpr_2<E>::vend(vectorizer) const
  {
    return {*this}; // not vectorizable anyway
  }
#endif

  template <class E>
  auto numpy_texpr_2<E>::operator[](long i) const -> decltype(this->fast(i))
  {
    if (i < 0)
      i += size();
    return fast(i);
  }

  template <class E>
  auto numpy_texpr_2<E>::operator[](long i) -> decltype(this->fast(i))
  {
    if (i < 0)
      i += size();
    return fast(i);
  }

  template <class E>
  template <class S>
  auto numpy_texpr_2<E>::operator[](S const &s0) const -> numpy_texpr<decltype(this->arg(
      fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), (s0.step, s0)))>
  {
    return {arg(fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), s0)};
  }

  template <class E>
  template <class S>
  auto numpy_texpr_2<E>::operator[](S const &s0) -> numpy_texpr<decltype(this->arg(
      fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), (s0.step, s0)))>
  {
    return {arg(fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), s0)};
  }

  /* element filtering */
  template <class E>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value == 1 && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_texpr_2<E>, ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::fast(F const &filter) const
  {
    long sz = filter.template shape<0>();
    long *raw = utils::allocate<long>(sz);
    long n = 0;
    for (long i = 0; i < sz; ++i)
      if (filter.fast(i))
        raw[n++] = i;
    // reallocate(raw, n);
    return this->fast(ndarray<long, pshape<long>>(raw, pshape<long>(n), types::ownership::owned));
  }
  template <class E>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value != 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<typename numpy_texpr_2<E>::dtype, pshape<long>>,
                               ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::fast(F const &filter) const
  {
    return numpy::functor::array{}(*this)
        .flat()[ndarray<typename F::dtype, typename F::shape_t>(filter).flat()];
  }

  template <class E>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value == 1 && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_texpr_2<E>, ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::operator[](F const &filter) const
  {
    return fast(filter);
  }

  template <class E>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value != 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<typename numpy_texpr_2<E>::dtype, pshape<long>>,
                               ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::operator[](F const &filter) const
  {
    return fast(filter);
  }

  template <class E>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !std::is_same<bool, typename F::dtype>::value &&
                       !is_pod_array<F>::value,
                   numpy_vexpr<numpy_texpr_2<E>, ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::operator[](F const &filter) const
  {
    static_assert(F::value == 1, "advanced indexing only supporint with 1D index");

    return {*this, filter};
  }

  template <class E>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !std::is_same<bool, typename F::dtype>::value &&
                       !is_pod_array<F>::value,
                   numpy_vexpr<numpy_texpr_2<E>, ndarray<long, pshape<long>>>>
  numpy_texpr_2<E>::fast(F const &filter) const
  {
    static_assert(F::value == 1, "advanced indexing only supported with 1D index");
    return {*this, filter};
  }

  template <class E>
  template <class S0, class... S>
  auto numpy_texpr_2<E>::operator()(S0 const &s0, S const &...s) const -> std::enable_if_t<
      !is_numexpr_arg<S0>::value,
      decltype(this->_reverse_index(std::tuple<S0 const &, S const &...>{s0, s...},
                                    utils::make_reversed_index_sequence<1 + sizeof...(S)>()))>
  {
    return _reverse_index(std::tuple<S0 const &, S const &...>{s0, s...},
                          utils::make_reversed_index_sequence<1 + sizeof...(S)>());
  }
  template <class E>
  template <class S0, class... S>
  auto numpy_texpr_2<E>::operator()(S0 const &s0, S const &...s) const
      -> std::enable_if_t<is_numexpr_arg<S0>::value, decltype(this->copy()(s0, s...))>
  {
    return copy()(s0, s...);
  }

  template <class E>
  numpy_texpr_2<E>::operator bool() const
  {
    return (bool)arg;
  }

  template <class E>
  long numpy_texpr_2<E>::flat_size() const
  {
    return arg.flat_size();
  }

  template <class E>
  intptr_t numpy_texpr_2<E>::id() const
  {
    return arg.id();
  }

  template <class Arg>
  template <class Expr>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator=(Expr const &expr)
  {
    return utils::broadcast_copy < numpy_texpr_2 &, Expr, value, value - utils::dim_of<Expr>::value,
           is_vectorizable && std::is_same<dtype, typename dtype_of<Expr>::type>::value &&
               types::is_vectorizable<Expr>::value > (*this, expr);
  }
  template <class Arg>
  template <class Expr>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator=(numpy_texpr<Expr> const &expr)
  {
    arg = expr.arg;
    return *this;
  }

  template <class Arg>
  template <class Op, class Expr>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::update_(Expr const &expr)
  {
    using BExpr =
        std::conditional_t<std::is_scalar<Expr>::value, broadcast<Expr, dtype>, Expr const &>;
    BExpr bexpr = expr;
    utils::broadcast_update<
        Op, numpy_texpr_2 &, BExpr, value,
        value - (std::is_scalar<Expr>::value + utils::dim_of<Expr>::value),
        is_vectorizable &&
            types::is_vectorizable<std::remove_cv_t<std::remove_reference_t<BExpr>>>::value &&
            std::is_same<dtype, typename dtype_of<std::decay_t<BExpr>>::type>::value>(*this, bexpr);
    return *this;
  }

  template <class Arg>
  template <class Expr>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator+=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::iadd>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator-=(E const &expr)
  {
    return update_<pythonic::operator_::functor::isub>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator*=(E const &expr)
  {
    return update_<pythonic::operator_::functor::imul>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator/=(E const &expr)
  {
    return update_<pythonic::operator_::functor::idiv>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator&=(E const &expr)
  {
    return update_<pythonic::operator_::functor::iand>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator|=(E const &expr)
  {
    return update_<pythonic::operator_::functor::ior>(expr);
  }

  template <class Arg>
  template <class E>
  numpy_texpr_2<Arg> &numpy_texpr_2<Arg>::operator^=(E const &expr)
  {
    return update_<pythonic::operator_::functor::ixor>(expr);
  }

  // only implemented for N = 2

  template <class T, class S0, class S1>
  numpy_texpr<ndarray<T, pshape<S0, S1>>>::numpy_texpr(ndarray<T, pshape<S0, S1>> const &arg)
      : numpy_texpr_2<ndarray<T, pshape<S0, S1>>>{arg}
  {
  }

  template <class T>
  numpy_texpr<ndarray<T, array_tuple<long, 2>>>::numpy_texpr(
      ndarray<T, array_tuple<long, 2>> const &arg)
      : numpy_texpr_2<ndarray<T, array_tuple<long, 2>>>{arg}
  {
  }

  template <class E, class... S>
  numpy_texpr<numpy_gexpr<E, S...>>::numpy_texpr(numpy_gexpr<E, S...> const &arg)
      : numpy_texpr_2<numpy_gexpr<E, S...>>{arg}
  {
  }
} // namespace types
PYTHONIC_NS_END

#endif
