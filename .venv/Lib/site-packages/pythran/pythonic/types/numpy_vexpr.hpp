#ifndef PYTHONIC_TYPES_NUMPY_VEXPR_HPP
#define PYTHONIC_TYPES_NUMPY_VEXPR_HPP

#include "pythonic/utils/allocate.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T, class F>
  template <class E>
  std::enable_if_t<is_iterable<E>::value, numpy_vexpr<T, F> &>
  numpy_vexpr<T, F>::operator=(E const &expr)
  {
    // TODO: avoid the tmp copy when no aliasing
    typename assignable<E>::type tmp{expr};
    for (long i = 0, n = tmp.template shape<0>(); i < n; ++i)
      (*this).fast(i) = tmp.fast(i);
    return *this;
  }
  template <class T, class F>
  template <class E>
  std::enable_if_t<!is_iterable<E>::value, numpy_vexpr<T, F> &>
  numpy_vexpr<T, F>::operator=(E const &expr)
  {
    for (long i = 0, n = shape<0>(); i < n; ++i)
      (*this).fast(i) = expr;
    return *this;
  }
  template <class T, class F>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator=(numpy_vexpr<T, F> const &expr)
  {
    // TODO: avoid the tmp copy when no aliasing
    typename assignable<numpy_vexpr<T, F>>::type tmp{expr};
    for (long i = 0, n = tmp.template shape<0>(); i < n; ++i)
      (*this).fast(i) = tmp.fast(i);
    return *this;
  }

  template <class T, class F>
  typename numpy_vexpr<T, F>::iterator numpy_vexpr<T, F>::begin()
  {
    return {*this, 0};
  }
  template <class T, class F>
  typename numpy_vexpr<T, F>::iterator numpy_vexpr<T, F>::end()
  {
    return {*this, shape<0>()};
  }
  template <class T, class F>
  typename numpy_vexpr<T, F>::const_iterator numpy_vexpr<T, F>::begin() const
  {
    return {*this, 0};
  }
  template <class T, class F>
  typename numpy_vexpr<T, F>::const_iterator numpy_vexpr<T, F>::end() const
  {
    return {*this, shape<0>()};
  }
  template <class T, class F>
  template <class... S>
  auto numpy_vexpr<T, F>::operator()(S const &...slices) const
      -> decltype(ndarray<dtype, array_tuple<long, value>>{*this}(slices...))
  {
    return ndarray<dtype, array_tuple<long, value>>{*this}(slices...);
  }
#ifdef USE_XSIMD
  template <class T, class F>
  template <class vectorizer>
  typename numpy_vexpr<T, F>::simd_iterator numpy_vexpr<T, F>::vbegin(vectorizer) const
  {
    return {*this, 0};
  }

  template <class T, class F>
  template <class vectorizer>
  typename numpy_vexpr<T, F>::simd_iterator numpy_vexpr<T, F>::vend(vectorizer) const
  {
    return {*this, 0};
  }
#endif

  /* element filtering */
  template <class T, class F>
  template <class E> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<E>::value && std::is_same<bool, typename E::dtype>::value &&
                       !is_pod_array<F>::value,
                   numpy_vexpr<numpy_vexpr<T, F>, ndarray<long, pshape<long>>>>
  numpy_vexpr<T, F>::fast(E const &filter) const
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

  template <class T, class F>
  template <class E> // indexing through an array of boolean -- a mask
  std::enable_if_t<!is_slice<E>::value && is_numexpr_arg<E>::value &&
                       std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_vexpr<T, F>, ndarray<long, pshape<long>>>>
  numpy_vexpr<T, F>::operator[](E const &filter) const
  {
    return fast(filter);
  }

  template <class T, class F>
  template <class E> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<E>::value && !is_array_index<E>::value &&
                       !std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_vexpr<T, F>, E>>
  numpy_vexpr<T, F>::operator[](E const &filter) const
  {
    return {*this, filter};
  }

  template <class T, class F>
  template <class E> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<E>::value && !is_array_index<E>::value &&
                       !std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<numpy_vexpr<T, F>, E>>
  numpy_vexpr<T, F>::fast(E const &filter) const
  {
    return (*this)[filter];
  }
  template <class T, class F>
  template <class Op, class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::update_(Expr const &expr)
  {
    using BExpr =
        std::conditional_t<std::is_scalar<Expr>::value, broadcast<Expr, dtype>, Expr const &>;
    BExpr bexpr = expr;
    utils::broadcast_update<
        Op, numpy_vexpr &, BExpr, value,
        value - (std::is_scalar<Expr>::value + utils::dim_of<Expr>::value),
        is_vectorizable &&
            types::is_vectorizable<std::remove_cv_t<std::remove_reference_t<BExpr>>>::value &&
            std::is_same<dtype, typename dtype_of<std::decay_t<BExpr>>::type>::value>(*this, bexpr);
    return *this;
  }
  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator+=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::iadd>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator-=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::isub>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator*=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::imul>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator/=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::idiv>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator&=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::iand>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator|=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::ior>(expr);
  }

  template <class T, class F>
  template <class Expr>
  numpy_vexpr<T, F> &numpy_vexpr<T, F>::operator^=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::ixor>(expr);
  }
} // namespace types
PYTHONIC_NS_END

#endif
