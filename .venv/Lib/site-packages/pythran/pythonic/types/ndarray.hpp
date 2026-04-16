#ifndef PYTHONIC_TYPES_NDARRAY_HPP
#define PYTHONIC_TYPES_NDARRAY_HPP

#include "pythonic/include/types/ndarray.hpp"

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/attr.hpp"
#include "pythonic/types/empty_iterator.hpp"

#include "pythonic/builtins/ValueError.hpp"

#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/broadcast_copy.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/nested_container.hpp"
#include "pythonic/utils/reserve.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include "pythonic/types/list.hpp"
#include "pythonic/types/raw_array.hpp"
#include "pythonic/types/slice.hpp"
#include "pythonic/types/tuple.hpp"

#include "pythonic/numpy/bool_.hpp"
#include "pythonic/numpy/complex128.hpp"
#include "pythonic/numpy/complex64.hpp"
#include "pythonic/numpy/float32.hpp"
#include "pythonic/numpy/float64.hpp"
#include "pythonic/numpy/int16.hpp"
#include "pythonic/numpy/int32.hpp"
#include "pythonic/numpy/int64.hpp"
#include "pythonic/numpy/int8.hpp"
#include "pythonic/numpy/uint16.hpp"
#include "pythonic/numpy/uint32.hpp"
#include "pythonic/numpy/uint64.hpp"
#include "pythonic/numpy/uint8.hpp"

#include "pythonic/types/numpy_expr.hpp"
#include "pythonic/types/numpy_gexpr.hpp"
#include "pythonic/types/numpy_iexpr.hpp"
#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/types/numpy_texpr.hpp"
#include "pythonic/types/numpy_vexpr.hpp"
#include "pythonic/types/vectorizable_type.hpp"
#include "pythonic/utils/array_helper.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/builtins/len.hpp"
#include "pythonic/operator_/iadd.hpp"
#include "pythonic/operator_/iand.hpp"
#include "pythonic/operator_/idiv.hpp"
#include "pythonic/operator_/imul.hpp"
#include "pythonic/operator_/ior.hpp"
#include "pythonic/operator_/isub.hpp"
#include "pythonic/operator_/ixor.hpp"

#include <array>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>

#if !defined(HAVE_SSIZE_T) || !HAVE_SSIZE_T
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif
#endif

PYTHONIC_NS_BEGIN

namespace types
{
  template <class pS, size_t... Is>
  array_tuple<long, std::tuple_size<pS>::value> make_strides(pS const &shape,
                                                             std::index_sequence<Is...>)
  {
    array_tuple<long, std::tuple_size<pS>::value> out;
    out[std::tuple_size<pS>::value - 1] = 1;
    (void)std::initializer_list<long>{
        (out[std::tuple_size<pS>::value - Is - 2] =
             out[std::tuple_size<pS>::value - Is - 1] *
             std::get<std::tuple_size<pS>::value - Is - 1>(shape))...};
    return out;
  }

  template <class pS>
  array_tuple<long, std::tuple_size<pS>::value> make_strides(pS const &shape)
  {
    return make_strides(shape, std::make_index_sequence<std::tuple_size<pS>::value - 1>());
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pS>>::iterator
  type_helper<ndarray<T, pS>>::make_iterator(ndarray<T, pS> &n, long i)
  {
    return {n, i};
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pS>>::const_iterator
  type_helper<ndarray<T, pS>>::make_iterator(ndarray<T, pS> const &n, long i)
  {
    return {n, i};
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, pS>>::initialize_from_iterable(S &shape, T *from, Iter &&iter)
  {
    return type_helper<ndarray<T, pS> const &>::initialize_from_iterable(shape, from,
                                                                         std::forward<Iter>(iter));
  }

  template <class T, class pS>
  numpy_iexpr<ndarray<T, pS>> type_helper<ndarray<T, pS>>::get(ndarray<T, pS> &&self, long i)
  {
    return {std::move(self), i};
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pS> const &>::iterator
  type_helper<ndarray<T, pS> const &>::make_iterator(ndarray<T, pS> &n, long i)
  {
    return {n, i};
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pS> const &>::const_iterator
  type_helper<ndarray<T, pS> const &>::make_iterator(ndarray<T, pS> const &n, long i)
  {
    return {n, i};
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, pS> const &>::initialize_from_iterable(S &shape, T *from, Iter &&iter)
  {
    sutils::assign(std::get<std::tuple_size<S>::value - std::tuple_size<pS>::value>(shape),
                   iter.size());
    for (auto content : iter)
      from = type_helper<ndarray<T, sutils::pop_tail_t<pS>> const &>::initialize_from_iterable(
          shape, from, content);
    return from;
  }

  template <class T, class pS>
  numpy_iexpr<ndarray<T, pS> const &>
  type_helper<ndarray<T, pS> const &>::get(ndarray<T, pS> const &self, long i)
  {
    return numpy_iexpr<ndarray<T, pS> const &>(self, i);
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>>>::iterator
  type_helper<ndarray<T, pshape<pS>>>::make_iterator(ndarray<T, pshape<pS>> &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>>>::const_iterator
  type_helper<ndarray<T, pshape<pS>>>::make_iterator(ndarray<T, pshape<pS>> const &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, pshape<pS>>>::initialize_from_iterable(S &shape, T *from, Iter &&iter)
  {
    sutils::assign(std::get<std::tuple_size<S>::value - 1>(shape), iter.size());
    return std::copy(iter.begin(), iter.end(), from);
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>>>::type
  type_helper<ndarray<T, pshape<pS>>>::get(ndarray<T, pshape<pS>> &&self, long i)
  {
    return self.buffer[i];
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>> const &>::iterator
  type_helper<ndarray<T, pshape<pS>> const &>::make_iterator(ndarray<T, pshape<pS>> &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>> const &>::const_iterator
  make_iterator(ndarray<T, pshape<pS>> const &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, pshape<pS>> const &>::initialize_from_iterable(S &shape, T *from,
                                                                           Iter &&iter)
  {
    sutils::assign(std::get<std::tuple_size<S>::value - 1>(shape), iter.size());
    return std::copy(iter.begin(), iter.end(), from);
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, pshape<pS>> const &>::type &
  type_helper<ndarray<T, pshape<pS>> const &>::get(ndarray<T, pshape<pS>> const &self, long i)
  {
    return self.buffer[i];
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>>>::iterator
  type_helper<ndarray<T, array_tuple<pS, 1>>>::make_iterator(ndarray<T, array_tuple<pS, 1>> &n,
                                                             long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>>>::const_iterator
  type_helper<ndarray<T, array_tuple<pS, 1>>>::make_iterator(
      ndarray<T, array_tuple<pS, 1>> const &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, array_tuple<pS, 1>>>::initialize_from_iterable(S &shape, T *from,
                                                                           Iter &&iter)
  {
    sutils::assign(std::get<std::tuple_size<S>::value - 1>(shape), iter.size());
    return std::copy(iter.begin(), iter.end(), from);
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>>>::type
  type_helper<ndarray<T, array_tuple<pS, 1>>>::get(ndarray<T, array_tuple<pS, 1>> &&self, long i)
  {
    return self.buffer[i];
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>> const &>::iterator
  type_helper<ndarray<T, array_tuple<pS, 1>> const &>::make_iterator(
      ndarray<T, array_tuple<pS, 1>> &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>> const &>::const_iterator
  make_iterator(ndarray<T, array_tuple<pS, 1>> const &n, long i)
  {
    return n.buffer + i;
  }

  template <class T, class pS>
  template <class S, class Iter>
  T *type_helper<ndarray<T, array_tuple<pS, 1>> const &>::initialize_from_iterable(S &shape,
                                                                                   T *from,
                                                                                   Iter &&iter)
  {
    sutils::assign(std::get<std::tuple_size<S>::value - 1>(shape), iter.size());
    return std::copy(iter.begin(), iter.end(), from);
  }

  template <class T, class pS>
  typename type_helper<ndarray<T, array_tuple<pS, 1>> const &>::type &
  type_helper<ndarray<T, array_tuple<pS, 1>> const &>::get(
      ndarray<T, array_tuple<pS, 1>> const &self, long i)
  {
    return self.buffer[i];
  }

  template <class S>
  long patch_index(long index, S const &)
  {
    return index;
  }
  inline long patch_index(long index, std::integral_constant<long, 1> const &)
  {
    return 0;
  }

  template <size_t L>
  template <class S, class Ty, size_t M>
  long noffset<L>::operator()(S const &strides, array_tuple<Ty, M> const &indices) const
  {
    auto index = patch_index(indices[M - L], std::tuple_element_t<M - L, typename S::shape_t>());
    auto offset = noffset<L - 1>{}(strides, indices);
    auto stride = strides.template strides<M - L>();
    return offset + stride * index;
  }

  template <size_t L>
  template <class S, class Ty, size_t M, class pS>
  long noffset<L>::operator()(S const &strides, array_tuple<Ty, M> const &indices,
                              pS const &shape) const
  {
    auto index = patch_index(indices[M - L], std::tuple_element_t<M - L, typename S::shape_t>());
    if (index < 0)
      index += std::get<M - L>(shape);
    assert(0 <= index and index < std::get<M - L>(shape));
    auto offset = noffset<L - 1>{}(strides, indices, shape);
    auto stride = strides.template strides<M - L>();
    return offset + stride * index;
  }

  template <>
  template <class S, class Ty, size_t M>
  long noffset<1>::operator()(S const &strides, array_tuple<Ty, M> const &indices) const
  {
    auto index = patch_index(indices[M - 1], std::tuple_element_t<M - 1, typename S::shape_t>());
    return strides.template strides<M - 1>() * index;
  }

  template <>
  template <class S, class Ty, size_t M, class pS>
  long noffset<1>::operator()(S const &strides, array_tuple<Ty, M> const &indices,
                              pS const &shape) const
  {
    auto index = patch_index(indices[M - 1], std::tuple_element_t<M - 1, typename S::shape_t>());
    if (index < 0)
      index += std::get<M - 1>(shape);
    assert(0 <= index && index < std::get<M - 1>(shape));
    auto stride = strides.template strides<M - 1>();
    return stride * index;
  }

  /* constructors */
  template <class T, class pS>
  ndarray<T, pS>::ndarray() : mem(utils::no_memory()), buffer(nullptr), _shape(), _strides()
  {
  }

  /* from other memory */
  template <class T, class pS>
  ndarray<T, pS>::ndarray(utils::shared_ref<raw_array<T>> const &mem, pS const &shape)
      : mem(mem), buffer(mem->data), _shape(shape), _strides(make_strides(shape))
  {
  }
  template <class T, class pS>
  ndarray<T, pS>::ndarray(utils::shared_ref<raw_array<T>> &&mem, pS const &shape)
      : mem(std::move(mem)), buffer(this->mem->data), _shape(shape), _strides(make_strides(shape))
  {
  }

  /* from other array */
  template <class T, class pS>
  template <class Tp, class pSp>
  ndarray<T, pS>::ndarray(ndarray<Tp, pSp> const &other)
      : mem(other.flat_size()), buffer(mem->data), _shape(other._shape), _strides(other._strides)
  {
    static_assert(std::tuple_size<pS>::value == std::tuple_size<pSp>::value, "compatible shapes");
    std::copy(other.fbegin(), other.fend(), fbegin());
  }

  template <class T, class pS>
  template <class pSp>
  ndarray<T, pS>::ndarray(ndarray<T, pSp> const &other)
      : mem(other.mem), buffer(mem->data), _shape(other._shape), _strides(other._strides)
  {
    static_assert(std::tuple_size<pS>::value == std::tuple_size<pSp>::value, "compatible shapes");
  }

  /* from a seed */
  template <class T, class pS>
  ndarray<T, pS>::ndarray(pS const &shape, none_type init)
      : mem(sutils::sprod(shape)), buffer(mem->data), _shape(shape), _strides(make_strides(shape))
  {
  }

  template <class T, class pS>
  ndarray<T, pS>::ndarray(pS const &shape, T init) : ndarray(shape, none_type())
  {
    std::fill(fbegin(), fend(), init);
  }

  /* from a foreign pointer */
  template <class T, class pS>
  template <class S>
  ndarray<T, pS>::ndarray(T *data, S const *pshape, ownership o)
      : mem(data, o), buffer(mem->data), _shape(pshape)
  {
    _strides = make_strides(_shape);
  }
  template <class T, class pS>
  ndarray<T, pS>::ndarray(T *data, pS const &pshape, ownership o)
      : mem(data, o), buffer(mem->data), _shape(pshape)
  {
    _strides = make_strides(_shape);
  }

#ifdef ENABLE_PYTHON_MODULE
  template <class T, class pS>
  template <class S>
  ndarray<T, pS>::ndarray(T *data, S const *pshape, PyObject *obj_ptr)
      : ndarray(data, pshape, ownership::external)
  {
    mem.external(obj_ptr); // mark memory as external to decref at the end of
                           // its lifetime
  }
  template <class T, class pS>
  ndarray<T, pS>::ndarray(T *data, pS const &pshape, PyObject *obj_ptr)
      : ndarray(data, pshape, ownership::external)
  {
    mem.external(obj_ptr); // mark memory as external to decref at the end of
                           // its lifetime
  }

#endif

  template <class T, class pS>
  template <class Iterable, class>
  ndarray<T, pS>::ndarray(Iterable &&iterable)
      : mem(utils::nested_container_size<Iterable>::flat_size(std::forward<Iterable>(iterable))),
        buffer(mem->data), _shape()
  {
    type_helper<ndarray>::initialize_from_iterable(_shape, mem->data,
                                                   std::forward<Iterable>(iterable));
    _strides = make_strides(_shape);
  }

  /* from a  numpy expression */
  template <class T, class pS>
  template <class E>
  void ndarray<T, pS>::initialize_from_expr(E const &expr)
  {
    assert(buffer);
    utils::broadcast_copy<ndarray &, E, value, 0,
                          is_vectorizable && E::is_vectorizable &&
                              std::is_same<dtype, typename E::dtype>::value>(*this, expr);
  }

  template <class T, class pS>
  template <class Op, class... Args>
  ndarray<T, pS>::ndarray(numpy_expr<Op, Args...> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  template <class T, class pS>
  template <class Arg>
  ndarray<T, pS>::ndarray(numpy_texpr<Arg> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  template <class T, class pS>
  template <class Arg>
  ndarray<T, pS>::ndarray(numpy_texpr_2<Arg> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  template <class T, class pS>
  template <class Arg, class... S>
  ndarray<T, pS>::ndarray(numpy_gexpr<Arg, S...> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  template <class T, class pS>
  template <class Arg>
  ndarray<T, pS>::ndarray(numpy_iexpr<Arg> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  template <class T, class pS>
  template <class Arg, class F>
  ndarray<T, pS>::ndarray(numpy_vexpr<Arg, F> const &expr)
      : mem(expr.flat_size()), buffer(mem->data), _shape(sutils::getshape(expr)),
        _strides(make_strides(_shape))
  {
    initialize_from_expr(expr);
  }

  /* update operators */

  template <class T, class pS>
  template <class Op, class Expr>
  ndarray<T, pS> &ndarray<T, pS>::update_(Expr const &expr)
  {
    using BExpr = std::conditional_t<std::is_scalar<Expr>::value, broadcast<Expr, T>, Expr const &>;
    BExpr bexpr = expr;
    utils::broadcast_update<
        Op, ndarray &, BExpr, value,
        value - (std::is_scalar<Expr>::value + utils::dim_of<Expr>::value),
        is_vectorizable &&
            types::is_vectorizable<std::remove_cv_t<std::remove_reference_t<BExpr>>>::value &&
            std::is_same<dtype, typename dtype_of<std::decay_t<BExpr>>::type>::value>(*this, bexpr);
    return *this;
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator+=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::iadd>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator-=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::isub>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator*=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::imul>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator/=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::idiv>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator&=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::iand>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator|=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::ior>(expr);
  }

  template <class T, class pS>
  template <class Expr>
  ndarray<T, pS> &ndarray<T, pS>::operator^=(Expr const &expr)
  {
    return update_<pythonic::operator_::functor::ixor>(expr);
  }

  /* element indexing
   * differentiate const from non const, && r-value from l-value
   * */

  template <class T, class pS>
  template <class Ty>
  std::enable_if_t<std::is_integral<Ty>::value, T &>
  ndarray<T, pS>::fast(array_tuple<Ty, value> const &indices)
  {
    assert(inbound_indices(indices));
    return *(buffer + noffset<std::tuple_size<pS>::value>{}(*this, indices));
  }

  template <class T, class pS>
  template <class Ty>
  std::enable_if_t<std::is_integral<Ty>::value, T>
  ndarray<T, pS>::fast(array_tuple<Ty, value> const &indices) const
  {
    assert(inbound_indices(indices));
    return *(buffer + noffset<std::tuple_size<pS>::value>{}(*this, indices));
  }

  template <class T, class pS>
  template <class Ty, size_t M>
  auto ndarray<T, pS>::fast(array_tuple<Ty, M> const &indices) const & -> std::enable_if_t<
      std::is_integral<Ty>::value, decltype(nget<M - 1>().fast(*this, indices))>
  {
    return nget<M - 1>().fast(*this, indices);
  }

  template <class T, class pS>
  template <class Ty, size_t M>
  auto ndarray<T, pS>::fast(array_tuple<Ty, M> const &indices) && -> std::enable_if_t<
      std::is_integral<Ty>::value, decltype(nget<M - 1>().fast(std::move(*this), indices))>
  {
    return nget<M - 1>().fast(std::move(*this), indices);
  }

  template <class T, class pS>
  template <class Ty>
  std::enable_if_t<std::is_integral<Ty>::value, T const &>
  ndarray<T, pS>::operator[](array_tuple<Ty, value> const &indices) const
  {
    return *(buffer + noffset<std::tuple_size<pS>::value>{}(*this, indices, _shape));
  }

  template <class T, class pS>
  template <class Ty>
  std::enable_if_t<std::is_integral<Ty>::value, T &>
  ndarray<T, pS>::operator[](array_tuple<Ty, value> const &indices)
  {
    return *(buffer + noffset<std::tuple_size<pS>::value>{}(*this, indices, _shape));
  }

  template <class T, class pS>
  template <class Ty, size_t M>
  auto ndarray<T, pS>::operator[](array_tuple<Ty, M> const &indices) const
      & -> std::enable_if_t<std::is_integral<Ty>::value, decltype(nget<M - 1>()(*this, indices))>
  {
    return nget<M - 1>()(*this, indices);
  }

  template <class T, class pS>
  template <class Ty, size_t M>
  auto ndarray<T, pS>::operator[](array_tuple<Ty, M> const &indices) && -> std::enable_if_t<
      std::is_integral<Ty>::value, decltype(nget<M - 1>()(std::move(*this), indices))>
  {
    return nget<M - 1>()(std::move(*this), indices);
  }

#ifdef USE_XSIMD
  template <class T, class pS>
  template <class vectorizer>
  typename ndarray<T, pS>::simd_iterator ndarray<T, pS>::vbegin(vectorizer) const
  {
    return {buffer};
  }

  template <class T, class pS>
  template <class vectorizer>
  typename ndarray<T, pS>::simd_iterator ndarray<T, pS>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {buffer + long(std::get<0>(_shape) / vector_size * vector_size)};
  }

#endif

  /* slice indexing */
  template <class T, class pS>
  ndarray<T, sutils::push_front_t<pS, std::integral_constant<long, 1>>>
  ndarray<T, pS>::operator[](none_type) const
  {
    sutils::push_front_t<pS, std::integral_constant<long, 1>> new_shape;
    sutils::copy_shape<1, -1>(new_shape, *this,
                              std::make_index_sequence<std::tuple_size<pS>::value>());
    return reshape(new_shape);
  }

  template <class T, class pS>
  template <class S>
  std::enable_if_t<is_slice<S>::value, numpy_gexpr<ndarray<T, pS> const &, normalize_t<S>>>
  ndarray<T, pS>::operator[](S const &s) const &
  {
    return make_gexpr(*this, s);
  }

  template <class T, class pS>
  template <class S>
  std::enable_if_t<is_slice<S>::value, numpy_gexpr<ndarray<T, pS>, normalize_t<S>>>
  ndarray<T, pS>::operator[](S const &s) &&
  {
    return make_gexpr(std::move(*this), s);
  }

  template <class T, class pS>
  long ndarray<T, pS>::size() const
  {
    return std::get<0>(_shape);
  }

  /* extended slice indexing */
  template <class T, class pS>
  template <class S0, class... S>
  auto ndarray<T, pS>::operator()(S0 const &s0, S const &...s)
      const & -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...))
  {
    return extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...);
  }

  template <class T, class pS>
  template <class S0, class... S>
  auto ndarray<T, pS>::operator()(S0 const &s0, S const &...s)
      & -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...))
  {
    return extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...);
  }

  template <class T, class pS>
  template <class S0, class... S>
  auto ndarray<T, pS>::operator()(S0 const &s0, S const &...s)
      && -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}(std::move(*this), s0, s...))
  {
    return extended_slice<count_new_axis<S0, S...>::value>{}(std::move(*this), s0, s...);
  }

  /* element filtering */
  template <class T, class pS>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value == 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pS>, ndarray<long, pshape<long>>>>
  ndarray<T, pS>::fast(F const &filter) const
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

  template <class T, class pS>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value == 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pS>, ndarray<long, pshape<long>>>>
  ndarray<T, pS>::operator[](F const &filter) const
  {
    return fast(filter);
  }
  template <class T, class pS>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value != 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pshape<long>>, ndarray<long, pshape<long>>>>
  ndarray<T, pS>::fast(F const &filter) const
  {
    return flat()[ndarray<typename F::dtype, typename F::shape_t>(filter).flat()];
  }

  template <class T, class pS>
  template <class F> // indexing through an array of boolean -- a mask
  std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                       F::value != 1 && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pshape<long>>, ndarray<long, pshape<long>>>>
  ndarray<T, pS>::operator[](F const &filter) const
  {
    return fast(filter);
  }

  template <class T, class pS>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                       !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pS>, F>>
  ndarray<T, pS>::operator[](F const &filter) const
  {
    return {*this, filter};
  }

  template <class T, class pS>
  template <class F> // indexing through an array of indices -- a view
  std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                       !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                   numpy_vexpr<ndarray<T, pS>, F>>
  ndarray<T, pS>::fast(F const &filter) const
  {
    return {*this, filter};
  }

  template <class T, class pS>
  template <class Ty0, class Ty1, class... Tys, class _>
  auto ndarray<T, pS>::operator[](std::tuple<Ty0, Ty1, Tys...> const &indices) const
      -> std::enable_if_t<is_numexpr_arg<Ty0>::value,
                          decltype(this->_fwdindex(indices,
                                                   std::make_index_sequence<2 + sizeof...(Tys)>()))>
  {
    return _fwdindex(indices, std::make_index_sequence<2 + sizeof...(Tys)>());
  }

  /* through iterators */
  template <class T, class pS>
  typename ndarray<T, pS>::iterator ndarray<T, pS>::begin()
  {
    return type_helper<ndarray>::make_iterator(*this, 0);
  }

  template <class T, class pS>
  typename ndarray<T, pS>::const_iterator ndarray<T, pS>::begin() const
  {
    return type_helper<ndarray>::make_iterator(*this, 0);
  }

  template <class T, class pS>
  typename ndarray<T, pS>::iterator ndarray<T, pS>::end()
  {
    return type_helper<ndarray>::make_iterator(*this, std::get<0>(_shape));
  }

  template <class T, class pS>
  typename ndarray<T, pS>::const_iterator ndarray<T, pS>::end() const
  {
    return type_helper<ndarray>::make_iterator(*this, std::get<0>(_shape));
  }

  template <class T, class pS>
  typename ndarray<T, pS>::const_flat_iterator ndarray<T, pS>::fbegin() const
  {
    return buffer;
  }

  template <class T, class pS>
  typename ndarray<T, pS>::const_flat_iterator ndarray<T, pS>::fend() const
  {
    return buffer + flat_size();
  }

  template <class T, class pS>
  typename ndarray<T, pS>::flat_iterator ndarray<T, pS>::fbegin()
  {
    return buffer;
  }

  template <class T, class pS>
  typename ndarray<T, pS>::flat_iterator ndarray<T, pS>::fend()
  {
    return buffer + flat_size();
  }

  /* member functions */
  template <class T, class pS>
  long ndarray<T, pS>::flat_size() const
  {
    return sutils::prod(*this);
  }
  template <class T, class pS>
  bool ndarray<T, pS>::may_overlap(ndarray const &expr) const
  {
    return id() == expr.id();
  }

  template <class T, class pS>
  template <class qS>
  ndarray<T, qS> ndarray<T, pS>::reshape(qS const &shape) const &
  {
    return {mem, shape};
  }

  template <class T, class pS>
  template <class qS>
  ndarray<T, qS> ndarray<T, pS>::reshape(qS const &shape) &&
  {
    return {std::move(mem), shape};
  }

  template <class T, class pS>
  ndarray<T, pS>::operator bool() const
  {
    if (sutils::any_of(*this, [](long n) { return n != 1; }))
      throw ValueError("The truth value of an array with more than one element "
                       "is ambiguous. Use a.any() or a.all()");
    return *buffer;
  }

  template <class T, class pS>
  ndarray<T, pshape<long>> ndarray<T, pS>::flat() const
  {
    return {mem, pshape<long>{{flat_size()}}};
  }

  template <class T, class pS>
  ndarray<T, pS> ndarray<T, pS>::copy() const
  {
    ndarray<T, pS> res(_shape, builtins::None);
    std::copy(fbegin(), fend(), res.fbegin());
    return res;
  }

  template <class T, class pS>
  intptr_t ndarray<T, pS>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*mem));
  }

  /* pretty printing { */
  namespace impl
  {
    template <class T, class pS>
    size_t get_spacing(ndarray<T, pS> const &e)
    {
      std::ostringstream oss;
      if (e.flat_size()) {
        oss << *std::max_element(e.fbegin(), e.fend());
        size_t s = oss.str().length();
        for (auto iter = e.fbegin(), end = e.fend(); iter != end; ++iter) {
          oss.str("");
          oss.width(s);
          oss << *iter;
          size_t ts = oss.str().length();
          if (ts > s)
            s = ts;
        }
        return s;
      }
      return 0;
    }
    template <class T, class pS>
    size_t get_spacing(ndarray<std::complex<T>, pS> const &e)
    {
      std::ostringstream oss;
      if (e.flat_size())
        oss << *e.fbegin();
      return oss.str().length() + 2;
    }
  } // namespace impl

  template <class T, class pS>
  std::ostream &operator<<(std::ostream &os, ndarray<T, pS> const &e)
  {
    std::array<long, std::tuple_size<pS>::value> strides;
    auto shape = sutils::getshape(e);
    strides[std::tuple_size<pS>::value - 1] = std::get<std::tuple_size<pS>::value - 1>(shape);
    if (strides[std::tuple_size<pS>::value - 1] == 0)
      return os << "[]";
    std::transform(strides.rbegin(), strides.rend() - 1, shape.rbegin() + 1, strides.rbegin() + 1,
                   std::multiplies<long>());
    size_t depth = std::tuple_size<pS>::value;
    int step = -1;
    size_t size = impl::get_spacing(e);
    auto iter = e.fbegin();
    int max_modulo = 1000;

    os << "[";
    if (std::get<0>(shape) != 0)
      do {
        if (depth == 1) {
          os.width(size);
          os << *iter++;
          for (int i = 1; i < std::get<std::tuple_size<pS>::value - 1>(shape); i++) {
            os.width(size + 1);
            os << *iter++;
          }
          step = 1;
          depth++;
          max_modulo = std::lower_bound(strides.begin(), strides.end(), iter - e.buffer,
                                        [](int comp, int val) { return val % comp != 0; }) -
                       strides.begin();
        } else if (max_modulo + depth == std::tuple_size<pS>::value + 1) {
          depth--;
          step = -1;
          os << "]";
          for (size_t i = 0; i < depth; i++)
            os << std::endl;
          for (size_t i = 0; i < std::tuple_size<pS>::value - depth; i++)
            os << " ";
          os << "[";
        } else {
          depth += step;
          if (step == 1)
            os << "]";
          else
            os << "[";
        }
      } while (depth != std::tuple_size<pS>::value + 1);

    return os << "]";
  }

  template <class E>
  std::enable_if_t<is_array<E>::value, std::ostream &> operator<<(std::ostream &os, E const &e)
  {
    return os << ndarray<typename E::dtype, typename E::shape_t>{e};
  }

  /* } */
  template <class T>
  template <class pS>
  list<T> &list<T>::operator=(ndarray<T, pshape<pS>> const &other)
  {
    _data = utils::shared_ref<T>(other.begin(), other.end());
    return *this;
  }
} // namespace types
PYTHONIC_NS_END

/* std::get overloads */
namespace std
{

  template <size_t I, class E>
  auto get(E &&a) -> std::enable_if_t<
      pythonic::types::is_array<std::remove_cv_t<std::remove_reference_t<E>>>::value,
      decltype(std::forward<E>(a)[I])>
  {
    return std::forward<E>(a)[I];
  }
} // namespace std

/* pythran attribute system { */
#include "pythonic/numpy/transpose.hpp"
PYTHONIC_NS_BEGIN
namespace builtins
{
  namespace details
  {
    template <size_t N>
    template <class E, class... S>
    auto _build_gexpr<N>::operator()(E const &a, S const &...slices)
        -> decltype(_build_gexpr<N - 1>{}(a, types::cstride_slice<1>(), slices...))
    {
      return _build_gexpr<N - 1>{}(a, types::cstride_slice<1>(0, a.size()), slices...);
    }

    template <class E, class... S>
    auto _build_gexpr<1>::operator()(E const &a, S const &...slices) -> decltype(E(a)(slices...))
    {
      return E(a)(slices...);
    }

    template <class E>
    E _make_real(E const &a, utils::int_<0>)
    {
      return a;
    }

    template <class E>
    auto _make_real(E const &a, utils::int_<1>) -> decltype(_build_gexpr<E::value>{}(
        types::ndarray<typename types::is_complex<typename E::dtype>::type,
                       types::array_tuple<long, E::value>>{},
        types::slice()))
    {
      using stype = typename types::is_complex<typename E::dtype>::type;
      auto new_shape = sutils::getshape(a);
      std::get<E::value - 1>(new_shape) *= 2;
      // this is tricky and dangerous!
      auto translated_mem =
          reinterpret_cast<utils::shared_ref<types::raw_array<stype>> const &>(a.mem);
      types::ndarray<stype, types::array_tuple<long, E::value>> translated{translated_mem,
                                                                           new_shape};
      return _build_gexpr<E::value>{}(translated,
                                      types::slice{0, std::get<E::value - 1>(new_shape), 2});
    }

    template <class Op, class... Args>
    auto _make_real(types::numpy_expr<Op, Args...> const &a, utils::int_<1>)
        -> decltype(_make_real(types::ndarray<typename types::numpy_expr<Op, Args...>::dtype,
                                              typename types::numpy_expr<Op, Args...>::shape_t>(a),
                               utils::int_<1>{}))
    {
      return _make_real(types::ndarray<typename types::numpy_expr<Op, Args...>::dtype,
                                       typename types::numpy_expr<Op, Args...>::shape_t>(a),
                        utils::int_<1>{});
    }

    template <class E>
    auto _make_real(types::numpy_iexpr<E> const &a, utils::int_<1>)
        -> decltype(_build_gexpr<types::numpy_iexpr<E>::value>{}(
            std::declval<types::ndarray<
                typename types::is_complex<typename types::numpy_iexpr<E>::dtype>::type,
                types::array_tuple<long, types::numpy_iexpr<E>::value + 1>>>(),
            long(), types::slice()))
    {
      constexpr size_t value = types::numpy_iexpr<E>::value;
      using stype = typename types::is_complex<typename types::numpy_iexpr<E>::dtype>::type;
      auto new_shape = sutils::getshape(a.arg);
      std::get<value>(new_shape) *= 2;
      // this is tricky and dangerous!
      auto translated_mem =
          reinterpret_cast<utils::shared_ref<types::raw_array<stype>> const &>(a.arg.mem);
      types::ndarray<stype, types::array_tuple<long, value + 1>> translated{translated_mem,
                                                                            new_shape};
      long offset = (a.buffer - a.arg.buffer) / a.arg.template strides<0>();
      return _build_gexpr<value>{}(translated, offset,
                                   types::slice{0, std::get<value>(new_shape), 2});
    }

    template <class E>
    types::ndarray<typename E::dtype, typename E::shape_t> _make_imag(E const &a, utils::int_<0>)
    {
      // cannot use numpy.zero: forward declaration issue
      return {utils::callocate<typename E::dtype>(a.flat_size()), sutils::getshape(a),
              types::ownership::owned};
    }

    template <class Op, class... Args>
    auto _make_imag(types::numpy_expr<Op, Args...> const &a, utils::int_<1>)
        -> decltype(_make_imag(types::ndarray<typename types::numpy_expr<Op, Args...>::dtype,
                                              typename types::numpy_expr<Op, Args...>::shape_t>(a),
                               utils::int_<1>{}))
    {
      return _make_imag(types::ndarray<typename types::numpy_expr<Op, Args...>::dtype,
                                       typename types::numpy_expr<Op, Args...>::shape_t>(a),
                        utils::int_<1>{});
    }

    template <class E>
    auto _make_imag(types::numpy_iexpr<E> const &a, utils::int_<1>)
        -> decltype(_build_gexpr<types::numpy_iexpr<E>::value>{}(
            std::declval<types::ndarray<
                typename types::is_complex<typename types::numpy_iexpr<E>::dtype>::type,
                types::array_tuple<long, types::numpy_iexpr<E>::value + 1>>>(),
            long(), types::slice()))
    {
      constexpr size_t value = types::numpy_iexpr<E>::value;
      using stype = typename types::is_complex<typename types::numpy_iexpr<E>::dtype>::type;
      auto new_shape = sutils::getshape(a.arg);
      std::get<types::numpy_iexpr<E>::value>(new_shape) *= 2;
      // this is tricky and dangerous!
      auto translated_mem =
          reinterpret_cast<utils::shared_ref<types::raw_array<stype>> const &>(a.arg.mem);
      types::ndarray<stype, types::array_tuple<long, value + 1>> translated{translated_mem,
                                                                            new_shape};
      long offset = (a.buffer - a.arg.buffer) / a.arg.template strides<0>();
      return _build_gexpr<value>{}(translated, offset,
                                   types::slice{1, std::get<value>(new_shape), 2});
    }

    template <class E>
    auto _make_imag(E const &a, utils::int_<1>) -> decltype(_build_gexpr<E::value>{}(
        types::ndarray<typename types::is_complex<typename E::dtype>::type,
                       types::array_tuple<long, E::value>>{},
        types::slice()))
    {
      using stype = typename types::is_complex<typename E::dtype>::type;
      auto new_shape = sutils::getshape(a);
      std::get<E::value - 1>(new_shape) *= 2;
      // this is tricky and dangerous!
      auto translated_mem =
          reinterpret_cast<utils::shared_ref<types::raw_array<stype>> const &>(a.mem);
      types::ndarray<stype, types::array_tuple<long, E::value>> translated{translated_mem,
                                                                           new_shape};
      return _build_gexpr<E::value>{}(translated,
                                      types::slice{1, std::get<E::value - 1>(new_shape), 2});
    }
  } // namespace details

  template <class E>
  types::array_tuple<long, E::value> getattr(types::attr::SHAPE, E const &a)
  {
    return sutils::getshape(a);
  }

  template <class E>
  constexpr decltype(long(E::value)) getattr(types::attr::NDIM, E const &a)
  {
    return E::value;
  }

  template <class E>
  std::integral_constant<long, E::value> getattr(types::attr::NDIM, E *const &a)
  {
    return {};
  }

  template <class E>
  types::array_tuple<long, E::value> getattr(types::attr::STRIDES, E const &a)
  {
    types::array_tuple<long, E::value> strides;
    strides[E::value - 1] = sizeof(typename E::dtype);
    auto shape = sutils::getshape(a);
    std::transform(strides.rbegin(), strides.rend() - 1, shape.rbegin(), strides.rbegin() + 1,
                   std::multiplies<long>());
    return strides;
  }

  template <class E>
  long getattr(types::attr::SIZE, E const &a)
  {
    return a.flat_size();
  }

  template <class E>
  constexpr long getattr(types::attr::ITEMSIZE, E const &a)
  {
    return sizeof(typename E::dtype);
  }

  template <class E>
  std::integral_constant<long, sizeof(typename E::dtype)> getattr(types::attr::ITEMSIZE,
                                                                  E *const &a)
  {
    return {};
  }

  template <class E>
  long getattr(types::attr::NBYTES, E const &a)
  {
    return a.flat_size() * sizeof(typename E::dtype);
  }

  template <class E>
  auto getattr(types::attr::FLAT, E const &a) -> decltype(a.flat())
  {
    return a.flat();
  }

  template <class T, class pS>
  auto getattr(types::attr::REAL, types::ndarray<T, pS> const &a)
      -> decltype(details::_make_real(a, utils::int_<types::is_complex<T>::value>{}))
  {
    return details::_make_real(a, utils::int_<types::is_complex<T>::value>{});
  }

  template <class E>
  auto getattr(types::attr::REAL, types::numpy_iexpr<E> const &e) -> decltype(details::_make_real(
      e, utils::int_<types::is_complex<typename types::numpy_iexpr<E>::dtype>::value>{}))
  {
    return details::_make_real(
        e, utils::int_<types::is_complex<typename types::numpy_iexpr<E>::dtype>::value>{});
  }

  template <class Op, class... Args>
  auto getattr(types::attr::REAL, types::numpy_expr<Op, Args...> const &a)
      -> decltype(details::_make_real(
          a,
          utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{}))
  {
    return details::_make_real(
        a, utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{});
  }

  template <class E>
  auto getattr(types::attr::REAL, types::numpy_texpr<E> const &a)
      -> decltype(types::numpy_texpr<decltype(getattr(types::attr::REAL{}, a.arg))>{
          getattr(types::attr::REAL{}, a.arg)})
  {
    auto ta = getattr(types::attr::REAL{}, a.arg);
    return types::numpy_texpr<decltype(ta)>{ta};
  }

  template <class T, class pS>
  auto getattr(types::attr::IMAG, types::ndarray<T, pS> const &a)
      -> decltype(details::_make_imag(a, utils::int_<types::is_complex<T>::value>{}))
  {
    return details::_make_imag(a, utils::int_<types::is_complex<T>::value>{});
  }

  template <class Op, class... Args>
  auto getattr(types::attr::IMAG, types::numpy_expr<Op, Args...> const &a)
      -> decltype(details::_make_imag(
          a,
          utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{}))
  {
    return details::_make_imag(
        a, utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{});
  }

  template <class E>
  auto getattr(types::attr::IMAG, types::numpy_texpr<E> const &a)
      -> decltype(types::numpy_texpr<decltype(getattr(types::attr::IMAG{}, a.arg))>{
          getattr(types::attr::IMAG{}, a.arg)})
  {
    auto ta = getattr(types::attr::IMAG{}, a.arg);
    return types::numpy_texpr<decltype(ta)>{ta};
  }

  template <class E>
  types::dtype_t<typename types::dtype_of<E>::type> getattr(types::attr::DTYPE, E const &a)
  {
    return {};
  }
} // namespace builtins
PYTHONIC_NS_END

/* } */

#include "pythonic/types/numpy_operators.hpp"

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/types/int.hpp"

PYTHONIC_NS_BEGIN

/* wrapper around Python array creation
 * its purpose is to hide the difference between the shape stored in pythran
 * (aka long) && the shape stored in numpy (aka npy_intp)
 * it should work (with an extra copy) on 32 bit architecture && without copy
 * on 64 bits architecture
 */
template <class T, size_t N>
struct pyarray_new {

  static_assert(!std::is_same<T, npy_intp>::value, "correctly specialized");

  PyObject *from_descr(PyTypeObject *subtype, PyArray_Descr *descr, T *dims, void *data, int flags,
                       PyObject *obj)
  {
    npy_intp shape[N];
    std::copy(dims, dims + N, shape);
    return pyarray_new<npy_intp, N>{}.from_descr(subtype, descr, shape, data, flags, obj);
  }
  PyObject *from_data(T *dims, int typenum, void *data)
  {
    npy_intp shape[N];
    std::copy(dims, dims + N, shape);
    return pyarray_new<npy_intp, N>{}.from_data(shape, typenum, data);
  }
};

template <size_t N>
struct pyarray_new<npy_intp, N> {

  PyObject *from_descr(PyTypeObject *subtype, PyArray_Descr *descr, npy_intp *dims, void *data,
                       int flags, PyObject *obj)
  {
    return PyArray_NewFromDescr(subtype, descr, N, dims, nullptr, data, flags, obj);
  }

  PyObject *from_data(npy_intp *dims, int typenum, void *data)
  {

    return PyArray_SimpleNewFromData(N, dims, typenum, data);
  }
};

inline void wrapfree(PyObject *capsule)
{
  void *obj = PyCapsule_GetPointer(capsule, PyCapsule_GetName(capsule));
  free(obj);
};

template <class T, class pS>
PyObject *to_python<types::ndarray<T, pS>>::convert(types::ndarray<T, pS> const &cn, bool transpose)
{
  types::ndarray<T, pS> &n = const_cast<types::ndarray<T, pS> &>(cn);
  if (PyObject *p = n.mem.get_foreign()) {
    PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(p);
    auto const *pshape = PyArray_DIMS(arr);
    Py_INCREF(p);

    // handle complex trick :-/
    if ((long)sizeof(T) != PyArray_ITEMSIZE((PyArrayObject *)(arr))) {
      arr = (PyArrayObject *)PyArray_View(
          (PyArrayObject *)(arr), PyArray_DescrFromType(c_type_to_numpy_type<T>::value), nullptr);
    }

    if (sutils::equals(n, pshape)) {
      if (transpose && !(PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS)) {
        PyObject *Transposed = PyArray_Transpose(arr, nullptr);
        Py_DECREF(arr);
        return Transposed;
      } else
        return p;
    } else if (sutils::requals(n, pshape)) {
      if (transpose)
        return p;
      else {
        PyObject *Transposed = PyArray_Transpose(arr, nullptr);
        Py_DECREF(arr);
        return Transposed;
      }
    } else {
      Py_INCREF((PyObject *)PyArray_DESCR(arr));
      auto array = sutils::array(n._shape);
      auto *res = pyarray_new<long, std::tuple_size<pS>::value>{}.from_descr(
          Py_TYPE((PyObject *)arr), PyArray_DESCR(arr), array.data(), PyArray_DATA(arr),
          PyArray_FLAGS(arr) & ~NPY_ARRAY_OWNDATA, p);
      if (transpose && (PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS)) {
        PyObject *Transposed = PyArray_Transpose(reinterpret_cast<PyArrayObject *>(arr), nullptr);
        Py_DECREF(arr);
        return Transposed;
      } else
        return res;
    }
  } else {
    auto array = sutils::array(n._shape);
    PyObject *result = pyarray_new<long, std::tuple_size<pS>::value>{}.from_data(
        array.data(), c_type_to_numpy_type<T>::value, n.buffer);
    if (!result)
      return nullptr;
    // Take responsibility for n.buffer by wrapping it in a capsule and
    // setting result.base to the capsule
    PyObject *capsule = PyCapsule_New(n.buffer, "wrapped_data", (PyCapsule_Destructor)&wrapfree);
    if (!capsule) {
      Py_DECREF(result);
      return nullptr;
    }
    n.mark_memory_external(result);
    Py_INCREF(result); // because it's going to be decrefed when n is destroyed
    if (PyArray_SetBaseObject(reinterpret_cast<PyArrayObject *>(result), capsule) == -1) {
      Py_DECREF(result);
      Py_DECREF(capsule); // will free n.buffer
      return nullptr;
    }
    if (transpose) {
      PyObject *Transposed = PyArray_Transpose(reinterpret_cast<PyArrayObject *>(result), nullptr);
      Py_DECREF(result);
      return Transposed;
    } else
      return result;
  }
}

template <class Arg>
PyObject *to_python<types::numpy_iexpr<Arg>>::convert(types::numpy_iexpr<Arg> const &v,
                                                      bool transpose)
{
  PyObject *res = ::to_python(types::ndarray<typename types::numpy_iexpr<Arg>::dtype,
                                             typename types::numpy_iexpr<Arg>::shape_t>(v));
  if (transpose) {
    PyObject *Transposed = PyArray_Transpose(reinterpret_cast<PyArrayObject *>(res), nullptr);
    Py_DECREF(res);
    return Transposed;
  } else
    return res;
}

template <class Arg, class... S>
PyObject *to_python<types::numpy_gexpr<Arg, S...>>::convert(types::numpy_gexpr<Arg, S...> const &v,
                                                            bool transpose)
{
  PyObject *slices =
      (sizeof...(S) == 1) ? ::to_python(std::get<0>(v.slices)) : ::to_python(v.slices);
  PyObject *base = ::to_python(v.arg);
  PyObject *res = PyObject_GetItem(base, slices);
  Py_DECREF(slices);
  Py_DECREF(base);
  if (transpose) {
    PyObject *Transposed = PyArray_Transpose(reinterpret_cast<PyArrayObject *>(res), nullptr);
    Py_DECREF(res);
    return Transposed;
  } else
    return res;
}

namespace impl
{
  template <class T>
  struct is_integral_constant : std::false_type {
  };
  template <class T, T N>
  struct is_integral_constant<std::integral_constant<T, N>> : std::true_type {
  };

  template <class pS, class T, size_t... Is>
  bool check_shape(T const *dims, std::index_sequence<Is...>)
  {
    types::array_tuple<bool, sizeof...(Is)> dims_match = {
        (is_integral_constant<std::tuple_element_t<Is, pS>>::value
             ? (dims[Is] ==
                std::conditional_t<is_integral_constant<std::tuple_element_t<Is, pS>>::value,
                                   std::tuple_element_t<Is, pS>,
                                   std::integral_constant<long, 0>>::value)
             : true)...};
    return std::find(dims_match.begin(), dims_match.end(), false) == dims_match.end();
  }

  template <typename T, class pS>
  PyArrayObject *check_array_type_and_dims(PyObject *obj)
  {
    if (!PyArray_Check(obj))
      return nullptr;
    // the array must have the same dtype && the same number of dimensions
    PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
    if (PyArray_TYPE(arr) != c_type_to_numpy_type<T>::value)
      return nullptr;
    if (PyArray_NDIM(arr) != std::tuple_size<pS>::value)
      return nullptr;
    return arr;
  }

  template <class T, class Slice, class S>
  void fill_slice(Slice &slice, long const *strides, long const *offsets, S const *dims,
                  utils::int_<0>)
  {
  }

  template <long stride>
  inline void set_slice(types::cstride_normalized_slice<stride> &cs, long lower, long upper,
                        long step)
  {
    cs.lower = lower;
    cs.upper = upper;
    assert(cs.step == step && "consistent steps");
  }
  inline void set_slice(types::normalized_slice &s, long lower, long upper, long step)
  {
    s.lower = lower;
    s.upper = upper;
    s.step = step;
  }

  template <class T, class Slice, class S, size_t N>
  void fill_slice(Slice &slice, long const *strides, long const *offsets, S const *dims,
                  utils::int_<N>)
  {
    set_slice(std::get<std::tuple_size<Slice>::value - N>(slice), *offsets,
              *offsets + *dims * *strides, *strides);
    fill_slice<T>(slice, strides + 1, offsets + 1, dims + 1, utils::int_<N - 1>());
  }
} // namespace impl

template <typename T, class pS>
bool from_python<types::ndarray<T, pS>>::is_convertible(PyObject *obj)
{
  PyArrayObject *arr = impl::check_array_type_and_dims<T, pS>(obj);
  if (!arr)
    return false;
  auto const *stride = PyArray_STRIDES(arr);
  auto const *dims = PyArray_DIMS(arr);
  long current_stride = PyArray_ITEMSIZE(arr);
  if (PyArray_SIZE(arr)) {
    for (long i = std::tuple_size<pS>::value - 1; i >= 0; i--) {
      if (stride[i] == 0 && dims[i] == 1) {
        // happens when a new dim is added though None/newaxis
      } else if (stride[i] != current_stride && dims[i] > 1) {
        return false;
      }
      current_stride *= dims[i];
    }
    // this is supposed to be a texpr
    if ((PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS) &&
        ((PyArray_FLAGS(arr) & NPY_ARRAY_C_CONTIGUOUS) == 0) && (std::tuple_size<pS>::value > 1)) {
      return false;
    }
  }

  // check if dimension size match
  return impl::check_shape<pS>(dims, std::make_index_sequence<std::tuple_size<pS>::value>());
}
template <typename T, class pS>
types::ndarray<T, pS> from_python<types::ndarray<T, pS>>::convert(PyObject *obj)
{
  PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
  types::ndarray<T, pS> r((T *)PyArray_BYTES(arr), PyArray_DIMS(arr), obj);
  Py_INCREF(obj);
  return r;
}

template <typename T, class pS, class... S>
bool from_python<types::numpy_gexpr<types::ndarray<T, pS>, S...>>::is_convertible(PyObject *obj)
{
  PyArrayObject *arr = impl::check_array_type_and_dims<T, pS>(obj);
  if (!arr)
    return false;

  if ((PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS) &&
      ((PyArray_FLAGS(arr) & NPY_ARRAY_C_CONTIGUOUS) == 0) && (std::tuple_size<pS>::value > 1)) {
    return false;
  }

  PyObject *base_obj = PyArray_BASE(arr);
  if (!base_obj || !PyArray_Check(base_obj))
    return false;
  PyArrayObject *base_arr = reinterpret_cast<PyArrayObject *>(base_obj);

  auto const *stride = PyArray_STRIDES(arr);
  auto const *dims = PyArray_DIMS(arr);

  /* FIXME If we have at least one stride, we convert the whole
   * array to a numpy_gexpr, without trying to be smarter with
   * contiguous slices
   */
  long current_stride = PyArray_ITEMSIZE(arr);
  bool at_least_one_stride = false;
  for (long i = std::tuple_size<pS>::value - 1; i >= 0; i--) {
    if (stride[i] < 0) {
      return false;
    }
    if (stride[i] == 0 && dims[i] == 1) {
      // happens when a new dim is added though None/newaxis
    } else if (stride[i] != current_stride) {
      at_least_one_stride = true;
      break;
    }
    current_stride *= dims[i];
  }
  if (at_least_one_stride) {
    if (PyArray_NDIM(base_arr) != std::tuple_size<pS>::value) {
      return false;
    }
    return true;
  } else
    return false;
}

template <typename T, class pS, class... S>
types::numpy_gexpr<types::ndarray<T, pS>, S...>
from_python<types::numpy_gexpr<types::ndarray<T, pS>, S...>>::convert(PyObject *obj)
{
  PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
  PyArrayObject *base_arr = reinterpret_cast<PyArrayObject *>(PyArray_BASE(arr));

  /* from the base array pointer && this array pointer, we can recover the
   * full slice informations
   * unfortunately, the PyArray representation is different from our.
   * - PyArray_BYTES gives the start of the base pointer
   * - PyArray_Dims give the dimension array (the shape)
   * - PyArray_STRIDES gives the stride information, but relative to the
   * base
   * pointer and not relative to the lower dimension
   */
  long offsets[std::tuple_size<pS>::value];
  long strides[std::tuple_size<pS>::value];
  auto const *base_dims = PyArray_DIMS(base_arr);

  auto full_offset = (PyArray_BYTES(arr) - PyArray_BYTES(base_arr)) / sizeof(T);
  auto const *arr_strides = PyArray_STRIDES(arr);
  long accumulated_dim = 1;
  offsets[std::tuple_size<pS>::value - 1] = full_offset % base_dims[std::tuple_size<pS>::value - 1];
  strides[std::tuple_size<pS>::value - 1] = arr_strides[std::tuple_size<pS>::value - 1] / sizeof(T);
  for (ssize_t i = std::tuple_size<pS>::value - 2; i >= 0; --i) {
    accumulated_dim *= base_dims[i + 1];
    offsets[i] = full_offset / accumulated_dim;
    strides[i] = arr_strides[i] / sizeof(T) / accumulated_dim;
  }
  types::ndarray<T, pS> base_array((T *)PyArray_BYTES(base_arr), PyArray_DIMS(base_arr),
                                   (PyObject *)base_arr);
  std::tuple<S...> slices;
  impl::fill_slice<T>(slices, strides, offsets, PyArray_DIMS(arr), utils::int_<sizeof...(S)>());
  types::numpy_gexpr<types::ndarray<T, pS>, S...> r(base_array, slices);

  Py_INCREF((PyObject *)base_arr);
  return r;
}

template <typename E>
bool from_python<types::numpy_texpr<E>>::

    is_convertible(PyObject *obj)
{
  constexpr auto N = E::value;
  PyArrayObject *arr = impl::check_array_type_and_dims<typename E::dtype, typename E::shape_t>(obj);
  if (!arr)
    return false;
  // check strides. Note that because it's a texpr, the check is done in the
  // opposite direction compared to ndarrays
  auto const *stride = PyArray_STRIDES(arr);
  auto const *dims = PyArray_DIMS(arr);
  long current_stride = PyArray_ITEMSIZE(arr);
  for (size_t i = 0; i < N; i++) {
    if (stride[i] != current_stride)
      return false;
    current_stride *= dims[i];
  }

  return PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS && N > 1;
}

template <typename E>
types::numpy_texpr<E> from_python<types::numpy_texpr<E>>::convert(PyObject *obj)
{
  constexpr size_t N = E::value;
  using T = typename E::dtype;
  PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
  typename E::shape_t shape;
  auto const *dims = PyArray_DIMS(arr);
  static_assert(N == 2, "only support texpr of matrices");
  sutils::assign(std::get<0>(shape), std::get<1>(dims));
  sutils::assign(std::get<1>(shape), std::get<0>(dims));

  PyObject *tobj = PyArray_Transpose(arr, nullptr);
  types::ndarray<T, typename E::shape_t> base_array((T *)PyArray_BYTES(arr), shape, tobj);
  types::numpy_texpr<types::ndarray<T, typename E::shape_t>> r(base_array);
  return r;
}
PYTHONIC_NS_END

#endif

#endif
