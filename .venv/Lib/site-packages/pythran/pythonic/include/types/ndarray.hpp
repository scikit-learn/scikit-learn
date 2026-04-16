#ifndef PYTHONIC_INCLUDE_TYPES_NDARRAY_HPP
#define PYTHONIC_INCLUDE_TYPES_NDARRAY_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/attr.hpp"
#include "pythonic/include/types/empty_iterator.hpp"

#include "pythonic/include/utils/broadcast_copy.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/utils/reserve.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/raw_array.hpp"
#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/types/tuple.hpp"

#include "pythonic/include/numpy/bool_.hpp"
#include "pythonic/include/numpy/complex128.hpp"
#include "pythonic/include/numpy/complex256.hpp"
#include "pythonic/include/numpy/complex64.hpp"
#include "pythonic/include/numpy/float128.hpp"
#include "pythonic/include/numpy/float32.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/numpy/int16.hpp"
#include "pythonic/include/numpy/int32.hpp"
#include "pythonic/include/numpy/int64.hpp"
#include "pythonic/include/numpy/int8.hpp"
#include "pythonic/include/numpy/uint16.hpp"
#include "pythonic/include/numpy/uint32.hpp"
#include "pythonic/include/numpy/uint64.hpp"
#include "pythonic/include/numpy/uint8.hpp"

#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/types/numpy_expr.hpp"
#include "pythonic/include/types/numpy_gexpr.hpp"
#include "pythonic/include/types/numpy_iexpr.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/types/numpy_texpr.hpp"
#include "pythonic/include/types/numpy_vexpr.hpp"
#include "pythonic/include/types/pointer.hpp"
#include "pythonic/include/types/vectorizable_type.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include "pythonic/include/builtins/len.hpp"

#include <array>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <ostream>

#ifdef ENABLE_PYTHON_MODULE
// Cython 0.29.x still uses the deprecated API, so we can't set this macro in
// this case! Also avoid redefining it if already set by the Pythran user.
#if !defined(CYTHON_ABI) && !defined(NPY_NO_DEPRECATED_API)
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include "numpy/arrayobject.h"
#endif

#ifdef USE_XSIMD
#include <xsimd/xsimd.hpp>
#endif

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T, class pS>
  struct ndarray;

  template <class T>
  struct type_helper;

  /* Helper for dimension-specific part of ndarray
   *
   * Instead of specializing the whole ndarray class, the dimension-specific
   *behavior are stored here.
   * There are two specialization for this type:
   * - a specialization depending on the dimensionality (==1 || > 1)
   * - a specialization depending on the constness.
   *
   * The raw ndarray<T,pS> specialization implies a *swallow copy* of the
   *ndarray, && thus a refcount increase.
   * It is meant to be used when indexing an rvalue, as in
   *``np.zeros(10)[i]``.
   *
   * The ndarray<T,pS> const& specialization implies a *reference copy*. It is
   *used when indexing a lvalue, as in ``a[i]``
   */

  template <class T, class pS>
  struct type_helper<ndarray<T, pS>> {
    static_assert(std::tuple_size<pS>::value != 1, "matching ok");
    using type = numpy_iexpr<ndarray<T, pS>>;
    using iterator = nditerator<ndarray<T, pS>>;
    using const_iterator = const_nditerator<ndarray<T, pS>>;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, pS> &n, long i);
    static const_iterator make_iterator(ndarray<T, pS> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);

    static numpy_iexpr<ndarray<T, pS>> get(ndarray<T, pS> &&self, long i);
  };

  template <class T, class pS>
  struct type_helper<ndarray<T, pS> const &> {
    static_assert(std::tuple_size<pS>::value != 1, "matching ok");
    using type = numpy_iexpr<ndarray<T, pS> const &>;

    using iterator = nditerator<ndarray<T, pS>>;
    using const_iterator = const_nditerator<ndarray<T, pS>>;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, pS> &n, long i);
    static const_iterator make_iterator(ndarray<T, pS> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);

    static numpy_iexpr<ndarray<T, pS> const &> get(ndarray<T, pS> const &self, long i);
  };

  template <class T, class pS>
  struct type_helper<ndarray<T, pshape<pS>>> {
    using type = T;

    using iterator = T *;
    using const_iterator = T const *;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, pshape<pS>> &n, long i);
    static const_iterator make_iterator(ndarray<T, pshape<pS>> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);

    static type get(ndarray<T, pshape<pS>> &&self, long i);
  };

  template <class T, class pS>
  struct type_helper<ndarray<T, pshape<pS>> const &> {
    using type = T;

    using iterator = T *;
    using const_iterator = T const *;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, pshape<pS>> &n, long i);
    static const_iterator make_iterator(ndarray<T, pshape<pS>> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);
    static type &get(ndarray<T, pshape<pS>> const &self, long i);
  };

  template <class T, class pS>
  struct type_helper<ndarray<T, array_tuple<pS, 1>>> {
    using type = T;

    using iterator = T *;
    using const_iterator = T const *;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, array_tuple<pS, 1>> &n, long i);
    static const_iterator make_iterator(ndarray<T, array_tuple<pS, 1>> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);

    static type get(ndarray<T, array_tuple<pS, 1>> &&self, long i);
  };

  template <class T, class pS>
  struct type_helper<ndarray<T, array_tuple<pS, 1>> const &> {
    using type = T;

    using iterator = T *;
    using const_iterator = T const *;

    type_helper() = delete; // Not intended to be instantiated

    static iterator make_iterator(ndarray<T, array_tuple<pS, 1>> &n, long i);
    static const_iterator make_iterator(ndarray<T, array_tuple<pS, 1>> const &n, long i);

    template <class S, class Iter>
    static T *initialize_from_iterable(S &shape, T *from, Iter &&iter);
    static type &get(ndarray<T, array_tuple<pS, 1>> const &self, long i);
  };

  /* Multidimensional array of values
   *
   * An ndarray wraps a raw array pointers && manages multiple dimensions
   * casted overt the raw data.
   * The number of dimensions is fixed as well as the type of the underlying
   * data.
   * A shared pointer is used internally to mimic Python's behavior.
   *
   */
  template <class T, class pS>
  struct ndarray {
    static const bool is_vectorizable = types::is_vectorizable<T>::value;
    static const bool is_flat = true;
    static const bool is_strided = false;

    /* types */
    static constexpr size_t value = std::tuple_size<pS>::value;
    using dtype = T;
    using value_type = typename type_helper<ndarray>::type;
    using reference = value_type &;
    using const_reference = value_type const &;

    using iterator = typename type_helper<ndarray>::iterator;
    using const_iterator = typename type_helper<ndarray>::const_iterator;
    using flat_iterator = T *;
    using const_flat_iterator = T const *;

    using shape_t = pS;
    static_assert(std::tuple_size<shape_t>::value == value, "consistent shape size");

    /* members */
    utils::shared_ref<raw_array<T>> mem; // shared data pointer
    T *buffer;                           // pointer to the first data stored in the equivalent flat
                                         // array
    shape_t _shape;                      // shape of the multidimensional array
    sutils::concat_t<types::array_tuple<long, value - 1>, pshape<std::integral_constant<long, 1>>>
        _strides; // strides

    /* mem management */
    void mark_memory_external(extern_type obj)
    {
      mem.external(obj);
      mem->forget();
    }

    /* constructors */
    ndarray();
    ndarray(ndarray const &) = default;
    ndarray(ndarray &&) = default;

    /* assignment */
    ndarray &operator=(ndarray const &other) = default;

    /* from other memory */
    ndarray(utils::shared_ref<raw_array<T>> const &mem, pS const &shape);
    ndarray(utils::shared_ref<raw_array<T>> &&mem, pS const &shape);

    /* from other array */
    template <class Tp, class pSp>
    ndarray(ndarray<Tp, pSp> const &other);
    template <class pSp>
    ndarray(ndarray<T, pSp> const &other);

    /* from a seed */
    ndarray(pS const &shape, none_type init);
    ndarray(pS const &shape, T init);

    /* from a foreign pointer */
    template <class S>
    ndarray(T *data, S const *pshape, ownership o);
    ndarray(T *data, pS const &pshape, ownership o);

#ifdef ENABLE_PYTHON_MODULE
    template <class S>
    ndarray(T *data, S const *pshape, PyObject *obj);
    ndarray(T *data, pS const &pshape, PyObject *obj);
#endif

    template <class Iterable,
              class = std::enable_if_t<
                  !is_array<std::remove_cv_t<std::remove_reference_t<Iterable>>>::value &&
                      is_iterable<std::remove_cv_t<std::remove_reference_t<Iterable>>>::value,
                  void>>
    ndarray(Iterable &&iterable);

    /* from a  numpy expression */
    template <class E>
    void initialize_from_expr(E const &expr);

    template <class Op, class... Args>
    ndarray(numpy_expr<Op, Args...> const &expr);

    template <class Arg>
    ndarray(numpy_texpr<Arg> const &expr);

    template <class Arg>
    ndarray(numpy_texpr_2<Arg> const &expr);

    template <class Arg, class... S>
    ndarray(numpy_gexpr<Arg, S...> const &expr);

    template <class Arg>
    ndarray(numpy_iexpr<Arg> const &expr);

    template <class Arg, class F>
    ndarray(numpy_vexpr<Arg, F> const &expr);

    /* update operators */
    template <class Op, class Expr>
    ndarray &update_(Expr const &expr);
    template <class Expr>
    ndarray &operator+=(Expr const &expr);

    template <class Expr>
    ndarray &operator-=(Expr const &expr);

    template <class Expr>
    ndarray &operator*=(Expr const &expr);

    template <class Expr>
    ndarray &operator/=(Expr const &expr);

    template <class Expr>
    ndarray &operator&=(Expr const &expr);

    template <class Expr>
    ndarray &operator|=(Expr const &expr);

    template <class Expr>
    ndarray &operator^=(Expr const &expr);

    template <class E, class... Indices>
    void store(E elt, Indices... indices)
    {
      static_assert(is_dtype<E>::value, "valid store");
      *(buffer + noffset<std::tuple_size<pS>::value>{}(
                     *this, array_tuple<long, value>{{indices...}})) = static_cast<E>(elt);
    }
    template <class... Indices>
    dtype load(Indices... indices) const
    {
      return *(buffer + noffset<std::tuple_size<pS>::value>{}(
                            *this, array_tuple<long, value>{{indices...}}));
    }

    template <class Op, class E, class... Indices>
    void update(E elt, Indices... indices) const
    {
      static_assert(is_dtype<E>::value, "valid store");
      Op{}(*(buffer +
             noffset<std::tuple_size<pS>::value>{}(*this, array_tuple<long, value>{{indices...}})),
           static_cast<E>(elt));
    }

    /* element indexing
     * differentiate const from non const, && r-value from l-value
     * */
    auto fast(long i) const & -> decltype(type_helper<ndarray const &>::get(*this, i))
    {
      return type_helper<ndarray const &>::get(*this, i);
    }

    auto fast(long i) && -> decltype(type_helper<ndarray>::get(std::move(*this), i))
    {
      return type_helper<ndarray>::get(std::move(*this), i);
    }

    template <class Ty>
    std::enable_if_t<std::is_integral<Ty>::value, T &> fast(array_tuple<Ty, value> const &indices);
    template <class Ty>
    std::enable_if_t<std::is_integral<Ty>::value, T>
    fast(array_tuple<Ty, value> const &indices) const;

    template <class Ty, size_t M>
    auto fast(array_tuple<Ty, M> const &indices) const & -> std::enable_if_t<
        std::is_integral<Ty>::value, decltype(nget<M - 1>().fast(*this, indices))>;

    template <class Ty, size_t M>
    auto fast(array_tuple<Ty, M> const &indices) && -> std::enable_if_t<
        std::is_integral<Ty>::value, decltype(nget<M - 1>().fast(std::move(*this), indices))>;

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<ndarray>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

#ifndef NDEBUG
    template <class IndicesTy>
    bool inbound_indices(IndicesTy const &indices) const
    {
      auto const shp = sutils::getshape(*this);
      for (size_t i = 0, n = indices.size(); i < n; ++i) {
        auto const index = indices[i];
        auto const dim = shp[i];
        if (0 > index || index >= dim)
          return false;
      }
      return true;
    }
#endif

    /* slice indexing */
    ndarray<T, sutils::push_front_t<pS, std::integral_constant<long, 1>>>
    operator[](none_type) const;

    template <class S>
    std::enable_if_t<is_slice<S>::value, numpy_gexpr<ndarray const &, normalize_t<S>>>
    operator[](S const &s) const &;

    template <class S>
    std::enable_if_t<is_slice<S>::value, numpy_gexpr<ndarray, normalize_t<S>>>
    operator[](S const &s) &&;

    long size() const;

    /* extended slice indexing */
    template <class Ty>
    auto operator()(Ty s) const
        -> std::enable_if_t<std::is_integral<Ty>::value, decltype((*this)[s])>
    {
      return (*this)[s];
    }

    template <class S0, class... S>
    auto operator()(S0 const &s0, S const &...s)
        const & -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...));

    template <class S0, class... S>
    auto operator()(S0 const &s0, S const &...s)
        & -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}((*this), s0, s...));

    template <class S0, class... S>
    auto operator()(S0 const &s0, S const &...s)
        && -> decltype(extended_slice<count_new_axis<S0, S...>::value>{}(std::move(*this), s0,
                                                                         s...));

    /* element filtering */
    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value == 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value == 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value != 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray<T, pshape<long>>, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value != 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray<T, pshape<long>>, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray, F>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray, F>>
    fast(F const &filter) const;

    auto operator[](long i) const & -> decltype(this->fast(i))
    {
      if (i < 0)
        i += std::get<0>(_shape);
      assert(0 <= i && i < std::get<0>(_shape));
      return fast(i);
    }

    auto operator[](long i) && -> decltype(std::move(*this).fast(i))
    {
      if (i < 0)
        i += std::get<0>(_shape);
      assert(0 <= i && i < std::get<0>(_shape));
      return std::move(*this).fast(i);
    }

    template <class Ty>
    std::enable_if_t<std::is_integral<Ty>::value, T const &>
    operator[](array_tuple<Ty, value> const &indices) const;

    template <class Ty>
    std::enable_if_t<std::is_integral<Ty>::value, T &>
    operator[](array_tuple<Ty, value> const &indices);

    template <class Ty, size_t M>
    auto operator[](array_tuple<Ty, M> const &indices) const
        & -> std::enable_if_t<std::is_integral<Ty>::value, decltype(nget<M - 1>()(*this, indices))>;

    template <class Ty, size_t M>
    auto operator[](array_tuple<Ty, M> const &indices) && -> std::enable_if_t<
        std::is_integral<Ty>::value, decltype(nget<M - 1>()(std::move(*this), indices))>;

    template <class... Tys, size_t... Is>
    auto _fwdlongindex(std::tuple<Tys...> const &indices, std::index_sequence<Is...>) const
        -> decltype((*this)(static_cast<long>(Is)...))
    {
      return (*this)(static_cast<long>(std::get<Is>(indices))...);
    }

    template <class... Tys>
    auto operator[](std::tuple<Tys...> const &indices) const -> std::enable_if_t<
        utils::all_of<std::is_integral<Tys>::value...>::value,
        decltype(this->_fwdlongindex(indices, std::make_index_sequence<sizeof...(Tys)>()))>
    {
      return _fwdlongindex(indices, std::make_index_sequence<sizeof...(Tys)>());
    }

    template <class... Tys, size_t... Is>
    auto _fwdlongindex(std::tuple<Tys...> const &indices, std::index_sequence<Is...>)
        -> decltype((*this)(static_cast<long>(Is)...))
    {
      return (*this)(static_cast<long>(std::get<Is>(indices))...);
    }

    template <class... Tys>
    auto operator[](std::tuple<Tys...> const &indices) -> std::enable_if_t<
        utils::all_of<std::is_integral<Tys>::value...>::value,
        decltype(this->_fwdlongindex(indices, std::make_index_sequence<sizeof...(Tys)>()))>
    {
      return _fwdlongindex(indices, std::make_index_sequence<sizeof...(Tys)>());
    }

    template <class Ty0, class Ty1, class... Tys>
    auto operator[](std::tuple<Ty0, Ty1, Tys...> const &indices) const -> std::enable_if_t<
        std::is_integral<Ty0>::value &&
            !utils::all_of<std::is_integral<Ty1>::value, std::is_integral<Tys>::value...>::value,
        decltype((*this)[std::get<0>(indices)][tuple_tail(indices)])>
    {
      return (*this)[std::get<0>(indices)][tuple_tail(indices)];
    }

    template <class Slices, size_t... Is>
    auto
    _fwdindex(Slices const &indices,
              std::index_sequence<Is...>) const & -> decltype((*this)(std::get<Is>(indices)...))
    {
      return (*this)(std::get<Is>(indices)...);
    }

    template <class S, size_t... Is>
    auto
    _fwdindex(dynamic_tuple<S> const &indices,
              std::index_sequence<Is...>) const & -> decltype((*this)(std::get<Is>(indices)...))
    {
      return (*this)((indices.size() > Is ? std::get<Is>(indices) : cstride_slice<1>())...);
    }

    template <class Ty0, class Ty1, class... Tys,
              class _ = std::enable_if_t<is_numexpr_arg<Ty0>::value, void>>
    auto operator[](std::tuple<Ty0, Ty1, Tys...> const &indices) const -> std::enable_if_t<
        is_numexpr_arg<Ty0>::value,
        decltype(this->_fwdindex(indices, std::make_index_sequence<2 + sizeof...(Tys)>()))>;

    template <class Ty, size_t M, class _ = std::enable_if_t<!std::is_integral<Ty>::value, void>>
    auto operator[](array_tuple<Ty, M> const &indices)
        const & -> decltype(this->_fwdindex(indices, std::make_index_sequence<M>()))
    {
      return _fwdindex(indices, std::make_index_sequence<M>());
    }
    template <class S>
    auto operator[](dynamic_tuple<S> const &indices) const
        -> decltype(this->_fwdindex(indices, std::make_index_sequence<value>()))
    {
      return _fwdindex(indices, std::make_index_sequence<value>());
    }

    /* through iterators */
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

    const_flat_iterator fbegin() const;
    const_flat_iterator fend() const;
    flat_iterator fbegin();
    flat_iterator fend();

    /* member functions */
    T *data()
    {
      return buffer;
    }
    T const *data() const
    {
      return buffer;
    }
    long flat_size() const;
    bool may_overlap(ndarray const &) const;

    template <class qS>
    ndarray<T, qS> reshape(qS const &shape) const &;
    template <class qS>
    ndarray<T, qS> reshape(qS const &shape) &&;

    template <class OT>
    ndarray<OT, types::array_tuple<long, value>> recast()
    {
      auto new_shape = sutils::array(_shape);
      new_shape[value - 1] = new_shape[value - 1] * sizeof(T) / sizeof(OT);
      auto new_mem = mem.template recast<raw_array<OT>>();
      return ndarray<OT, types::array_tuple<long, value>>(new_mem, new_shape);
    }

    explicit operator bool() const;

    ndarray<T, pshape<long>> flat() const;
    ndarray<T, pS> copy() const;
    intptr_t id() const;
    intptr_t baseid() const
    {
      return id();
    }
    template <size_t I>
    auto shape() const -> decltype(std::get<I>(_shape))
    {
      return std::get<I>(_shape);
    }

    template <size_t I>
    auto strides() const -> decltype(std::get<I>(_strides))
    {
      return std::get<I>(_strides);
    }

    operator pointer<T>()
    {
      return {buffer};
    }
  };

  /* pretty printing { */
  template <class T, class pS>
  std::ostream &operator<<(std::ostream &os, ndarray<T, pS> const &e);

  template <class E>
  std::enable_if_t<is_array<E>::value, std::ostream &> operator<<(std::ostream &os, E const &e);
  /* } */
} // namespace types
PYTHONIC_NS_END

/* std::get overloads */
namespace std
{

  template <size_t I, class E>
  auto get(E &&a) -> std::enable_if_t<
      pythonic::types::is_array<std::remove_cv_t<std::remove_reference_t<E>>>::value,
      decltype(std::forward<E>(a)[I])>;

  template <size_t I, class T, class pS>
  struct tuple_element<I, pythonic::types::ndarray<T, pS>> {
    using type = typename pythonic::types::ndarray<T, pS>::value_type;
  };

  template <size_t I, class Op, class... Args>
  struct tuple_element<I, pythonic::types::numpy_expr<Op, Args...>> {
    using type = typename pythonic::types::numpy_expr<Op, Args...>::dtype;
  };

  template <size_t I, class E>
  struct tuple_element<I, pythonic::types::numpy_iexpr<E>> {
    using type = decltype(std::declval<pythonic::types::numpy_iexpr<E>>()[0]);
  };

  template <size_t I, class E>
  struct tuple_element<I, pythonic::types::numpy_texpr<E>> {
    using type = decltype(std::declval<pythonic::types::numpy_texpr<E>>()[0]);
  };

  template <size_t I, class E, class... S>
  struct tuple_element<I, pythonic::types::numpy_gexpr<E, S...>> {
    using type = decltype(std::declval<pythonic::types::numpy_gexpr<E, S...>>()[0]);
  };

  template <size_t I, class T, class F>
  struct tuple_element<I, pythonic::types::numpy_vexpr<T, F>> {
    using type = decltype(std::declval<pythonic::types::numpy_vexpr<T, F>>()[0]);
  };

} // namespace std

/* pythran attribute system { */
#include "pythonic/include/numpy/transpose.hpp"
PYTHONIC_NS_BEGIN

namespace types
{

  namespace details
  {
    using dtype_table =
        std::tuple<void, pythonic::numpy::functor::int8, pythonic::numpy::functor::int16, void,
                   pythonic::numpy::functor::int32, void, void, void,
                   pythonic::numpy::functor::int64>;
    using dtype_utable =
        std::tuple<void, pythonic::numpy::functor::uint8, pythonic::numpy::functor::uint16, void,
                   pythonic::numpy::functor::uint32, void, void, void,
                   pythonic::numpy::functor::uint64>;

    template <class T>
    struct dtype_helper {
      using table = std::conditional_t<std::is_signed<T>::value, dtype_table, dtype_utable>;
      using type =
          std::tuple_element_t<(sizeof(T) < std::tuple_size<table>::value) ? sizeof(T) : 0, table>;
    };

    template <>
    struct dtype_helper<bool> {
      using type = pythonic::numpy::functor::bool_;
    };

    template <>
    struct dtype_helper<float> {
      using type = pythonic::numpy::functor::float32;
    };
    template <>
    struct dtype_helper<double> {
      using type = pythonic::numpy::functor::float64;
    };
    template <>
    struct dtype_helper<long double> {
      using type = pythonic::numpy::functor::float128;
    };
    template <>
    struct dtype_helper<std::complex<float>> {
      using type = pythonic::numpy::functor::complex64;
    };
    template <>
    struct dtype_helper<std::complex<double>> {
      using type = pythonic::numpy::functor::complex128;
    };
    template <>
    struct dtype_helper<std::complex<long double>> {
      using type = pythonic::numpy::functor::complex256;
    };
  } // namespace details
  template <class T>
  using dtype_t = typename details::dtype_helper<T>::type;
} // namespace types
namespace builtins
{
  namespace details
  {
    template <size_t N>
    struct _build_gexpr {
      template <class E, class... S>
      auto operator()(E const &a, S const &...slices)
          -> decltype(_build_gexpr<N - 1>{}(a, types::cstride_slice<1>(), slices...));
    };

    template <>
    struct _build_gexpr<1> {
      template <class E, class... S>
      auto operator()(E const &a, S const &...slices) -> decltype(E(a)(slices...));
    };

    template <class E>
    E _make_real(E const &a, utils::int_<0>);

    template <class E>
    auto _make_real(E const &a, utils::int_<1>) -> decltype(_build_gexpr<E::value>{}(
        types::ndarray<typename types::is_complex<typename E::dtype>::type,
                       types::array_tuple<long, E::value>>{},
        types::slice()));
    template <class T, class Ss, size_t... Is>
    auto real_get(T &&expr, Ss const &indices, std::index_sequence<Is...>)
        -> decltype(std::forward<T>(expr)(std::get<Is>(indices)...))
    {
      return std::forward<T>(expr)(std::get<Is>(indices)...);
    }
    template <class E>
    types::ndarray<typename E::dtype, typename E::shape_t> _make_imag(E const &a, utils::int_<0>);

    template <class E>
    auto _make_imag(E const &a, utils::int_<1>) -> decltype(_build_gexpr<E::value>{}(
        types::ndarray<typename types::is_complex<typename E::dtype>::type,
                       types::array_tuple<long, E::value>>{},
        types::slice()));
    template <class T, class Ss, size_t... Is>
    auto imag_get(T &&expr, Ss const &indices, std::index_sequence<Is...>)
        -> decltype(std::forward<T>(expr)(std::get<Is>(indices)...))
    {
      return std::forward<T>(expr)(std::get<Is>(indices)...);
    }
  } // namespace details

  template <class E>
  types::array_tuple<long, E::value> getattr(types::attr::SHAPE, E const &a);

  inline types::pshape<> getattr(types::attr::SHAPE, ...)
  {
    return {};
  }

  template <class E>
  constexpr decltype(long(E::value)) getattr(types::attr::NDIM, E const &a);

  template <class E>
  std::integral_constant<long, E::value> getattr(types::attr::NDIM, E *const &a);

  inline long getattr(types::attr::NDIM, ...)
  {
    return 0;
  }

  template <class E>
  types::array_tuple<long, E::value> getattr(types::attr::STRIDES, E const &a);

  inline std::tuple<> getattr(types::attr::STRIDES, ...)
  {
    return {};
  }

  template <class E>
  long getattr(types::attr::SIZE, E const &a);

  template <class E>
  constexpr long getattr(types::attr::ITEMSIZE, E const &a);

  template <class E>
  std::integral_constant<long, sizeof(typename E::dtype)> getattr(types::attr::ITEMSIZE,
                                                                  E *const &a);

  template <class E>
  long getattr(types::attr::NBYTES, E const &a);

  template <class E>
  auto getattr(types::attr::FLAT, E const &a) -> decltype(a.flat());

  template <class E>
  auto getattr(types::attr::T, E const &a) -> decltype(numpy::transpose(a))
  {
    return numpy::transpose(a);
  }

  template <class T, class pS>
  auto getattr(types::attr::REAL, types::ndarray<T, pS> const &a)
      -> decltype(details::_make_real(a, utils::int_<types::is_complex<T>::value>{}));

  template <class Op, class... Args>
  auto getattr(types::attr::REAL, types::numpy_expr<Op, Args...> const &a)
      -> decltype(details::_make_real(
          a,
          utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{}));

  template <class E>
  auto getattr(types::attr::REAL, types::numpy_texpr<E> const &a)
      -> decltype(types::numpy_texpr<decltype(getattr(types::attr::REAL{}, a.arg))>{
          getattr(types::attr::REAL{}, a.arg)});

  template <class E>
  auto getattr(types::attr::REAL, types::numpy_iexpr<E> const &a)
      -> decltype(types::numpy_iexpr<decltype(getattr(types::attr::REAL{}, a.arg))>{
          getattr(types::attr::REAL{}, a.arg)})
  {
    return {getattr(types::attr::REAL{}, a.arg)};
  }

  template <class T, class F>
  auto getattr(types::attr::REAL, types::numpy_vexpr<T, F> const &a)
      -> decltype(types::numpy_vexpr<decltype(getattr(types::attr::REAL{}, a.data_)), F>{
          getattr(types::attr::REAL{}, a.data_), a.view_})
  {
    return {getattr(types::attr::REAL{}, a.data_), a.view_};
  }

  template <class E, class... S>
  auto getattr(types::attr::REAL, types::numpy_gexpr<E, S...> const &a)
      -> decltype(details::real_get(
          getattr(types::attr::REAL{}, a.arg), a.slices,
          std::make_index_sequence<std::tuple_size<decltype(a.slices)>::value>()))
  {
    return details::real_get(
        getattr(types::attr::REAL{}, a.arg), a.slices,
        std::make_index_sequence<std::tuple_size<decltype(a.slices)>::value>());
  }

  template <class T, class pS>
  auto getattr(types::attr::IMAG, types::ndarray<T, pS> const &a)
      -> decltype(details::_make_imag(a, utils::int_<types::is_complex<T>::value>{}));

  template <class Op, class... Args>
  auto getattr(types::attr::IMAG, types::numpy_expr<Op, Args...> const &a)
      -> decltype(details::_make_imag(
          a,
          utils::int_<types::is_complex<typename types::numpy_expr<Op, Args...>::dtype>::value>{}));

  template <class E>
  auto getattr(types::attr::IMAG, types::numpy_texpr<E> const &a)
      -> decltype(types::numpy_texpr<decltype(getattr(types::attr::IMAG{}, a.arg))>{
          getattr(types::attr::IMAG{}, a.arg)});

  template <class E>
  auto geatttr(types::attr::IMAG, types::numpy_iexpr<E> const &a)
      -> decltype(types::numpy_iexpr<decltype(getattr(types::attr::IMAG{}, a.arg))>{
          getattr(types::attr::IMAG{}, a.arg)})
  {
    return {getattr(types::attr::IMAG{}, a.arg)};
  }

  template <class T, class F>
  auto getattr(types::attr::IMAG, types::numpy_vexpr<T, F> const &a)
      -> decltype(types::numpy_vexpr<decltype(getattr(types::attr::IMAG{}, a.data_)), F>{
          getattr(types::attr::IMAG{}, a.data_), a.view_})
  {
    return {getattr(types::attr::IMAG{}, a.data_), a.view_};
  }

  template <class E, class... S>
  auto getattr(types::attr::IMAG, types::numpy_gexpr<E, S...> const &a)
      -> decltype(details::imag_get(
          getattr(types::attr::IMAG{}, a.arg), a.slices,
          std::make_index_sequence<std::tuple_size<decltype(a.slices)>::value>()))
  {
    return details::imag_get(
        getattr(types::attr::IMAG{}, a.arg), a.slices,
        std::make_index_sequence<std::tuple_size<decltype(a.slices)>::value>());
  }

  template <class E>
  types::dtype_t<typename types::dtype_of<E>::type> getattr(types::attr::DTYPE, E const &);
} // namespace builtins
PYTHONIC_NS_END

/* } */

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class T1, class T2, class pS1, class pS2>
struct __combined<pythonic::types::ndarray<T1, pS1>, pythonic::types::ndarray<T2, pS2>> {
  using type = pythonic::types::ndarray<
      typename __combined<T1, T2>::type,
      pythonic::sutils::common_shapes_t<std::tuple_size<pS1>::value, pS1, pS2>>;
};

template <class pS, class T, class... Tys>
struct __combined<pythonic::types::ndarray<T, pS>, pythonic::types::numpy_expr<Tys...>> {
  using expr_type = pythonic::types::numpy_expr<Tys...>;
  using type =
      pythonic::types::ndarray<typename __combined<T, typename expr_type::dtype>::type,
                               pythonic::sutils::common_shapes_t<std::tuple_size<pS>::value, pS,
                                                                 typename expr_type::shape_t>>;
};

template <class pS, class T, class O>
struct __combined<pythonic::types::ndarray<T, pS>, O> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class pS, class T, class O>
struct __combined<pythonic::types::ndarray<T, pS>, pythonic::types::none<O>> {
  using type = pythonic::types::none<typename __combined<pythonic::types::ndarray<T, pS>, O>::type>;
};

template <class pS, class T, class O>
struct __combined<pythonic::types::none<O>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::none<typename __combined<O, pythonic::types::ndarray<T, pS>>::type>;
};

template <class pS, class T>
struct __combined<pythonic::types::ndarray<T, pS>, pythonic::types::none_type> {
  using type = pythonic::types::none<pythonic::types::ndarray<T, pS>>;
};

template <class pS, class T>
struct __combined<pythonic::types::none_type, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::none<pythonic::types::ndarray<T, pS>>;
};

template <class pS, class T, class C, class I>
struct __combined<indexable_container<C, I>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class pS, class T, class C>
struct __combined<indexable<C>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class pS, class T, class C>
struct __combined<container<C>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

/* } */

#include "pythonic/include/types/numpy_operators.hpp"
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <class T, class pS>
struct to_python<types::ndarray<T, pS>> {
  static PyObject *convert(types::ndarray<T, pS> const &n, bool transpose = false);
};

template <class Arg>
struct to_python<types::numpy_iexpr<Arg>> {
  static PyObject *convert(types::numpy_iexpr<Arg> const &v, bool transpose = false);
};

template <class Arg, class... S>
struct to_python<types::numpy_gexpr<Arg, S...>> {
  static PyObject *convert(types::numpy_gexpr<Arg, S...> const &v, bool transpose = false);
};

template <class E>
struct to_python<types::numpy_texpr<E>> {
  static PyObject *convert(types::numpy_texpr<E> const &t, bool transpose = false)
  {
    auto const &n = t.arg;
    PyObject *result = to_python<E>::convert(n, !transpose);
    return result;
  }
};

template <typename T, class pS>
struct from_python<types::ndarray<T, pS>> {
  static bool is_convertible(PyObject *obj);
  static types::ndarray<T, pS> convert(PyObject *obj);
};

template <typename T, class pS, class... S>
struct from_python<types::numpy_gexpr<types::ndarray<T, pS>, S...>> {
  static bool is_convertible(PyObject *obj);

  static types::numpy_gexpr<types::ndarray<T, pS>, S...> convert(PyObject *obj);
};
template <typename T, class pS, class... S>
struct from_python<types::numpy_gexpr<types::ndarray<T, pS> const &, S...>>
    : from_python<types::numpy_gexpr<types::ndarray<T, pS>, S...>> {
};

template <typename E>
struct from_python<types::numpy_texpr<E>> {

  static bool is_convertible(PyObject *obj);

  static types::numpy_texpr<E> convert(PyObject *obj);
};
PYTHONIC_NS_END

/* specialization of std::copy to avoid the multiple calls implied by the
 * recursive calls to std::copy */
namespace std
{
  template <class T, class pS>
  typename pythonic::types::nditerator<pythonic::types::ndarray<T, pS>>
  copy(typename pythonic::types::const_nditerator<pythonic::types::ndarray<T, pS>> begin,
       typename pythonic::types::const_nditerator<pythonic::types::ndarray<T, pS>> end,
       typename pythonic::types::nditerator<pythonic::types::ndarray<T, pS>> out)
  {
    const long offset = pythonic::sutils::prod_tail(begin.data);
    std::copy(begin.data.buffer + begin.index * offset, end.data.buffer + end.index * offset,
              out.data.buffer + out.index * offset);
    return out + (end - begin);
  }
} // namespace std

#endif

#endif
