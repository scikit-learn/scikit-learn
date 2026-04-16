#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_VEXPR_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_VEXPR_HPP

#include "pythonic/include/types/nditerator.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T, class F>
  struct numpy_vexpr {

    static constexpr size_t value = T::value;
    static const bool is_vectorizable = false;
    static const bool is_flat = false;
    using dtype = typename dtype_of<T>::type;
    using value_type = T;
    static constexpr bool is_strided = T::is_strided;

    using iterator = nditerator<numpy_vexpr>;
    using const_iterator = const_nditerator<numpy_vexpr>;

    T data_;
    F view_;

    numpy_vexpr() = default;

    numpy_vexpr(T const &data, F const &view) : data_(data), view_(view)
    {
    }

    long flat_size() const
    {
      return sutils::prod_tail(data_) * view_.template shape<0>();
    }

    long size() const
    {
      return view_.size();
    }

    template <class E>
    std::enable_if_t<is_iterable<E>::value, numpy_vexpr &> operator=(E const &);
    template <class E>
    std::enable_if_t<!is_iterable<E>::value, numpy_vexpr &> operator=(E const &expr);

    numpy_vexpr &operator=(numpy_vexpr const &);

    using shape_t = array_tuple<long, value>;
    template <size_t I>
    long shape() const
    {
      if (I == 0)
        return view_.template shape<0>();
      else
        return data_.template shape<I>();
    }

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;
#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<numpy_vexpr>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    template <class... Indices>
    dtype load(long i, Indices... indices) const
    {
      return data_.load(view_.fast(i), indices...);
    }
    template <class Elt, class... Indices>
    void store(Elt elt, long i, Indices... indices) const
    {
      data_.store(elt, view_.fast(i), indices...);
    }
    template <class Op, class Elt, class... Indices>
    void update(Elt elt, long i, Indices... indices) const
    {
      data_.template update<Op>(elt, view_.fast(i), indices...);
    }

    auto fast(long i) -> decltype(data_.fast(i))
    {
      return data_.fast(view_.fast(i));
    }
    auto fast(long i) const -> decltype(data_.fast(i))
    {
      return data_.fast(view_.fast(i));
    }
    template <class... S>
    auto operator()(S const &...slices) const
        -> decltype(ndarray<dtype, array_tuple<long, value>>{*this}(slices...));

    auto operator[](long i) const -> decltype(data_[i])
    {
      return data_.fast(view_[i]);
    }

    template <class S>
    std::enable_if_t<is_slice<S>::value,
                     numpy_gexpr<numpy_vexpr, decltype(std::declval<S>().normalize(1))>>
    operator[](S s) const
    {
      return {*this, s.normalize(size())};
    }

    /* element filtering */
    template <class E> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<E>::value && std::is_same<bool, typename E::dtype>::value &&
                         !is_pod_array<F>::value,
                     numpy_vexpr<numpy_vexpr, ndarray<long, pshape<long>>>>
    fast(E const &filter) const;

    template <class E> // indexing through an array of boolean -- a mask
    std::enable_if_t<!is_slice<E>::value && is_numexpr_arg<E>::value &&
                         std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_vexpr, ndarray<long, pshape<long>>>>
    operator[](E const &filter) const;

    template <class E> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<E>::value && !is_array_index<E>::value &&
                         !std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_vexpr, E>>
    operator[](E const &filter) const;

    template <class E> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<E>::value && !is_array_index<E>::value &&
                         !std::is_same<bool, typename E::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_vexpr, E>>
    fast(E const &filter) const;

    template <class Op, class Expr>
    numpy_vexpr &update_(Expr const &expr);

    template <class E>
    numpy_vexpr &operator+=(E const &expr);

    template <class E>
    numpy_vexpr &operator-=(E const &expr);

    template <class E>
    numpy_vexpr &operator*=(E const &expr);

    template <class E>
    numpy_vexpr &operator/=(E const &expr);

    template <class E>
    numpy_vexpr &operator&=(E const &expr);

    template <class E>
    numpy_vexpr &operator|=(E const &expr);

    template <class E>
    numpy_vexpr &operator^=(E const &expr);
  };
} // namespace types

template <class T, class F>
struct assignable<types::numpy_vexpr<T, F>> {
  using type =
      types::ndarray<typename types::dtype_of<T>::type, typename types::numpy_vexpr<T, F>::shape_t>;
};

template <class T, class F>
struct lazy<types::numpy_vexpr<T, F>> {
  using type = types::numpy_vexpr<typename lazy<T>::type, typename lazy<F>::type>;
};

PYTHONIC_NS_END

/* combined are sorted such that the assigned type comes first */
template <class E, class F, class T, class pS>
struct __combined<pythonic::types::numpy_vexpr<E, F>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class E, class F, class T, class pS>
struct __combined<pythonic::types::ndarray<T, pS>, pythonic::types::numpy_vexpr<E, F>> {
  using type = pythonic::types::ndarray<T, pS>;
};

#endif
