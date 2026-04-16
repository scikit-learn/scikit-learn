#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_TEXPR_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_TEXPR_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/transpose.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  template <class Arg, class... S>
  struct numpy_gexpr;

  /* expression template for Transposed matrix */
  template <class Arg>
  struct numpy_texpr;

  // wrapper around numpy expression for 2D transposed matrix using gexpr
  // representation
  // >>> b = a.transpose
  // >>> b[i] == a[:,i]
  // True
  //
  // for N = 2
  template <class E>
  struct numpy_texpr_2 {
    static_assert(E::value == 2, "texpr only implemented for matrices");
    static const bool is_vectorizable = false;
    static const bool is_flat = false;
    static const bool is_strided = true;
    using Arg = E;

    using iterator = nditerator<numpy_texpr_2<Arg>>;
    using const_iterator = const_nditerator<numpy_texpr_2<Arg>>;

    static constexpr size_t value = Arg::value;
    using value_type = numpy_gexpr<Arg, cstride_normalized_slice<1>, long>;
    using dtype = typename E::dtype;

    Arg arg;
    using shape_t = sutils::transpose_t<typename E::shape_t>;
    template <size_t I>
    auto shape() const -> decltype(arg.template shape<I == 0 ? 1 : 0>())
    {
      return arg.template shape<I == 0 ? 1 : 0>();
    }

    numpy_texpr_2();
    numpy_texpr_2(numpy_texpr_2 const &) = default;
    numpy_texpr_2(numpy_texpr_2 &&) = default;
    numpy_texpr_2 &operator=(numpy_texpr_2 const &) = default;

    numpy_texpr_2(Arg const &arg);
    const_iterator begin() const;
    const_iterator end() const;

    iterator begin();
    iterator end();

    long size() const
    {
      return this->template shape<0>();
    }

    auto fast(long i) const -> decltype(this->arg(fast_contiguous_slice(pythonic::builtins::None,
                                                                        pythonic::builtins::None),
                                                  i));

    auto fast(long i) -> decltype(this->arg(fast_contiguous_slice(pythonic::builtins::None,
                                                                  pythonic::builtins::None),
                                            i));
    auto fast(array_tuple<long, value> const &indices)
        -> decltype(arg.fast(array_tuple<long, 2>{{indices[1], indices[0]}}))
    {
      return arg.fast(array_tuple<long, 2>{{indices[1], indices[0]}});
    }
    auto fast(array_tuple<long, value> const &indices) const
        -> decltype(arg.fast(array_tuple<long, 2>{{indices[1], indices[0]}}))
    {
      return arg.fast(array_tuple<long, 2>{{indices[1], indices[0]}});
    }

    auto load(long i, long j) const -> decltype(arg.load(j, i))
    {
      return arg.load(j, i);
    }
    template <class Elt>
    void store(Elt elt, long i, long j)
    {
      arg.store(elt, j, i);
    }
    template <class Op, class Elt>
    void update(Elt elt, long i, long j) const
    {
      arg.template update<Op>(elt, j, i);
    }

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<numpy_texpr_2>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    /* element filtering */
    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value == 1 && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_texpr_2, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;
    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value != 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray<dtype, pshape<long>>, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !std::is_same<bool, typename F::dtype>::value &&
                         !is_pod_array<F>::value,
                     numpy_vexpr<numpy_texpr_2, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value == 1 && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_texpr_2, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;
    template <class F> // indexing through an array of boolean -- a mask
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         F::value != 1 && !is_pod_array<F>::value,
                     numpy_vexpr<ndarray<dtype, pshape<long>>, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !std::is_same<bool, typename F::dtype>::value &&
                         !is_pod_array<F>::value,
                     numpy_vexpr<numpy_texpr_2, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;
    auto operator[](long i) const -> decltype(this->fast(i));
    auto operator[](long i) -> decltype(this->fast(i));
    template <class T>
    auto operator[](array_tuple<T, value> const &indices)
        -> decltype(arg[array_tuple<T, 2>{{indices[1], indices[0]}}])
    {
      return arg[array_tuple<T, 2>{{indices[1], indices[0]}}];
    }
    template <class T>
    auto operator[](array_tuple<T, value> const &indices) const
        -> decltype(arg[array_tuple<T, 2>{{indices[1], indices[0]}}])
    {
      return arg[array_tuple<T, 2>{{indices[1], indices[0]}}];
    }
    template <class T0, class T1>
    auto operator[](std::tuple<T0, T1> const &indices)
        -> decltype(arg[std::tuple<T1, T0>{std::get<1>(indices), std::get<0>(indices)}])
    {
      return arg[std::tuple<T1, T0>{std::get<1>(indices), std::get<0>(indices)}];
    }

    template <class T0, class T1>
    auto operator[](std::tuple<T0, T1> const &indices) const
        -> decltype(arg[std::tuple<T1, T0>{std::get<1>(indices), std::get<0>(indices)}])
    {

      return arg[std::tuple<T1, T0>{std::get<1>(indices), std::get<0>(indices)}];
    }

    template <class S>
    auto operator[](S const &s0) const -> numpy_texpr<decltype(this->arg(
        fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), (s0.step, s0)))>;
    template <class S>
    auto operator[](S const &s0) -> numpy_texpr<decltype(this->arg(
        fast_contiguous_slice(pythonic::builtins::None, pythonic::builtins::None), (s0.step, s0)))>;

    template <class S, size_t... I>
    auto _reverse_index(S const &indices, std::index_sequence<I...>) const
        -> decltype(numpy::functor::transpose{}(this->arg(std::get<I>(indices)...)))
    {
      return numpy::functor::transpose{}(arg(std::get<I>(indices)...));
    }
    ndarray<dtype, typename E::shape_t> copy() const
    {
      return *this;
    }

    template <class Tp, size_t... Is>
    auto recast() -> decltype(numpy::functor::transpose{}(arg.template recast<Tp>()))
    {
      return numpy::functor::transpose{}(arg.template recast<Tp>());
    }

    template <class S0, class... S>
    auto operator()(S0 const &s0, S const &...s) const -> std::enable_if_t<
        !is_numexpr_arg<S0>::value,
        decltype(this->_reverse_index(std::tuple<S0 const &, S const &...>{s0, s...},
                                      utils::make_reversed_index_sequence<1 + sizeof...(S)>()))>;

    template <class S0, class... S>
    auto operator()(S0 const &s0, S const &...s) const
        -> std::enable_if_t<is_numexpr_arg<S0>::value, decltype(this->copy()(s0, s...))>;

    explicit operator bool() const;
    long flat_size() const;
    intptr_t id() const;
    intptr_t baseid() const
    {
      return arg.baseid();
    }
    template <class Expr>
    numpy_texpr_2 &operator=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator=(numpy_texpr<Expr> const &expr);

    template <class Op, class Expr>
    numpy_texpr_2 &update_(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator+=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator-=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator*=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator/=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator&=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator|=(Expr const &expr);
    template <class Expr>
    numpy_texpr_2 &operator^=(Expr const &expr);

    template <class NewShape>
    ndarray<dtype, NewShape> reshape(NewShape const &shape) const
    {
      return copy().reshape(shape);
    }
  };

  // only implemented for N = 2
  template <class T, class S0, class S1>
  struct numpy_texpr<ndarray<T, pshape<S0, S1>>> : numpy_texpr_2<ndarray<T, pshape<S0, S1>>> {
    numpy_texpr() = default;
    numpy_texpr(numpy_texpr const &) = default;
    numpy_texpr(numpy_texpr &&) = default;
    numpy_texpr(ndarray<T, pshape<S0, S1>> const &arg);

    numpy_texpr &operator=(numpy_texpr const &) = default;

    using numpy_texpr_2<ndarray<T, pshape<S0, S1>>>::operator=;
  };
  template <class T>
  struct numpy_texpr<ndarray<T, array_tuple<long, 2>>>
      : numpy_texpr_2<ndarray<T, array_tuple<long, 2>>> {
    numpy_texpr() = default;
    numpy_texpr(numpy_texpr const &) = default;
    numpy_texpr(numpy_texpr &&) = default;
    numpy_texpr(ndarray<T, array_tuple<long, 2>> const &arg);

    numpy_texpr &operator=(numpy_texpr const &) = default;

    using numpy_texpr_2<ndarray<T, array_tuple<long, 2>>>::operator=;
  };

  template <class E, class... S>
  struct numpy_texpr<numpy_gexpr<E, S...>> : numpy_texpr_2<numpy_gexpr<E, S...>> {
    numpy_texpr() = default;
    numpy_texpr(numpy_texpr const &) = default;
    numpy_texpr(numpy_texpr &&) = default;
    numpy_texpr(numpy_gexpr<E, S...> const &arg);
    template <class F>
    numpy_texpr(numpy_texpr<F> const &other) : numpy_texpr(numpy_gexpr<E, S...>(other.arg))
    {
    }

    numpy_texpr &operator=(numpy_texpr const &) = default;

    using numpy_texpr_2<numpy_gexpr<E, S...>>::operator=;
  };
  template <class E>
  struct numpy_texpr<broadcasted<E>> {
    static constexpr auto value = broadcasted<E>::value;
    using value_type = broadcast<typename E::dtype, typename E::dtype>;
    using dtype = typename broadcasted<E>::dtype;
    using shape_t = types::array_tuple<long, value>;
    using iterator = nditerator<numpy_texpr<broadcasted<E>>>;
    using const_iterator = const_nditerator<numpy_texpr<broadcasted<E>>>;
    // FIXME: I've got the feeling that this could be improved
    static constexpr bool is_vectorizable = false;
    static constexpr bool is_strided = true;

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<numpy_texpr<broadcasted<E>>>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const
    {
      return {*this}; // not vectorized anyway
    }
    template <class vectorizer>
    simd_iterator vend(vectorizer) const
    {
      return {*this}; // not vectorized anyway
    }
#endif

    broadcasted<E> arg;
    numpy_texpr() = default;
    numpy_texpr(numpy_texpr const &) = default;
    numpy_texpr(numpy_texpr &&) = default;
    numpy_texpr(broadcasted<E> const &arg) : arg(arg)
    {
    }
    value_type fast(long i) const
    {
      return arg.ref.fast(i);
    }
    template <size_t I>
    long shape() const
    {
      return arg.template shape<I == 0 ? 1 : 0>();
    }
    auto load(long i, long j) const -> decltype(arg.ref.load(i))
    {
      return arg.ref.load(i);
    }
    template <class Elt>
    void store(Elt elt, long i, long j)
    {
      arg.ref.store(elt, i);
    }
    const_iterator begin() const
    {
      return {*this, 0};
    }
    const_iterator end() const
    {
      return {*this, shape<0>()};
    }

    iterator begin()
    {
      return {*this, 0};
    }
    iterator end()
    {
      return {*this, shape<0>()};
    }
  };
} // namespace types

template <class Arg>
struct assignable_noescape<types::numpy_texpr<Arg>> {
  using type = types::numpy_texpr<Arg>;
};

template <class Arg>
struct assignable<types::numpy_texpr<Arg>> {
  using type = types::numpy_texpr<typename assignable<Arg>::type>;
};

template <class Arg>
struct returnable<types::numpy_texpr<Arg>> {
  using type = types::numpy_texpr<typename returnable<Arg>::type>;
};

template <class Arg>
struct lazy<types::numpy_texpr<Arg>> {
  using type = types::numpy_texpr<typename lazy<Arg>::type>;
};
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"
template <class E>
struct __combined<pythonic::types::numpy_texpr<E>, pythonic::types::numpy_texpr<E>> {
  using type = pythonic::types::numpy_texpr<E>;
};
template <class E0, class E1>
struct __combined<pythonic::types::numpy_texpr<E0>, pythonic::types::numpy_texpr<E1>> {
  using type = pythonic::types::numpy_texpr<typename __combined<E0, E1>::type>;
};

template <class E, class K>
struct __combined<pythonic::types::numpy_texpr<E>, K> {
  using type = pythonic::types::numpy_texpr<E>;
};

template <class E0, class E1, class... S>
struct __combined<pythonic::types::numpy_texpr<E0>, pythonic::types::numpy_gexpr<E1, S...>> {
  using type = pythonic::types::numpy_texpr<E0>;
};

template <class E, class O>
struct __combined<pythonic::types::numpy_texpr<E>, pythonic::types::none<O>> {
  using type = pythonic::types::none<typename __combined<pythonic::types::numpy_texpr<E>, O>::type>;
};

template <class E, class O>
struct __combined<pythonic::types::none<O>, pythonic::types::numpy_texpr<E>> {
  using type = pythonic::types::none<typename __combined<O, pythonic::types::numpy_texpr<E>>::type>;
};

template <class E>
struct __combined<pythonic::types::numpy_texpr<E>, pythonic::types::none_type> {
  using type = pythonic::types::none<pythonic::types::numpy_texpr<E>>;
};

template <class E>
struct __combined<pythonic::types::none_type, pythonic::types::numpy_texpr<E>> {
  using type = pythonic::types::none<pythonic::types::numpy_texpr<E>>;
};

/*}*/
#endif
