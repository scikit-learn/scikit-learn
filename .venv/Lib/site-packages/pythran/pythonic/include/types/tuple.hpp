#ifndef PYTHONIC_INCLUDE_TYPES_TUPLE_HPP
#define PYTHONIC_INCLUDE_TYPES_TUPLE_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/utils/seq.hpp"

#include <algorithm>
#include <tuple>

#if !defined(HAVE_SSIZE_T) || !HAVE_SSIZE_T
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif
#endif

// Equality comparison between pair && tuple
namespace std
{
  template <class F0, class S0, class F1, class S1>
  bool operator==(pair<F0, S0> const &self, tuple<F1, S1> const &other);

  template <class F0, class S0, class F1, class S1>
  bool operator==(pair<const F0, S0> const &self, tuple<F1, S1> const &other);
} // namespace std

// Tuple concatenation with operator+
template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> const &t0,
                                           std::tuple<Types1...> const &t1);

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> &&t0,
                                           std::tuple<Types1...> const &t1);

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> const &t0,
                                           std::tuple<Types1...> &&t1);

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> &&t0, std::tuple<Types1...> &&t1);

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct iterator {
    using type = T;
  };

  template <typename T>
  struct dynamic_tuple;

  template <typename T, size_t N, class V>
  struct array_base;

  struct tuple_version {
  };
  struct list_version {
  };

  template <class T, size_t N>
  using array_tuple = array_base<T, N, tuple_version>;
  template <class T, size_t N>
  using static_list = array_base<T, N, list_version>;

  template <class T>
  struct is_pod_array {
    static constexpr bool value = false;
  };
  template <typename T, size_t N, class V>
  struct is_pod_array<types::array_base<T, N, V>> {
    static constexpr bool value = true;
  };

  template <class... Tys>
  struct pshape;
  template <class T, class pS>
  struct ndarray;

  class str;

  struct slice;
  template <long stride>
  struct cstride_slice;

  /* helper to extract the tail of a tuple, && pop the head */
  template <int Offset, class T, size_t... N>
  auto make_tuple_tail(T const &t, std::index_sequence<N...>)
      -> decltype(std::make_tuple(std::get<Offset + 1 + N>(t)...))
  {
    return std::make_tuple(std::get<Offset + 1 + N>(t)...);
  }

  template <class S, class... Stail>
  std::tuple<Stail...> tuple_tail(std::tuple<S, Stail...> const &t);

  template <class... S>
  struct count_trailing_long : std::integral_constant<size_t, 0> {
  };

  template <class... S>
  struct count_trailing_long<long, S...>
      : std::integral_constant<size_t, 1 + count_trailing_long<S...>::value> {
  };

  template <class S, class... Stail>
  auto tuple_pop(std::tuple<S, Stail...> const &t)
      -> decltype(make_tuple_tail<count_trailing_long<Stail...>::value>(
          t, std::make_index_sequence<sizeof...(Stail) - count_trailing_long<Stail...>::value>{}))
  {
    return make_tuple_tail<count_trailing_long<Stail...>::value>(
        t, std::make_index_sequence<sizeof...(Stail) - count_trailing_long<Stail...>::value>{});
  }

  template <class A, size_t... I, class... Types>
  std::tuple<Types...> array_to_tuple(A const &a, std::index_sequence<I...>,
                                      utils::type_sequence<Types...>)
  {
    return std::tuple<Types...>(a[I]...);
  }

  template <class... Tys>
  struct pshape;

  template <class... Tys>
  struct iterator<pshape<Tys...>> {
    using type = array_tuple<long, sizeof...(Tys)>;
  };

  template <long N>
  long check_type(long, std::integral_constant<long, N>)
  {
    return N;
  }
  inline long check_type(long, long value)
  {
    return value;
  }

  template <long N, long P>
  std::integral_constant<long, N> check_type(std::integral_constant<long, N>,
                                             std::integral_constant<long, P>)
  {
    assert(N == P && "consistent init");
    return {};
  }
  template <long N>
  std::integral_constant<long, N> check_type(std::integral_constant<long, N>, long v)
  {
    assert(N == v && "consistent init");
    return {};
  }

  template <class T>
  struct is_pshape_element : std::is_integral<T> {
  };
  template <long N>
  struct is_pshape_element<std::integral_constant<long, N>> : std::true_type {
  };

  template <class... Tys>
  struct pshape {
    static_assert(utils::all_of<is_pshape_element<Tys>::value...>::value, "valid pshape");

    std::tuple<Tys...> values;

    template <class... Args, size_t... Is>
    pshape(std::tuple<Args...> const &v, std::index_sequence<Is...>)
        : values{check_type(std::get<Is>(values), std::get<Is>(v))...}
    {
    }
    template <class... Args>
    pshape(std::tuple<Args...> const &v) : pshape(v, std::make_index_sequence<sizeof...(Args)>())
    {
    }

    template <class... Args>
    pshape(long arg, Args... args)
        : pshape(std::make_tuple(arg, args...), std::make_index_sequence<1 + sizeof...(args)>())
    {
    }
    template <class T, T N, class... Args>
    pshape(std::integral_constant<T, N> arg, Args... args)
        : pshape(std::make_tuple(arg, args...), std::make_index_sequence<1 + sizeof...(args)>())
    {
    }

    template <class S, size_t... Is>
    pshape(S const *buffer, std::index_sequence<Is...>)
        : values{check_type(std::get<Is>(values), buffer[Is])...}
    {
    }
    template <class S>
    pshape(S const *buffer) : pshape(buffer, std::make_index_sequence<sizeof...(Tys)>())
    {
    }
    template <class... TyOs>
    pshape(pshape<TyOs...> other)
        : pshape(other.values, std::make_index_sequence<sizeof...(TyOs)>())
    {
      static_assert(sizeof...(TyOs) == sizeof...(Tys), "compatible sizes");
    }

    template <class S, class V>
    pshape(pythonic::types::array_base<S, sizeof...(Tys), V> data) : pshape(data.data())
    {
    }

    pshape() = default;
    pshape(pshape const &) = default;
    pshape(pshape &&) = default;
    pshape &operator=(pshape const &) = default;
    pshape &operator=(pshape &&) = default;

    template <size_t... Is>
    types::array_tuple<long, sizeof...(Tys)> array(std::index_sequence<Is...>) const
    {
      return {{get<Is>()...}};
    }

    types::array_tuple<long, sizeof...(Tys)> array() const
    {
      return array(std::make_index_sequence<sizeof...(Tys)>());
    }
    operator types::array_tuple<long, sizeof...(Tys)>() const
    {
      return array();
    }

    template <size_t I>
    long get() const
    {
      return std::get<I>(values);
    }
    template <size_t I>
    auto get() -> decltype(std::get<I>(values))
    {
      return std::get<I>(values);
    }
  };
  template <class P, size_t M, class... Ss>
  struct shape_builder;

  template <class P, size_t I, class V, size_t M, class... Ss>
  struct shape_builder<array_base<P, I, V>, M, Ss...>
      : shape_builder<P, M - 1, Ss..., std::integral_constant<long, I>> {
  };

  template <class P, class... Ss>
  struct shape_builder<P, 0, Ss...> {
    using type = pshape<Ss...>;
  };

  template <class P, size_t M, class... Ss>
  struct shape_builder : shape_builder<typename P::value_type, M - 1, Ss..., long> {
  };

  struct array_base_slicer {
    template <class T, size_t N>
    dynamic_tuple<T> operator()(array_tuple<T, N> const &b, slice const &s);

    template <class T, size_t N, long stride>
    dynamic_tuple<T> operator()(array_tuple<T, N> const &b, cstride_slice<stride> const &s);

    template <class T, size_t N>
    dynamic_tuple<T> operator()(array_tuple<T, N> const &b, fast_contiguous_slice const &s);

    template <class T, size_t N, class S>
    std::enable_if_t<is_slice<S>::value, sliced_list<T, S>> operator()(static_list<T, N> const &b,
                                                                       S const &s)
    {
      return {b, s};
    }
  };

  namespace details
  {
    template <class E>
    auto extract_shape(E const &e, utils::int_<0>) -> decltype(e.size())
    {
      return e.size();
    }
    template <class E, size_t L>
    auto extract_shape(E const &e, utils::int_<L>)
        -> decltype(extract_shape(e[0], utils::int_<L - 1>{}))
    {
      return extract_shape(e[0], utils::int_<L - 1>{});
    }
  } // namespace details

  /* inspired by std::array implementation */
  template <typename T, size_t N, typename Version>
  struct array_base {
    using value_type = T;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator = value_type *;
    using const_iterator = const value_type *;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    // minimal ndarray interface
    using dtype = typename utils::nested_container_value_type<array_base>::type;
    static const size_t value = utils::nested_container_depth<array_base>::value;
    static const bool is_vectorizable = true;
    static const bool is_flat = true;
    static const bool is_strided = false;

    // flat_size implementation
    template <class E>
    long _flat_size(E const &e, utils::int_<1>) const;
    template <class E, size_t L>
    long _flat_size(E const &e, utils::int_<L>) const;

    long flat_size() const;

    // Support for zero-sized arrays mandatory.
    value_type buffer[N ? N : 1];

    // No explicit construct/copy/destroy for aggregate type.

    void fill(const value_type &__u);

    long count(value_type const &u) const
    {
      return std::count(begin(), end(), u);
    }

    // Iterators.
    iterator begin() noexcept;
    const_iterator begin() const noexcept;

    iterator end() noexcept;
    const_iterator end() const noexcept;

    reverse_iterator rbegin() noexcept;
    const_reverse_iterator rbegin() const noexcept;

    reverse_iterator rend() noexcept;
    const_reverse_iterator rend() const noexcept;

    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

    const_reverse_iterator crbegin() const noexcept;
    const_reverse_iterator crend() const noexcept;

    // Capacity.
    constexpr size_type size() const noexcept;
    constexpr size_type max_size() const noexcept;
    constexpr bool empty() const noexcept;

    intptr_t id() const;

    // Element access.
    reference fast(long n);
    const_reference fast(long n) const noexcept;
#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<array_base>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif
    template <class... Indices>
    dtype load(long index0, long index1, Indices... indices) const
    {
      return fast(index0).load(index1, indices...);
    }

    dtype load(long index) const
    {
      return fast(index);
    }

    reference operator[](long __n);

    const_reference operator[](long __n) const noexcept;

    template <class S>
    auto operator[](S s) const -> decltype(array_base_slicer{}(*this, (s.lower, s)))
    {
      return array_base_slicer{}(*this, s);
    }

    reference front();
    const_reference front() const;

    reference back();
    const_reference back() const;

    pointer data() noexcept;
    const_pointer data() const noexcept;

    // operator

    // for conversion to dict item type
    template <class K, class V>
    operator std::pair<const K, V>() const
    {
      static_assert(std::is_same<K, T>::value && std::is_same<V, T>::value && N == 2,
                    "compatible conversion");
      return {data()[0], data()[1]};
    }

    template <size_t M>
    bool operator==(array_base<T, M, Version> const &other) const;

    template <size_t M>
    bool operator!=(array_base<T, M, Version> const &other) const;

    template <size_t M>
    bool operator<(array_base<T, M, Version> const &other) const;

    template <class Tp, size_t M>
    array_base<typename __combined<T, Tp>::type, N + M, Version>
    operator+(array_base<Tp, M, Version> const &other) const;

    // tuple conversion
    template <class... Types>
    operator std::tuple<Types...>() const;

    template <class Tp>
    operator array_base<Tp, N, Version>() const;

    auto to_tuple() const -> decltype(array_to_tuple(*this, std::make_index_sequence<N>{},

                                                     utils::make_repeated_type<T, N>()));

    template <class W>
    array_base<T, N, W> to_array() const;

    template <class W>
    explicit operator array_base<T, N, W>() const
    {
      return to_array<W>();
    }

    template <class S>
    auto operator()(S const &s) const -> decltype((*this)[s])
    {
      return (*this)[s];
    }

    bool operator!() const
    {
      return N == 0;
    }

    /* array */
    template <class T1, size_t N1, class Version1>
    friend std::ostream &operator<<(std::ostream &os, types::array_base<T1, N1, Version1> const &v);

    using shape_t = typename shape_builder<array_base, value>::type;
    template <size_t I>
    auto shape() const -> decltype(details::extract_shape(*this, utils::int_<I>{}))
    {
      return details::extract_shape(*this, utils::int_<I>{});
    }
  };

  // Implementation for detection of "same type".
  // With this information, we know if we must create a real tuple || a
  // static sized array
  namespace details
  {

    template <class... Types>
    struct alike;

    template <>
    struct alike<> {
      static bool const value = false;
      using type = void;
    };

    template <class T>
    struct alike<T> {
      static bool const value = true;
      using type = std::remove_cv_t<std::remove_reference_t<T>>;
    };
    template <class A, class... S>
    struct alike<numpy_gexpr<A, S...>, numpy_gexpr<A const &, S...>> {
      static bool const value = true;
      using type = numpy_gexpr<A, S...>;
    };

    template <class T0, class T1>
    struct alike<T0, T1> {
      static bool const value = std::is_same<T0, T1>::value;
      using type = std::conditional_t<value, T0, void>;
    };

    // specialization to make static string alike types::str
    template <size_t N>
    struct alike<char[N], str> {
      static bool const value = true;
      using type = str;
    };

    template <size_t N>
    struct alike<str, char[N]> {
      static bool const value = true;
      using type = str;
    };

    template <size_t N, size_t M>
    struct alike<char[M], char[N]> {
      static bool const value = true;
      using type = str;
    };

    template <class T, size_t N, class V, class... Types>
    struct alike<std::tuple<Types...>, array_base<T, N, V>> {
      static bool const value =
          sizeof...(Types) == N &&
          alike<T, std::remove_cv_t<std::remove_reference_t<Types>>...>::value;
      using type = std::conditional_t<
          value, typename alike<T, std::remove_cv_t<std::remove_reference_t<Types>>...>::type,
          void>;
    };

    template <class T, size_t N, class V, class... Types>
    struct alike<array_base<T, N, V>, std::tuple<Types...>>
        : alike<std::tuple<Types...>, array_base<T, N, V>> {
    };

    template <class T, class... Types>
    struct alike<T, Types...> {
      static bool const value =
          alike<Types...>::value && alike<T, typename alike<Types...>::type>::value;
      using type = typename alike<T, typename alike<Types...>::type>::type;
    };
  } // namespace details

  template <class... Types>
  struct alike : details::alike<std::remove_cv_t<std::remove_reference_t<Types>>...> {
  };

  // Pythonic implementation for make_tuple to have the best return type
  // (static array for sames types || real tuple otherwise)
  template <bool Same, class... Types>
  struct _make_tuple {
    auto operator()(Types &&...types) -> decltype(std::make_tuple(std::forward<Types>(types)...))
    {
      return std::make_tuple(std::forward<Types>(types)...);
    }
  };

  template <class... Types>
  struct _make_tuple<true, Types...> {
    types::array_tuple<typename alike<Types...>::type, sizeof...(Types)>
    operator()(Types &&...types)
    {
      return {{std::forward<Types>(types)...}};
    }
  };

  template <class... Types>
  auto make_tuple(Types &&...types)
#if !_MSC_VER || __clang__
      -> decltype(_make_tuple<alike<Types...>::value, Types...>()(std::forward<Types>(types)...))
#endif
  {
    return _make_tuple<alike<Types...>::value, Types...>()(std::forward<Types>(types)...);
  }

  template <class... Tys>
  using make_tuple_t = decltype(types::make_tuple(std::declval<Tys>()...));

  template <class T, class Tuple, size_t... S>
  types::array_tuple<T, sizeof...(S)> _to_array(Tuple const &t, std::index_sequence<S...>)
  {
    return {{static_cast<T>(std::get<S>(t))...}};
  }

  template <class T, class... Tys>
  types::array_tuple<T, sizeof...(Tys)> to_array(std::tuple<Tys...> const &t)
  {
    return _to_array<T>(t, std::make_index_sequence<sizeof...(Tys)>());
  }

  // Tuple concatenation for array && tuple
  template <class T, size_t N, class V, class... Types>
  auto operator+(std::tuple<Types...> const &t, types::array_base<T, N, V> const &lt)
      -> decltype(std::tuple_cat(t, lt.to_tuple()));

  template <class T, size_t N, class V, class... Types>
  auto operator+(types::array_base<T, N, V> const &lt, std::tuple<Types...> const &t)
      -> decltype(std::tuple_cat(lt.to_tuple(), t));
} // namespace types

template <class... Types>
struct assignable<std::tuple<Types...>> {
  using type = std::tuple<typename assignable<Types>::type...>;
};

template <typename T, size_t N, class V>
struct assignable<pythonic::types::array_base<T, N, V>> {
  using type = pythonic::types::array_base<typename assignable<T>::type, N, V>;
};

template <class... Types>
struct returnable<std::tuple<Types...>> {
  using type = std::tuple<typename returnable<Types>::type...>;
};

template <typename T, size_t N, class V>
struct returnable<pythonic::types::array_base<T, N, V>> {
  using type = pythonic::types::array_base<typename returnable<T>::type, N, V>;
};
PYTHONIC_NS_END

/* specialize std::get */
namespace std
{

  template <size_t I, class T, size_t N, class V>
  typename pythonic::types::array_base<T, N, V>::reference
  get(pythonic::types::array_base<T, N, V> &t)
  {
    return t[I];
  }

  template <size_t I, class T, size_t N, class V>
  typename pythonic::types::array_base<T, N, V>::const_reference
  get(pythonic::types::array_base<T, N, V> const &t)
  {
    return t[I];
  }

  template <size_t I, class T, size_t N, class V>
  typename pythonic::types::array_base<T, N, V>::value_type
  get(pythonic::types::array_base<T, N, V> &&t)
  {
    return std::move(t)[I];
  }

  template <size_t I, class T, size_t N, class V>
  struct tuple_element<I, pythonic::types::array_base<T, N, V>> {
    using type = typename pythonic::types::array_base<T, N, V>::value_type;
  };

  template <typename T, size_t N, class V>
  struct tuple_size<pythonic::types::array_base<T, N, V>> {
    static const size_t value = N;
  };
} // namespace std

/* hashable tuples, as proposed in
 * http://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
 */
namespace
{
  inline size_t hash_combiner(size_t left, size_t right); // replacable

  template <size_t index, class... types>
  struct hash_impl {
    size_t operator()(size_t a, const std::tuple<types...> &t) const;
  };

  template <class... types>
  struct hash_impl<0, types...> {
    size_t operator()(size_t a, const std::tuple<types...> &t) const;
  };
} // namespace

/* specialize std::hash */
namespace std
{
  template <class... Types>
  struct hash<std::tuple<Types...>> {
    size_t operator()(std::tuple<Types...> const &t) const;
  };

  template <typename T, size_t N, class V>
  struct hash<pythonic::types::array_base<T, N, V>> {
    size_t operator()(pythonic::types::array_base<T, N, V> const &l) const;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"
template <class K, class... Types>
struct __combined<indexable<K>, std::tuple<Types...>> {
  using type = std::tuple<Types...>;
};

template <class K, class... Types>
struct __combined<std::tuple<Types...>, indexable<K>> {
  using type = std::tuple<Types...>;
};

template <class T, size_t N>
struct __combined<pythonic::types::static_list<T, N>, pythonic::types::static_list<T, N>> {
  using type = pythonic::types::static_list<T, N>;
};

template <class T, size_t N>
struct __combined<pythonic::types::array_tuple<T, N>, pythonic::types::array_tuple<T, N>> {
  using type = pythonic::types::array_tuple<T, N>;
};

template <class T0, class T1, size_t N, class V>
struct __combined<pythonic::types::array_base<T0, N, V>, pythonic::types::array_base<T1, N, V>> {
  using type = pythonic::types::array_base<typename __combined<T0, T1>::type, N, V>;
};

template <class T0, class T1, size_t N>
struct __combined<pythonic::types::static_list<T0, N>, pythonic::types::static_list<T1, N>> {
  using type = pythonic::types::static_list<typename __combined<T0, T1>::type, N>;
};

template <class T0, class T1, size_t N0, size_t N1>
struct __combined<pythonic::types::static_list<T0, N0>, pythonic::types::static_list<T1, N1>> {
  using type = pythonic::types::list<typename __combined<T0, T1>::type>;
};

template <class K, class T, size_t N, class V>
struct __combined<indexable<K>, pythonic::types::array_base<T, N, V>> {
  using type = pythonic::types::array_base<T, N, V>;
};

template <class K, class T, size_t N, class V>
struct __combined<pythonic::types::array_base<T, N, V>, indexable<K>> {
  using type = pythonic::types::array_base<T, N, V>;
};

template <class K, class T, size_t N, class V>
struct __combined<container<K>, pythonic::types::array_base<T, N, V>> {
  using type = pythonic::types::array_base<typename __combined<T, K>::type, N, V>;
};

template <class K, class T, size_t N, class V>
struct __combined<pythonic::types::array_base<T, N, V>, container<K>> {
  using type = pythonic::types::array_base<typename __combined<T, K>::type, N, V>;
};

template <class K, class V, class T, size_t N, class AV>
struct __combined<indexable_container<K, V>, pythonic::types::array_base<T, N, AV>> {
  using type = pythonic::types::array_base<typename __combined<V, T>::type, N, AV>;
};

template <class K, class V, class T, size_t N, class AV>
struct __combined<pythonic::types::array_base<T, N, AV>, indexable_container<K, V>> {
  using type = pythonic::types::array_base<typename __combined<T, V>::type, N, AV>;
};

template <class... t0, class... t1>
struct __combined<std::tuple<t0...>, std::tuple<t1...>> {
  using type = std::tuple<typename __combined<t0, t1>::type...>;
};

template <class t, class... t0>
struct __combined<std::tuple<t0...>, container<t>> {
  using type = std::tuple<t0...>;
};

template <class t, class... t0>
struct __combined<container<t>, std::tuple<t0...>> {
  using type = std::tuple<t0...>;
};

PYTHONIC_NS_BEGIN
namespace details
{
  template <class T, class P, bool Combine>
  struct pick_combined;
  template <class T, class P>
  struct pick_combined<T, P, true> {
    using type = typename __combined<T, P>::type;
  };
  template <class T, class P>
  struct pick_combined<T, P, false> {
    using type = T;
  };
} // namespace details
PYTHONIC_NS_END

template <long I, class t, class... t0>
struct __combined<std::tuple<t0...>, indexable_container<std::integral_constant<long, I>, t>> {
  using holder = std::tuple<t0...>;
  template <size_t... Is>
  static std::tuple<typename pythonic::details::pick_combined<std::tuple_element_t<Is, holder>, t,
                                                              I == Is>::type...>
      make_type(std::index_sequence<Is...>);
  static auto make_type() -> decltype(make_type(std::make_index_sequence<sizeof...(t0)>()));
  using type = decltype(make_type());
};

template <class k, class t, class... t0>
struct __combined<indexable_container<k, t>, std::tuple<t0...>>
    : __combined<std::tuple<t0...>, indexable_container<k, t>> {
};

template <class t, size_t n, class... types>
struct __combined<pythonic::types::array_tuple<t, n>, std::tuple<types...>> {
  using type = std::tuple<typename __combined<t, types>::type...>;
};

template <class t, size_t n, class... types>
struct __combined<pythonic::types::array_tuple<t, n>, pythonic::types::pshape<types...>> {
  using type = pythonic::types::array_tuple<t, n>;
};
template <class t, size_t n, class... types>
struct __combined<pythonic::types::pshape<types...>, pythonic::types::array_tuple<t, n>> {
  using type = pythonic::types::array_tuple<t, n>;
};

template <class t, size_t n, class... types>
struct __combined<std::tuple<types...>, pythonic::types::array_tuple<t, n>> {
  using type = std::tuple<typename __combined<types, t>::type...>;
};
template <class t00, class t01, class t10, class t11>
struct __combined<std::pair<t00, t01>, std::pair<t10, t11>> {
  using type = std::pair<typename __combined<t00, t10>::type,
                         typename __combined<t01, t11>::type>; // no further combination
};

/* } */

PYTHONIC_NS_BEGIN

namespace types
{

  template <class Tuple, size_t I>
  void print_tuple(std::ostream &os, Tuple const &t, utils::int_<I>);

  template <class Tuple>
  void print_tuple(std::ostream &os, Tuple const &t, utils::int_<0>);

  template <typename T, size_t N, class V>
  struct len_of<array_base<T, N, V>> {
    static constexpr long value = N;
  };
  template <typename T, long I, class... Is>
  struct len_of<ndarray<T, pshape<std::integral_constant<long, I>, Is...>>> {
    static constexpr long value = I;
  };

  template <class... Types>
  struct len_of<std::tuple<Types...>> {
    static constexpr long value = sizeof...(Types);
  };
} // namespace types
PYTHONIC_NS_END

namespace std
{
  template <class... Args>
  ostream &operator<<(ostream &os, tuple<Args...> const &t);
  template <size_t I, class... Tys>
  long get(pythonic::types::pshape<Tys...> const &s)
  {
    return s.template get<I>();
  }
  template <size_t I, class... Tys>
  auto get(pythonic::types::pshape<Tys...> &s) -> decltype(s.template get<I>())
  {
    return s.template get<I>();
  }
  template <size_t I, class T>
  auto get(T *s) -> decltype(s[I])
  {
    return s[I];
  }
  template <size_t I, class T>
  long get(T const *s)
  {
    return s[I];
  }

  template <class... Tys>
  struct tuple_size<pythonic::types::pshape<Tys...>>
      : public std::integral_constant<std::size_t, sizeof...(Tys)> {
  };

  template <size_t I, class... Tys>
  struct tuple_element<I, pythonic::types::pshape<Tys...>> {
    using type = std::tuple_element_t < I<sizeof...(Tys) ? I : 0, std::tuple<Tys...>>;
  };
} // namespace std
PYTHONIC_NS_BEGIN
namespace sutils
{

  template <class T>
  struct make_shape {
    using type = T;
  };

  template <typename T, size_t N, class V>
  struct make_shape<types::array_base<T, N, V>> {
    using type = types::array_tuple<long, N>;
  };

  template <class T>
  using shape_t = std::enable_if_t<!std::is_integral<T>::value, typename make_shape<T>::type>;

  template <class Curr, class... Ss>
  struct shape_merger;
  template <class Curr>
  struct shape_merger<Curr> {
    using type = Curr;
  };

  template <class Curr, class... Ss>
  struct shape_merger<Curr, long, Ss...> {
    using type = long;
  };
  template <long N0, long N1, class... Ss>
  struct shape_merger<std::integral_constant<long, N0>, std::integral_constant<long, N1>, Ss...>
      : shape_merger<std::integral_constant<long, (N0 > N1 ? N0 : N1)>, Ss...> {
  };
  template <long N, class... Ss>
  struct shape_merger<long, std::integral_constant<long, N>, Ss...> {
    using type = long;
  };

  template <size_t I, class Ss>
  struct shape_selecter
      : std::conditional<(I < std::tuple_size<Ss>::value),
                         std::tuple_element_t<(I < std::tuple_size<Ss>::value ? I : 0L), Ss>,
                         std::integral_constant<long, 1>> {
  };

  template <size_t I, class Ss>
  struct merge_shape;
  template <size_t I, class... Ss>
  struct merge_shape<I, std::tuple<Ss...>> {
    using type = typename shape_merger<typename shape_selecter<I, Ss>::type...>::type;
  };
  template <class Ss, class T>
  struct merged_shapes;

  template <class Ss, size_t... Is>
  struct merged_shapes<Ss, std::index_sequence<Is...>> {
    using type = types::pshape<typename merge_shape<Is, Ss>::type...>;
  };

  template <size_t N, class... Ss>
  using merged_shapes_t =
      typename merged_shapes<std::tuple<Ss...>, std::make_index_sequence<N>>::type;

  template <class... Ss>
  struct shape_commonifier;
  template <class Ss>
  struct shape_commonifier<Ss> {
    using type = Ss;
  };
  template <class S1, class... Ss>
  struct shape_commonifier<long, S1, Ss...> {
    using type = long;
  };
  template <long N, class... Ss>
  struct shape_commonifier<std::integral_constant<long, N>, long, Ss...> {
    using type = long;
  };
  template <long N0, long N1, class... Ss>
  struct shape_commonifier<std::integral_constant<long, N0>, std::integral_constant<long, N1>,
                           Ss...> {
    using type = std::conditional_t<
        N0 == N1, typename shape_commonifier<std::integral_constant<long, N0>, Ss...>::type, long>;
  };

  template <size_t I, class Ss>
  struct common_shape;
  template <size_t I, class... Ss>
  struct common_shape<I, std::tuple<Ss...>> {
    using type = typename shape_commonifier<std::tuple_element_t<I, Ss>...>::type;
  };

  template <class Ss, class T>
  struct common_shapes;

  template <class Ss, size_t... Is>
  struct common_shapes<Ss, std::index_sequence<Is...>> {
    using type = types::pshape<typename common_shape<Is, Ss>::type...>;
  };
  template <size_t N, class... Ss>
  using common_shapes_t =
      typename common_shapes<std::tuple<Ss...>, std::make_index_sequence<N>>::type;

  template <class T>
  struct transpose;
  template <class T>
  struct transpose<types::array_tuple<T, 2>> {
    using type = types::array_tuple<T, 2>;
  };

  template <class T0, class T1>
  struct transpose<types::pshape<T0, T1>> {
    using type = types::pshape<T1, T0>;
  };
  template <class T>
  using transpose_t = typename transpose<T>::type;

  template <class T0, class T1>
  void assign(T0 &t0, T1 t1)
  {
    t0 = (T0)t1;
  }
  template <class T0, T0 N, class T1>
  void assign(std::integral_constant<T0, N> &t0, T1 t1)
  {
    assert((long)t0 == (long)t1 && "consistent");
  }

  template <size_t Start, ssize_t Offset, class T0, class T1, size_t... Is>
  void copy_shape(T0 &shape0, T1 const &shape1, std::index_sequence<Is...>)
  {
    (void)std::initializer_list<int>{
        (assign(std::get<Start + Is>(shape0), shape1.template shape<Is + Start + Offset>()), 1)...};
  }
  template <size_t Start, ssize_t Offset, class T0, class T1, size_t... Is>
  void scopy_shape(T0 &shape0, T1 const &shape1, std::index_sequence<Is...>)
  {
    (void)std::initializer_list<int>{
        (assign(std::get<Start + Is>(shape0), std::get<Is + Start + Offset>(shape1)), 1)...};
  }
  template <size_t Start, ssize_t Offset, class T0, class T1, size_t... Is>
  void copy_strides(T0 &stride0, T1 const &stride1, std::index_sequence<Is...>)
  {
    (void)std::initializer_list<int>{
        (assign(std::get<Start + Is>(stride0), stride1.template strides<Is + Start + Offset>()),
         1)...};
  }
  template <class P, class... Tys>
  struct pop_type;

  template <class... Ps, class Ty>
  struct pop_type<types::pshape<Ps...>, Ty> {
    using type = types::pshape<Ps...>;
  };
  template <class... Ps, class Ty, class... Tys>
  struct pop_type<types::pshape<Ps...>, Ty, Tys...> : pop_type<types::pshape<Ps..., Ty>, Tys...> {
  };

  template <class T>
  struct pop_tail;

  template <class... Tys>
  struct pop_tail<types::pshape<Tys...>> {
    using type = typename pop_type<types::pshape<>, Tys...>::type;
  };
  template <typename T, size_t N, class V>
  struct pop_tail<types::array_base<T, N, V>> {
    using type = types::array_tuple<T, N - 1>;
  };

  template <class T>
  struct pop_head;

  template <class Ty, class... Tys>
  struct pop_head<types::pshape<Ty, Tys...>> {
    using type = types::pshape<Tys...>;
  };
  template <typename T, size_t N, class V>
  struct pop_head<types::array_base<T, N, V>> {
    using type = types::array_tuple<T, N - 1>;
  };

  template <class T>
  struct head;

  template <class Ty, class... Tys>
  struct head<types::pshape<Ty, Tys...>> {
    using type = Ty;
  };
  template <typename T, size_t N, class V>
  struct head<types::array_base<T, N, V>> {
    using type = T;
  };

  template <class T>
  using pop_head_t = typename pop_head<T>::type;

  template <class T>
  using pop_tail_t = typename pop_tail<T>::type;

  template <class T>
  using head_t = typename head<T>::type;

  template <class... Tys>
  types::array_tuple<long, sizeof...(Tys)> array(types::pshape<Tys...> const &pS)
  {
    return pS.array();
  }

  template <typename T, size_t N, class V>
  types::array_base<T, N, V> array(types::array_base<T, N, V> const &pS)
  {
    return pS;
  }
  template <class E, size_t... Is>
  types::array_tuple<long, sizeof...(Is)> getshape(E const &e, std::index_sequence<Is...>)
  {
    return {(long)(e.template shape<Is>())...};
  }
  template <class E>
  auto getshape(E const &e) -> decltype(getshape(e, std::make_index_sequence<E::value>()))
  {
    return getshape(e, std::make_index_sequence<E::value>());
  }

  inline std::tuple<> getshape(...)
  {
    return {};
  }

  template <class pS0, class pS1>
  struct concat;

  template <class... Ty0s, class... Ty1s>
  struct concat<types::pshape<Ty0s...>, types::pshape<Ty1s...>> {
    using type = types::pshape<Ty0s..., Ty1s...>;
  };

  template <class... Tys>
  struct concat<types::pshape<Tys...>, types::array_tuple<long, 0>> {
    using type = types::pshape<Tys...>;
  };
  template <class... Tys, size_t N>
  struct concat<types::pshape<Tys...>, types::array_tuple<long, N>>
      : concat<types::pshape<Tys..., long>, types::array_tuple<long, N - 1>> {
  };

  template <class... Ty1s>
  struct concat<types::array_tuple<long, 0>, types::pshape<Ty1s...>> {
    using type = types::pshape<Ty1s...>;
  };

  template <size_t N, class... Ty1s>
  struct concat<types::array_tuple<long, N>, types::pshape<Ty1s...>>
      : concat<types::array_tuple<long, N - 1>, types::pshape<long, Ty1s...>> {
  };

  template <size_t N, size_t M>
  struct concat<types::array_tuple<long, N>, types::array_tuple<long, M>> {
    using type = types::array_tuple<long, N + M>;
  };

  template <class... Tys>
  using concat_t = typename concat<Tys...>::type;

  template <class P, class T>
  using push_front_t = concat_t<types::pshape<T>, P>;

  template <class S>
  long find(S &s, long v, std::integral_constant<size_t, 0>, long start, bool comp(long, long))
  {
    return comp(s.template shape<0>(), v) && 0 < start ? 0 : -1;
  }
  template <class S, size_t I>
  long find(S &s, long v, std::integral_constant<size_t, I>, long start, bool comp(long, long))
  {
    return comp(s.template shape<I>(), v) && I < start
               ? I
               : find(s, v, std::integral_constant<size_t, I - 1>(), start, comp);
  }

  template <class S>
  long find(
      S &s, long v, long start = S::value,
      bool comp(long, long) = [](long a, long b) { return (a == b); })
  {
    return find(s, v, std::integral_constant<size_t, S::value - 1>(), start, comp);
  }
  template <class S>
  long sfind(S &s, long v, std::integral_constant<size_t, 0>, long start, bool comp(long, long))
  {
    return comp(std::get<0>(s), v) && 0 < start ? 0 : -1;
  }
  template <class S, size_t I>
  long sfind(S &s, long v, std::integral_constant<size_t, I>, long start, bool comp(long, long))
  {
    return comp(std::get<I>(s), v) && (long)I < start
               ? (long)I
               : sfind(s, v, std::integral_constant<size_t, I - 1>(), start, comp);
  }

  template <class S>
  long sfind(
      S &s, long v, long start = std::tuple_size<S>::value,
      bool comp(long, long) = [](long a, long b) { return (a == b); })
  {
    return sfind(s, v, std::integral_constant<size_t, std::tuple_size<S>::value - 1>(), start,
                 comp);
  }

  template <class S, class B>
  bool equals(S const &s, B const &other, std::integral_constant<size_t, 0>)
  {
    return std::get<0>(other) == s.template shape<0>();
  }
  template <class S, class B, size_t I>
  bool equals(S const &s, B const &other, std::integral_constant<size_t, I>)
  {
    return std::get<I>(other) == s.template shape<I>() &&
           equals(s, other, std::integral_constant<size_t, I - 1>());
  }

  template <class S, class B>
  std::enable_if_t<S::value == std::tuple_size<B>::value, bool> equals(S const &s, B const &other)
  {
    return equals(s, other, std::integral_constant<size_t, S::value - 1>());
  }
  template <class S, class B>
  std::enable_if_t<std::tuple_size<S>::value != std::tuple_size<B>::value, bool>
  equals(S const &s, B const &other)
  {
    return false;
  }
  template <class S, class B>
  bool equals(S const &s, B *other)
  {
    return equals(s, other, std::integral_constant<size_t, S::value - 1>());
  }
  template <class S, class B>
  bool requals(S const &s, B const *other, std::integral_constant<size_t, 0>)
  {
    return other[S::value - 1] == s.template shape<0>();
  }
  template <class S, class B, size_t I>
  bool requals(S const &s, B const *other, std::integral_constant<size_t, I>)
  {
    return other[S::value - I - 1] == s.template shape<I>() &&
           requals(s, other, std::integral_constant<size_t, I - 1>());
  }
  template <class S, class B>
  bool requals(S const &s, B const *other)
  {
    return requals(s, other, std::integral_constant<size_t, S::value - 1>());
  }

  template <class S, class P>
  bool any_of(S const &s, P pred, std::integral_constant<size_t, 0>)
  {
    return pred(s.template shape<0>());
  }
  template <class S, class P, size_t I>
  bool any_of(S const &s, P pred, std::integral_constant<size_t, I>)
  {
    return pred(s.template shape<I>()) || any_of(s, pred, std::integral_constant<size_t, I - 1>());
  }
  template <class S, class Pred>
  bool any_of(S const &s, Pred pred)
  {
    return any_of(s, pred, std::integral_constant<size_t, S::value - 1>());
  }

  template <class S>
  long min(long curr, S const &s, std::integral_constant<size_t, 0>)
  {
    return std::min(curr, s.template shape<0>());
  }
  template <class S, size_t I>
  long min(long curr, S const &s, std::integral_constant<size_t, I>)
  {
    return min(std::min(curr, s.template shape<I>()), s, std::integral_constant<size_t, I - 1>());
  }
  template <class S>
  long min(S const &s)
  {
    return min(s.template shape<S::value - 1>(), s, std::integral_constant<size_t, S::value - 1>());
  }

  template <class S>
  long prod(S const &s, std::integral_constant<size_t, 0>)
  {
    return s.template shape<0>();
  }
  template <class S, size_t I>
  long prod(S const &s, std::integral_constant<size_t, I>)
  {
    return s.template shape<I>() * prod(s, std::integral_constant<size_t, I - 1>());
  }
  template <class S>
  long prod(S const &s)
  {
    return prod(s, std::integral_constant<size_t, S::value - 1>());
  }
  template <class S>
  long sprod(S const &s, std::integral_constant<size_t, 0>)
  {
    return std::get<0>(s);
  }
  template <class S, size_t I>
  long sprod(S const &s, std::integral_constant<size_t, I>)
  {
    return std::get<I>(s) * sprod(s, std::integral_constant<size_t, I - 1>());
  }
  template <class S>
  long sprod(S const &s)
  {
    return sprod(s, std::integral_constant<size_t, std::tuple_size<S>::value - 1>());
  }

  template <class S>
  long prod_tail(S, std::integral_constant<size_t, 0>)
  {
    return 1;
  }
  template <class S, size_t I>
  long prod_tail(S const &s, std::integral_constant<size_t, I>)
  {
    return s.template shape<I>() * prod_tail(s, std::integral_constant<size_t, I - 1>());
  }
  template <class S>
  long prod_tail(S const &s)
  {
    return prod_tail(s, std::integral_constant<size_t, S::value - 1>());
  }
  template <class S>
  long prod_head(S, std::integral_constant<size_t, 0>)
  {
    return 1;
  }
  template <class S, size_t I>
  long prod_head(S const &s, std::integral_constant<size_t, I>)
  {
    return s.template shape<S::value - 1 - I>() *
           prod_head(s, std::integral_constant<size_t, I - 1>());
  }
  template <class S>
  long prod_head(S const &s)
  {
    return prod_head(s, std::integral_constant<size_t, S::value - 1>());
  }

  template <size_t I, class P>
  struct safe_tuple_element {
    using type = std::tuple_element_t<(I < std::tuple_size<P>::value ? I : 0), P>;
  };

  template <size_t I>
  struct copy_new_axis_helper;

  template <>
  struct copy_new_axis_helper<0> {
    template <class S0, class S1, class S2, size_t J>
    std::enable_if_t<(0 != std::tuple_size<S2>::value) && std::tuple_element_t<0, S2>::value,
                     sutils::push_front_t<S0, std::integral_constant<long, 1>>>
    doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
    {
      return {std::tuple_cat(std::tuple<std::integral_constant<long, 1>>(), s.values)};
    }
    template <class S0, class S1, class S2, size_t J>
    std::enable_if_t<(0 != std::tuple_size<S2>::value) && !std::tuple_element_t<0, S2>::value,
                     sutils::push_front_t<S0, std::tuple_element_t<0, typename S1::shape_t>>>
    doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
    {
      return {std::tuple_cat(std::make_tuple(shape.template shape<0>()), s.values)};
    }

    template <class S0, class S1, class S2, size_t J>
    std::enable_if_t<(0 == std::tuple_size<S2>::value),
                     sutils::push_front_t<S0, std::tuple_element_t<J, typename S1::shape_t>>>
    doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
    {
      return {std::tuple_cat(std::make_tuple(shape.template shape<J>()), s.values)};
    }
  };

  template <size_t I>
  struct copy_new_axis_helper {
    template <class S0, class S1, class S2, size_t J>
    auto doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
        -> std::enable_if_t<(I < std::tuple_size<S2>::value) &&
                                safe_tuple_element<I, S2>::type::value,
                            decltype(copy_new_axis_helper<I - 1>{}.doit(
                                sutils::push_front_t<S0, std::integral_constant<long, 1>>(), shape,
                                new_axis, std::integral_constant<size_t, J>()))>
    {
      return copy_new_axis_helper<I - 1>{}.doit(
          sutils::push_front_t<S0, std::integral_constant<long, 1>>(
              std::tuple_cat(std::tuple<std::integral_constant<long, 1>>(), s.values)),
          shape, new_axis, std::integral_constant<size_t, J>());
    }

    template <class S0, class S1, class S2, size_t J>
    auto doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
        -> std::enable_if_t<
            (I >= std::tuple_size<S2>::value),
            decltype(copy_new_axis_helper<I - 1>{}.doit(
                sutils::push_front_t<S0, std::tuple_element_t<J, typename S1::shape_t>>(), shape,
                new_axis, std::integral_constant<size_t, J == 0 ? J : J - 1>()))>
    {
      return copy_new_axis_helper<I - 1>{}.doit(
          sutils::push_front_t<S0, std::tuple_element_t<J, typename S1::shape_t>>(
              std::tuple_cat(std::make_tuple(shape.template shape<J>()), s.values)),
          shape, new_axis, std::integral_constant<size_t, J == 0 ? J : J - 1>());
    }
    template <class S0, class S1, class S2, size_t J>
    auto doit(S0 s, S1 const &shape, S2 const &new_axis, std::integral_constant<size_t, J>)
        -> std::enable_if_t<
            (I < std::tuple_size<S2>::value) && !safe_tuple_element<I, S2>::type::value,
            decltype(copy_new_axis_helper<I - 1>{}.doit(
                sutils::push_front_t<S0, std::tuple_element_t<J, typename S1::shape_t>>(), shape,
                new_axis, std::integral_constant<size_t, J == 0 ? J : J - 1>()))>
    {
      return copy_new_axis_helper<I - 1>{}.doit(
          sutils::push_front_t<S0, std::tuple_element_t<J, typename S1::shape_t>>(
              std::tuple_cat(std::make_tuple(shape.template shape<J>()), s.values)),
          shape, new_axis, std::integral_constant<size_t, J == 0 ? J : J - 1>());
    }
  };

  template <size_t N, class S1, class S2>
  auto copy_new_axis(S1 const &shape, S2 const &new_axis)
      -> decltype(copy_new_axis_helper<N - 1>{}.doit(
          types::pshape<>(), shape, new_axis, std::integral_constant<size_t, S1::value - 1>()))
  {
    return copy_new_axis_helper<N - 1>{}.doit(types::pshape<>(), shape, new_axis,
                                              std::integral_constant<size_t, S1::value - 1>());
  }
} // namespace sutils

namespace types
{
  namespace details
  {
    template <class E, class S>
    void init_shape(S &res, E const &e, utils::int_<1>)
    {
      sutils::assign(std::get<std::tuple_size<S>::value - 1>(res), e.size());
    }
    template <class E, class S, size_t L>
    void init_shape(S &res, E const &e, utils::int_<L>)
    {
      sutils::assign(std::get<std::tuple_size<S>::value - L>(res), e.size());
      init_shape(res, e[0], utils::int_<L - 1>{});
    }
  } // namespace details
  template <class T, class... Tys>
  bool operator==(T const &self, pshape<Tys...> const &other)
  {
    return sutils::equals(self, other);
  }
  template <class T, class... Tys>
  bool operator==(pshape<Tys...> const &self, T const &other)
  {
    return sutils::equals(self, other);
  }
  template <class... Ty0s, class... Ty1s>
  bool operator==(pshape<Ty0s...> const &self, pshape<Ty1s...> const &other)
  {
    return sutils::equals(self, other);
  }
  template <class T, class... Tys>
  bool operator!=(T const &self, pshape<Tys...> const &other)
  {
    return !sutils::equals(self, other);
  }
  template <class T, class... Tys>
  bool operator!=(pshape<Tys...> const &self, T const &other)
  {
    return !sutils::equals(self, other);
  }
  template <class... Ty0s, class... Ty1s>
  bool operator!=(pshape<Ty0s...> const &self, pshape<Ty1s...> const &other)
  {
    return !sutils::equals(self, other);
  }
} // namespace types

PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/include/utils/fwd.hpp"
#include "pythonic/include/utils/seq.hpp"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename K, typename V>
struct to_python<std::pair<K, V>> {
  static PyObject *convert(std::pair<K, V> const &t);
};

template <typename... Tys>
struct to_python<types::pshape<Tys...>> {
  static PyObject *convert(types::pshape<Tys...> const &t);
};

template <typename... Types>
struct to_python<std::tuple<Types...>> {

  template <size_t... S>
  static PyObject *do_convert(std::tuple<Types...> const &t, std::index_sequence<S...>);

  static PyObject *convert(std::tuple<Types...> const &t);
};

template <typename T, size_t N>
struct to_python<types::array_tuple<T, N>> {
  template <size_t... S>
  static PyObject *do_convert(types::array_tuple<T, N> const &t, std::index_sequence<S...>);

  static PyObject *convert(types::array_tuple<T, N> const &t);
};

template <typename T, size_t N>
struct to_python<types::static_list<T, N>> {
  template <size_t... S>
  static PyObject *do_convert(types::static_list<T, N> const &t, std::index_sequence<S...>);

  static PyObject *convert(types::static_list<T, N> const &t);
};

template <typename... Types>
struct from_python<std::tuple<Types...>> {

  template <size_t... S>
  static bool do_is_convertible(PyObject *obj, typename std::index_sequence<S...>);

  static bool is_convertible(PyObject *obj);

  template <size_t... S>
  static std::tuple<Types...> do_convert(PyObject *obj, typename std::index_sequence<S...>);
  static std::tuple<Types...> convert(PyObject *obj);
};

template <typename T, size_t N>
struct from_python<types::array_tuple<T, N>> {

  static bool is_convertible(PyObject *obj);

  template <size_t... S>
  static types::array_tuple<T, N> do_convert(PyObject *obj, typename std::index_sequence<S...>);
  static types::array_tuple<T, N> convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
