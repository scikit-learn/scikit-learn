#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_BROADCAST_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_BROADCAST_HPP

#ifdef USE_XSIMD
#include <xsimd/xsimd.hpp>
#endif

#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/types/vectorizable_type.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct broadcasted_iterator
      : std::iterator<std::random_access_iterator_tag, std::remove_reference_t<T>> {
    T value_;

    broadcasted_iterator(T const &value) : value_(value)
    {
    }

    T const &operator*() const
    {
      return value_;
    }

    broadcasted_iterator &operator++()
    {
      return *this;
    }

    broadcasted_iterator &operator+=(long i)
    {
      return *this;
    }

    long operator-(broadcasted_iterator const &other) const
    {
      return 0;
    }

    bool operator!=(broadcasted_iterator const &other) const
    {
      return false;
    }
    bool operator==(broadcasted_iterator const &other) const
    {
      return true;
    }
    bool operator<(broadcasted_iterator const &other) const
    {
      return false;
    }
  };

  /* Type adaptor for broadcasted array values
   *
   * Used when the args of a binary operator do not have the same dimensions:
   * in that case their first dimension always yields a copy
   */
  template <class T>
  struct broadcasted {

    template <class Tp, bool is_dtype>
    struct broadcast_or_broadcasted {
      using type = broadcasted<Tp>;
    };
    template <class Tp>
    struct broadcast_or_broadcasted<Tp, true> {
      using type = broadcast<Tp, Tp>;
    };
    template <class Tp>
    using broadcast_or_broadcasted_t =
        typename broadcast_or_broadcasted<Tp, is_dtype<Tp>::value>::type;
    static const bool is_vectorizable = true;
    static const bool is_flat = false;
    static const bool is_strided = false;
    using dtype = typename std::remove_reference_t<T>::dtype;
    using value_type = typename std::remove_reference_t<T>::value_type;
    static constexpr size_t value = std::remove_reference_t<T>::value + 1;
    using const_iterator = broadcasted_iterator<T>;
    using iterator = const_iterator;

    T ref;
    using shape_t = types::array_tuple<long, value>;

    template <size_t I>
    long shape() const
    {
      return I == 0 ? 1 : (long)(ref.template shape<I == 0 ? 0 : (I - 1)>());
    }

    broadcasted() = default;

    template <class E>
    broadcasted(E &&other) : ref(std::forward<E>(other))
    {
    }

    const_iterator begin() const
    {
      return {ref};
    }
    const_iterator end() const
    {
      return {ref};
    }

    T const &operator[](long i) const;

    template <class S>
    std::enable_if_t<is_slice<S>::value, broadcasted const &> operator[](S s) const
    {
      return *this;
    }
    T const &fast(long i) const;
    template <class... Indices>
    dtype load(long i, Indices... indices) const
    {
      return ref.load(indices...);
    }

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator_nostep<broadcasted>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    template <class S>
    std::enable_if_t<is_slice<S>::value, broadcasted const &> operator()(S s) const
    {
      return *this;
    }

    T operator()(long) const
    {
      return ref;
    }

    template <class Arg1, class... Args>
    auto operator()(long arg0, Arg1 &&arg1, Args &&...args) const
        -> decltype(ref(std::forward<Arg1>(arg1), std::forward<Args>(args)...));

    template <class S, class Arg1, class... Args>
    auto operator()(S arg0, Arg1 &&arg1, Args &&...args) const -> broadcast_or_broadcasted_t<
        std::decay_t<decltype(ref(std::forward<Arg1>(arg1), std::forward<Args>(args)...))>>;

    long flat_size() const;
  };

  /* Type adaptor for scalar values
   *
   * Have them behave like infinite arrays of that value
   *
   * B is the original type of the broadcast value, and T is the type of the
   *expression it is combined with
   * if both B and T are integer types, we choose T instead of B to prevent
   * automatic conversion into larger types
   *
   * That way, np.ones(10, dtype=np.uint8) + 1 yields an array of np.uint8,
   * although 1 is of type long
   */
  template <class dtype, bool is_vectorizable>
  struct broadcast_base {
    dtype _value;
    struct ignored {
    } _splated;
    broadcast_base() = default;
    template <class V>
    broadcast_base(V v);
    template <class I>
    void load(I) const;
  };

#ifdef USE_XSIMD
  template <class dtype>
  struct broadcast_base<dtype, true> {
    dtype _value;
    xsimd::batch<dtype> _splated;
    broadcast_base() = default;

    template <class V>
    broadcast_base(V v);
    template <class I>
    auto load(I) const -> decltype(this->_splated);
  };
#endif

  template <class T>
  struct const_broadcast_iterator : public std::iterator<std::random_access_iterator_tag, T> {
    T value;
    const_broadcast_iterator(T data) : value{data}
    {
    }

    T operator*() const
    {
      return value;
    }
    const_broadcast_iterator &operator++()
    {
      return *this;
    }
    const_broadcast_iterator &operator--()
    {
      return *this;
    }
    const_broadcast_iterator &operator+=(long i)
    {
      return *this;
    }
    const_broadcast_iterator &operator-=(long i)
    {
      return *this;
    }
    const_broadcast_iterator operator+(long i) const
    {
      return *this;
    }
    const_broadcast_iterator operator-(long i) const
    {
      return *this;
    }
    long operator-(const_broadcast_iterator const &other) const
    {
      return 0;
    }
    bool operator!=(const_broadcast_iterator const &other) const
    {
      return false;
    }
    bool operator==(const_broadcast_iterator const &other) const
    {
      return true;
    }
    bool operator<(const_broadcast_iterator const &other) const
    {
      return false;
    }
    const_broadcast_iterator &operator=(const_broadcast_iterator const &other)
    {
      return *this;
    }
  };

  template <class T, class B>
  struct broadcast_dtype {
    using type = std::conditional_t<(std::is_integral<T>::value && std::is_integral<B>::value) ||
                                        (std::is_floating_point<T>::value &&
                                         std::is_floating_point<B>::value),
                                    T, typename __combined<T, B>::type>;
  };
#ifndef USE_XSIMD
  template <class T, class B>
  struct broadcast_dtype<std::complex<T>, B> {
    using type = T;
  };
  template <class T0, class T1>
  struct broadcast_dtype<std::complex<T0>, std::complex<T1>> {
    using type = std::complex<typename __combined<T0, T1>::type>;
  };
#endif

  template <class T, class B>
  struct broadcast {
    // Perform the type conversion here if it seems valid (although it is not
    // always)
    using dtype = typename broadcast_dtype<T, B>::type;
    static const bool is_vectorizable = types::is_vectorizable<dtype>::value;
    static const bool is_flat = false;
    static const bool is_strided = false;
    using value_type = dtype;
    using const_iterator = const_broadcast_iterator<dtype>;
    using iterator = const_iterator;
    static constexpr size_t value = 1;

    broadcast_base<dtype, is_vectorizable> _base;
    operator dtype() const
    {
      return _base._value;
    }

    broadcast() = default;
    template <class V>
    broadcast(V v);

    dtype operator[](long) const;
    template <size_t N>
    dtype operator[](array_tuple<long, N>) const;

    template <class S>
    std::enable_if_t<is_slice<S>::value, broadcast const &> operator[](S) const
    {
      return *this;
    }
    dtype fast(long) const;

    template <class... Indices>
    dtype load(long i, Indices... indices) const
    {
      return _base._value;
    }

    template <class I>
    auto load(I i) const -> decltype(this->_base.load(i));
    template <class... Args>
    dtype operator()(Args &&...) const;
    using shape_t = types::pshape<std::integral_constant<long, 1>>;
    template <size_t I>
    std::integral_constant<long, 1> shape() const;
    long flat_size() const;
    const_iterator begin() const
    {
      return {_base._value};
    }
    const_iterator end() const
    {
      return {_base._value};
    }
#ifdef USE_XSIMD
    using simd_iterator = const_broadcast_iterator<decltype(_base._splated)>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const
    {
      return {_base._splated};
    }
    template <class vectorizer>
    simd_iterator vend(vectorizer) const
    {
      return {_base._splated};
    }
#endif
  };
} // namespace types
PYTHONIC_NS_END

#endif
