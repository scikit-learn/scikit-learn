#ifndef PYTHONIC_INCLUDE_TYPES_VECTORIZABLE_TYPE_HPP
#define PYTHONIC_INCLUDE_TYPES_VECTORIZABLE_TYPE_HPP

PYTHONIC_NS_BEGIN
namespace types
{
  /* types used during vectorization specialization
   */
  struct vectorize {
  };
  struct novectorize {
  };
  struct novectorize_nobroadcast {
  };
  struct vectorizer {
    template <class E>
    static auto vbegin(E &&expr) -> decltype(std::forward<E>(expr).vbegin(vectorize{}))
    {
      return std::forward<E>(expr).vbegin(vectorize{});
    }
    template <class E>
    static auto vend(E &&expr) -> decltype(std::forward<E>(expr).vend(vectorize{}))
    {
      return std::forward<E>(expr).vend(vectorize{});
    }
  };
  struct vectorize_nobroadcast {
  };
  struct vectorizer_nobroadcast {
    template <class E>
    static auto vbegin(E &&expr) -> decltype(std::forward<E>(expr).vbegin(vectorize_nobroadcast{}))
    {
      return std::forward<E>(expr).vbegin(vectorize_nobroadcast{});
    }
    template <class E>
    static auto vend(E &&expr) -> decltype(std::forward<E>(expr).vend(vectorize_nobroadcast{}))
    {
      return std::forward<E>(expr).vend(vectorize_nobroadcast{});
    }
  };

  template <class T>
  struct is_vectorizable_dtype {
    static const bool value = is_dtype<T>::value && !std::is_same<T, bool>::value &&
                              !std::is_same<T, long double>::value &&
                              !std::is_same<T, std::complex<long double>>::value;
  };

  /* trait to check if is T is an array-like type that supports vectorization
   */
  template <class T, bool scalar = has_vectorizable<T>::value>
  struct is_vectorizable_array;

  template <class T>
  struct is_vectorizable_array<T, false> : std::false_type {
  };

  template <class T>
  struct is_vectorizable_array<T, true> {
    static const bool value = T::is_vectorizable;
  };

  template <class T>
  struct is_vectorizable {
    static const bool value = std::conditional_t<is_dtype<T>::value, is_vectorizable_dtype<T>,
                                                 is_vectorizable_array<T>>::value;
  };

  template <class O, class... Args>
  struct is_vector_op;

  template <class Op, class... Args>
  struct numpy_expr;
} // namespace types

namespace utils
{
  template <class Op, class... Args>
  bool no_broadcast(types::numpy_expr<Op, Args...> const &arg)
  {
    return arg.no_broadcast();
  }
  template <class Op, class... Args>
  bool no_broadcast_ex(types::numpy_expr<Op, Args...> const &arg)
  {
    return arg.no_broadcast_ex();
  }
  template <class Op, class... Args>
  bool no_broadcast_vectorize(types::numpy_expr<Op, Args...> const &arg)
  {
    return arg.no_broadcast_vectorize();
  }
  template <class Arg>
  constexpr bool no_broadcast(Arg const &arg)
  {
    return true;
  }
  template <class Arg>
  constexpr bool no_broadcast_ex(Arg const &arg)
  {
    return true;
  }
  template <class Arg>
  constexpr bool no_broadcast_vectorize(Arg const &arg)
  {
    return true;
  }
} // namespace utils
PYTHONIC_NS_END
#endif
