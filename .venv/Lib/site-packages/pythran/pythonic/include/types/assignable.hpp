#ifndef PYTHONIC_INCLUDE_TYPES_ASSIGNABLE_HPP
#define PYTHONIC_INCLUDE_TYPES_ASSIGNABLE_HPP

#include <type_traits>
#include <utility>

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T>
  constexpr T as_const(T &&t) noexcept
  {
    return std::forward<T>(t);
  }

  // Pass all scalars by value when called through pythonic::types::call
  template <class T, bool is_integral>
  struct by_val {
    using type = T;
  };
  template <class T>
  struct by_val<T &, true> {
    using type = T;
  };
  template <class T>
  struct by_val<T &&, true> {
    using type = T;
  };
  template <class T>
  struct by_val<T const &, true> {
    using type = T;
  };

  template <class T>
  using by_val_t = typename by_val<T, std::is_integral<std::decay_t<T>>::value>::type;

  template <class F, class... Args>
  static inline auto call(F &&f, Args &&...args)
      -> decltype(std::forward<F>(f).template operator()<by_val_t<Args>...>(
          static_cast<by_val_t<Args>>(args)...))
  {
    return std::forward<F>(f).template operator()<by_val_t<Args>...>(
        static_cast<by_val_t<Args>>(args)...);
  }

} // namespace types

struct dummy {
};

template <class T>
struct assignable {
  using type = T;
};

template <class T>
struct assignable<T const> : assignable<T> {
};

template <class T>
struct assignable<T const &> : assignable<T> {
};

template <class T>
struct assignable<T &> : assignable<T> {
};

template <class T>
struct assignable<T &&> : assignable<T> {
};

template <class T>
struct lazy : assignable<T> {
}; // very conservative

template <class T>
struct assignable_noescape : assignable<T> {
};

template <class T>
struct assignable_noescape<T const> : assignable_noescape<T> {
};

template <class T>
struct assignable_noescape<T const &> : assignable_noescape<T> {
};

template <class T>
struct assignable_noescape<T &> : assignable_noescape<T> {
};

template <class T>
struct assignable_noescape<T &&> : assignable_noescape<T> {
};

template <class T>
struct returnable : assignable<T> {
};

template <class T>
struct returnable<T const &> : assignable<typename returnable<T>::type> {
};

template <class T>
struct returnable<T &> : assignable<typename returnable<T>::type> {
};

template <class T>
struct returnable<T &&> : assignable<typename returnable<T>::type> {
};
PYTHONIC_NS_END

#endif
