#ifndef PYTHONIC_INCLUDE_UTILS_FUNCTOR_HPP
#define PYTHONIC_INCLUDE_UTILS_FUNCTOR_HPP

#include <utility>

// create a function named `name' using function `f'

#define DEFINE_FUNCTOR_2(name, f)                                                                  \
  namespace functor                                                                                \
  {                                                                                                \
    struct name {                                                                                  \
      using callable = void;                                                                       \
      template <typename... Types>                                                                 \
      auto operator()(Types &&...types) const -> decltype(f(std::forward<Types>(types)...))        \
      {                                                                                            \
        return f(std::forward<Types>(types)...);                                                   \
      }                                                                                            \
                                                                                                   \
      friend std::ostream &operator<<(std::ostream &os, name)                                      \
      {                                                                                            \
        return os << #name;                                                                        \
      }                                                                                            \
    };                                                                                             \
  }

// create a functor named `f' using function `ns::f'
#define DEFINE_FUNCTOR(ns, f) DEFINE_FUNCTOR_2(f, ns::f)

#define USING_FUNCTOR(f, alias)                                                                    \
  namespace functor                                                                                \
  {                                                                                                \
    using f = alias;                                                                               \
  }

#endif
