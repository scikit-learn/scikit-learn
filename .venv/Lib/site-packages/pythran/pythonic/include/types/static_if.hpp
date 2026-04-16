#ifndef PYTHONIC_INCLUDE_TYPES_STATIC_IF_HPP
#define PYTHONIC_INCLUDE_TYPES_STATIC_IF_HPP

#include <tuple>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct StaticIfReturn {
    T arg;
  };

  template <class T>
  struct StaticIfNoReturn {
    T arg;
    StaticIfNoReturn() = default;
    StaticIfNoReturn(T const &arg) : arg(arg)
    {
    }
    long get(std::integral_constant<size_t, 0>)
    {
      return 0;
    }
    StaticIfNoReturn &get(std::integral_constant<size_t, 1>)
    {
      return *this;
    }
    T &get(std::integral_constant<size_t, 2>)
    {
      return arg;
    }

    template <class OT>
    operator OT() const
    {
      return OT();
    }
  };

  template <class T>
  struct StaticIfBreak {
    T arg;
    StaticIfBreak() = default;
    StaticIfBreak(T const &arg) : arg(arg)
    {
    }
    long get(std::integral_constant<size_t, 0>)
    {
      return 1;
    }
    T &get(std::integral_constant<size_t, 1>)
    {
      return arg;
    }
    T &get(std::integral_constant<size_t, 2>)
    {
      return arg;
    }

    template <class OT>
    operator OT() const
    {
      return OT();
    }
  };

  template <class T>
  struct StaticIfCont {
    T arg;
    StaticIfCont() = default;
    StaticIfCont(T const &arg) : arg(arg)
    {
    }
    long get(std::integral_constant<size_t, 0>)
    {
      return 2;
    }
    StaticIfCont &get(std::integral_constant<size_t, 1>)
    {
      return *this;
    }
    T &get(std::integral_constant<size_t, 2>)
    {
      return arg;
    }

    template <class OT>
    operator OT() const
    {
      return OT();
    }
  };

  template <class T0, class T1>
  struct StaticIfReturnHolder {
    std::tuple<int, T0, T1> args;

    StaticIfReturnHolder() : args(0, T0(), T1())
    {
    }
    StaticIfReturnHolder(StaticIfReturnHolder<T0, T1> const &) = default;

    template <class Tp0, class Tp1>
    StaticIfReturnHolder(StaticIfReturnHolder<Tp0, Tp1> const &other) : args(other.args)
    {
    }

    template <class Tp0>
    StaticIfReturnHolder(StaticIfReturn<Tp0> const &arg) : args(1, arg.arg, T1())
    {
    }
    StaticIfReturnHolder(StaticIfNoReturn<T1> const &arg) : args(0, T0(), arg.arg)
    {
    }
    StaticIfReturnHolder(StaticIfBreak<T1> const &arg) : args(2, T0(), arg.arg)
    {
    }
    StaticIfReturnHolder(StaticIfCont<T1> const &arg) : args(3, T0(), arg.arg)
    {
    }
  };
} // namespace types

PYTHONIC_NS_END

namespace std
{
  template <size_t I, class T0, class T1>
  struct tuple_element<I, pythonic::types::StaticIfReturnHolder<T0, T1>> {
    using type = std::conditional_t<I == 0, bool, std::conditional_t<I == 1, T0, T1>>;
  };
  template <size_t I, class T0, class T1>
  auto get(pythonic::types::StaticIfReturnHolder<T0, T1> &t) -> decltype(std::get<I>(t.args))
  {
    return std::get<I>(t.args);
  }

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::StaticIfNoReturn<T>> {
    using type = decltype(std::declval<pythonic::types::StaticIfNoReturn<T>>().get(
        std::integral_constant<size_t, I>{}));
  };

  template <size_t I, class T>
  auto get(pythonic::types::StaticIfNoReturn<T> &t)
      -> decltype(t.get(std::integral_constant<size_t, I>{}))
  {
    return t.get(std::integral_constant<size_t, I>{});
  }

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::StaticIfBreak<T>> {
    using type = decltype(std::declval<pythonic::types::StaticIfBreak<T>>().get(
        std::integral_constant<size_t, I>{}));
  };

  template <size_t I, class T>
  auto get(pythonic::types::StaticIfBreak<T> &t)
      -> decltype(t.get(std::integral_constant<size_t, I>{}))
  {
    return t.get(std::integral_constant<size_t, I>{});
  }

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::StaticIfCont<T>> {
    using type = decltype(std::declval<pythonic::types::StaticIfCont<T>>().get(
        std::integral_constant<size_t, I>{}));
  };

  template <size_t I, class T>
  auto get(pythonic::types::StaticIfCont<T> &t)
      -> decltype(t.get(std::integral_constant<size_t, I>{}))
  {
    return t.get(std::integral_constant<size_t, I>{});
  }
} // namespace std

/* type inference stuff { */
#include "pythonic/include/types/combined.hpp"

template <class T0, class T1>
struct __combined<pythonic::types::StaticIfReturn<T0>, pythonic::types::StaticIfNoReturn<T1>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};
template <class T0, class T1>
struct __combined<pythonic::types::StaticIfReturn<T0>, pythonic::types::StaticIfBreak<T1>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};
template <class T0, class T1>
struct __combined<pythonic::types::StaticIfReturn<T0>, pythonic::types::StaticIfCont<T1>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};

template <class T0, class T1>
struct __combined<pythonic::types::StaticIfNoReturn<T1>, pythonic::types::StaticIfReturn<T0>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};

template <class T0>
struct __combined<pythonic::types::StaticIfNoReturn<T0> &, pythonic::types::none_type> {
  using type = pythonic::types::none_type;
};
template <class T0>
struct __combined<pythonic::types::StaticIfNoReturn<T0>, pythonic::types::none_type> {
  using type = pythonic::types::none_type;
};
template <class T0>
struct __combined<pythonic::types::none_type, pythonic::types::StaticIfNoReturn<T0>> {
  using type = pythonic::types::none_type;
};
template <class T0, class T1>
struct __combined<pythonic::types::StaticIfBreak<T1>, pythonic::types::StaticIfReturn<T0>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};
template <class T0, class T1>
struct __combined<pythonic::types::StaticIfCont<T1>, pythonic::types::StaticIfReturn<T0>> {
  using type = pythonic::types::StaticIfReturnHolder<T0, T1>;
};

template <class T0, class T1, class T2, class T3>
struct __combined<pythonic::types::StaticIfReturnHolder<T0, T1>,
                  pythonic::types::StaticIfReturnHolder<T2, T3>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type,
                                                     typename __combined<T1, T3>::type>;
};
template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfReturnHolder<T0, T1>,
                  pythonic::types::StaticIfCont<T2>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};
template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfReturnHolder<T0, T1>,
                  pythonic::types::StaticIfBreak<T2>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};
template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfReturnHolder<T0, T1>,
                  pythonic::types::StaticIfReturn<T2>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};

template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfCont<T2>,
                  pythonic::types::StaticIfReturnHolder<T0, T1>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};

template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfBreak<T2>,
                  pythonic::types::StaticIfReturnHolder<T0, T1>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};

template <class T0, class T1, class T2>
struct __combined<pythonic::types::StaticIfReturn<T2>,
                  pythonic::types::StaticIfReturnHolder<T0, T1>> {
  using type = pythonic::types::StaticIfReturnHolder<typename __combined<T0, T2>::type, T1>;
};

template <class T0, class T1>
struct __combined<T0, pythonic::types::StaticIfNoReturn<T1>> {
  using type = T0;
};

template <class T0, class T1>
struct __combined<pythonic::types::StaticIfNoReturn<T1>, T0> {
  using type = T0;
};

/* } */
#endif
