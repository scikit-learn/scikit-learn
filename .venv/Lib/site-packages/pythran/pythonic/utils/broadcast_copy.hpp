#ifndef PYTHONIC_UTILS_BROADCAST_COPY_HPP
#define PYTHONIC_UTILS_BROADCAST_COPY_HPP

#include "pythonic/include/utils/broadcast_copy.hpp"

#include "pythonic/types/tuple.hpp"

PYTHONIC_NS_BEGIN

namespace utils
{

  /* helper for specialization of the broadcasting, vectorizing copy operator
   * due to expression templates, this may also triggers a lot of
   *computations!
   *
   * ``vector_form'' is set to true if the operation can be done using
   *Boost.SIMD
   *
   * the call operator has four template parameters:
   *
   * template <class E, class F, size_t N>
   * void operator()(E &&self, F const &other, utils::int_<N>, utils::int_<M>)
   *
   * ``E'' is the type of the object to which the data are copied
   *
   * ``F'' is the type of the object from which the data are copied
   *
   * ``N'' is the depth of the loop nest. When it reaches ``1'', we have a raw
   *loop
   *       that may be vectorizable
   *
   * ``D'' is the delta between the number of dimensions of E && F. When set
   *to a
   *       value greater than ``0'', some broadcasting is needed
   */

  template <typename vector_form, size_t N, size_t D>
  struct _broadcast_copy;

  struct fast_novectorize {
  };

  template <>
  struct _broadcast_copy<fast_novectorize, 0, 0> {
    template <class E, class F, class SelfIndices, class OtherIndices, size_t... Is>
    void helper(E &&self, F const &other, SelfIndices &&self_indices, OtherIndices &&other_indices,
                std::index_sequence<Is...>)
    {
      std::forward<E>(self).store(
          (typename std::decay_t<E>::dtype)other.load((long)std::get<Is>(other_indices)...),
          (long)std::get<Is>(self_indices)...);
    }
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }

    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      helper(std::forward<E>(self), other, self_indices, other_indices,
             std::make_index_sequence<std::tuple_size<std::decay_t<SelfIndices>>::value>());
    }
  };
  template <size_t N>
  struct _broadcast_copy<fast_novectorize, N, 0> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }
    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      long const other_size = other.template shape<std::decay<E>::type::value - N>();
      long const self_size = self.template shape<std::decay<E>::type::value - N>();
      if (self_size == other_size)
        for (long i = 0; i < self_size; ++i)
          _broadcast_copy<fast_novectorize, N - 1, 0>{}(
              std::forward<E>(self), other, std::tuple_cat(self_indices, std::make_tuple(i)),
              std::tuple_cat(other_indices, std::make_tuple(i)));
      else
        for (long i = 0; i < self_size; ++i)
          _broadcast_copy<fast_novectorize, N - 1, 0>{}(
              std::forward<E>(self), other, std::tuple_cat(self_indices, std::make_tuple(i)),
              std::tuple_cat(other_indices, std::make_tuple(0)));
    }
  };

  template <size_t N, size_t D>
  struct _broadcast_copy<fast_novectorize, N, D> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }
    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      using broadcaster = std::conditional_t<types::is_dtype<F>::value,
                                             types::broadcast<F, typename std::decay_t<E>::dtype>,
                                             types::broadcasted<F>>;
      _broadcast_copy<fast_novectorize, N, D - 1>{}(std::forward<E>(self), broadcaster(other),
                                                    std::forward<SelfIndices>(self_indices),
                                                    std::forward<OtherIndices>(other_indices));
    }
  };

  template <size_t N, class vectorizer>
  struct _broadcast_copy<vectorizer, N, 0> {
    template <class E, class F, class... Indices>
    void operator()(E &&self, F const &other, Indices... indices)
    {
      long self_size = self.size(), other_size = other.size();
      std::copy(other.begin(), other.end(), self.begin());

      // eventually repeat the pattern
      for (long i = other_size; i < self_size; i += other_size)
        std::copy_n(self.begin(), other_size, self.begin() + i);
    }
  };

  // ``D'' is not ``0'' so we should broadcast
  template <class vectorizer, size_t N, size_t D>
  struct _broadcast_copy {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      if (types::is_dtype<F>::value) {
        std::fill(self.begin(), self.end(), other);
      } else {
        auto sfirst = self.begin();
        *sfirst = other;
        std::fill(self.begin() + 1, self.end(), *sfirst);
      }
    }
    template <class E, class F, class ES, class FS>
    void operator()(E &&self, F const &other, ES, FS)
    {
      if (types::is_dtype<F>::value) {
        std::fill(self.begin(), self.end(), other);
      } else {
        auto sfirst = self.begin();
        *sfirst = other;
        std::fill(self.begin() + 1, self.end(), *sfirst);
      }
    }
  };

#ifdef USE_XSIMD
  // specialize for SIMD only if available
  // otherwise use the std::copy fallback
  template <class vectorizer, class E, class F>
  void vbroadcast_copy(E &&self, F const &other)
  {
    using T = typename F::dtype;
    using vT = xsimd::batch<T>;

    static const std::size_t vN = vT::size;

    auto oiter = vectorizer::vbegin(other);
    const long other_size = other.size();
    const long vbound = other_size / vN;

    for (auto iter = vectorizer::vbegin(self), end = iter + vbound; iter != end; ++iter, ++oiter) {
      iter.store(*oiter);
    }

    // tail
    const long bound = vbound * vN;
    if (other_size != bound) {
      auto siter = self.begin() + bound;
      auto oiter = other.begin() + bound;
      for (auto oend = other.end(); oiter < oend; ++siter, ++oiter) {
        *siter = *oiter;
      }
    }

    if (other_size != self.size())
      for (auto siter = self.begin(), send = self.end(); siter != send;)
        siter = std::copy_n(self.begin(), other_size, siter);
  }

  template <>
  struct _broadcast_copy<types::vectorizer, 1, 0> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return vbroadcast_copy<types::vectorizer>(std::forward<E>(self), other);
    }
  };
  template <>
  struct _broadcast_copy<types::vectorizer_nobroadcast, 1, 0> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return vbroadcast_copy<types::vectorizer_nobroadcast>(std::forward<E>(self), other);
    }
  };

#endif
  template <class E, class F, size_t N, size_t D, bool vector_form>
  struct broadcast_copy_dispatcher;

  template <class E, class F, size_t N, size_t D>
  struct broadcast_copy_dispatcher<E, F, N, D, false> {
    void operator()(E &self, F const &other)
    {
      if (utils::no_broadcast_ex(other))
        _broadcast_copy<fast_novectorize, N, D>{}(self, other);
      else
        _broadcast_copy<types::novectorize, N, D>{}(self, other);
    }
  };
  template <class E, class F, size_t N, size_t D>
  struct broadcast_copy_dispatcher<E, F, N, D, true> {
    void operator()(E &self, F const &other)
    {
      if (utils::no_broadcast_vectorize(other)) {
        _broadcast_copy<types::vectorizer_nobroadcast, N, D>{}(self, other);
      } else
        _broadcast_copy<types::vectorizer, N, D>{}(self, other);
    }
  };

  template <class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_copy_helper(E &self, F const &other, std::integral_constant<bool, true>,
                           std::integral_constant<bool, false>)
  {
    static_assert(D >= 0, "downcasting already happened");
    if (self.size()) {
#ifdef USE_XSIMD
      constexpr bool vectorize = vector_form;
#else
      constexpr bool vectorize = false;
      ;
#endif
      broadcast_copy_dispatcher<E, F, N, (size_t)D, vectorize>{}(self, other);
    }
    return self;
  }

  template <class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_copy_helper(E &self, F const &other, std::integral_constant<bool, true>,
                           std::integral_constant<bool, true>)
  {
    if (D == 0) {
      std::copy(other.data(), other.data() + other.flat_size(), self.data());
      return self;
    } else {
      return broadcast_copy_helper<E, F, N, D, vector_form>(
          self, other, std::integral_constant<bool, true>(), std::integral_constant<bool, false>{});
    }
  }

  template <class E, class F, size_t N, int D, bool vector_form, bool plain>
  E &broadcast_copy_helper(E &self, F const &other, std::integral_constant<bool, false>,
                           std::integral_constant<bool, plain> is_plain)
  {
    auto reshaped = other.reshape(sutils::getshape(self));
    return broadcast_copy_helper<E, decltype(reshaped), N, 0, vector_form>(
        self, reshaped, std::true_type(), is_plain);
  }

  template <class T, bool = types::is_dtype<T>::value>
  struct is_flat {
    static const bool value = T::is_flat;
  };
  template <class T>
  struct is_flat<T, true> {
    static const bool value = false;
  };

  template <class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_copy(E &self, F const &other)
  {
    return broadcast_copy_helper<E, F, N, D, vector_form>(
        self, other, std::integral_constant<bool, (D >= 0)>(), std::integral_constant < bool,
        std::decay<E>::type::is_flat &&is_flat<std::decay_t<F>>::value > {});
  }

  /* update
   */
  // ``D'' is not ``0'' so we should broadcast
  template <class Op, typename vector_form, size_t N, size_t D>
  struct _broadcast_update {

    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      long n = self.template shape<0>();
      auto siter = self.begin();
      for (long i = 0; i < n; ++i)
        Op{}(*(siter + i), other);
    }
  };

  template <class Op, size_t N, class vector_form>
  struct _broadcast_update<Op, vector_form, N, 0> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      long other_size = other.size();
      auto siter = self.begin();
      auto oiter = other.begin();
      if (other_size == 1) {
        auto value = *oiter;
        for (auto send = self.end(); siter != send; ++siter)
          Op{}(*siter, value);
      } else
        for (auto send = self.end(); siter != send;) {
          auto ooiter = oiter;
          for (long i = 0; i < other_size; ++i, ++siter, ++ooiter)
            Op{}(*siter, *ooiter);
        }
    }

    template <class E, class F0, class F1>
    void operator()(E &&self, types::broadcast<F0, F1> const &other)
    {
      auto value = *other.begin();
      for (auto siter = self.begin(), send = self.end(); siter != send; ++siter)
        Op{}(*siter, value);
    }

    template <class E, class F>
    void operator()(E &&self, types::broadcasted<F> const &other)
    {
      auto value = *other.end();
      for (auto siter = self.begin(), send = self.end(); siter != send; ++siter)
        Op{}(*siter, value);
    }
  };

  template <class Op>
  struct _broadcast_update<Op, fast_novectorize, 0, 0> {
    template <class E, class F, class SelfIndices, class OtherIndices, size_t... Is>
    void helper(E &&self, F const &other, SelfIndices &&self_indices, OtherIndices &&other_indices,
                std::index_sequence<Is...>)
    {
      self.template update<Op>(other.load((long)std::get<Is>(other_indices)...),
                               (long)std::get<Is>(self_indices)...);
    }
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }

    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      helper(std::forward<E>(self), other, self_indices, other_indices,
             std::make_index_sequence<std::tuple_size<std::decay_t<SelfIndices>>::value>());
    }
  };

  template <class Op, size_t N>
  struct _broadcast_update<Op, fast_novectorize, N, 0> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }
    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      auto const other_size = other.template shape<std::decay<E>::type::value - N>();
      auto const self_size = self.template shape<std::decay<E>::type::value - N>();
      if (self_size == other_size)
        for (long i = 0; i < self_size; ++i)
          _broadcast_update<Op, fast_novectorize, N - 1, 0>{}(
              std::forward<E>(self), other, std::tuple_cat(self_indices, std::make_tuple(i)),
              std::tuple_cat(other_indices, std::make_tuple(i)));
      else
        for (long i = 0; i < self_size; ++i)
          _broadcast_update<Op, fast_novectorize, N - 1, 0>{}(
              std::forward<E>(self), other, std::tuple_cat(self_indices, std::make_tuple(i)),
              std::tuple_cat(other_indices, std::make_tuple(0)));
    }
  };
  template <class Op, size_t N, size_t D>
  struct _broadcast_update<Op, fast_novectorize, N, D> {
    template <class E, class F>
    void operator()(E &&self, F const &other)
    {
      return (*this)(std::forward<E>(self), other, std::tuple<>(), std::tuple<>());
    }
    template <class E, class F, class SelfIndices, class OtherIndices>
    void operator()(E &&self, F const &other, SelfIndices &&self_indices,
                    OtherIndices &&other_indices)
    {
      using broadcaster = std::conditional_t<types::is_dtype<F>::value,
                                             types::broadcast<F, typename std::decay_t<E>::dtype>,
                                             types::broadcasted<F>>;
      _broadcast_update<Op, fast_novectorize, N, D - 1>{}(
          std::forward<E>(self), broadcaster(other), std::forward<SelfIndices>(self_indices),
          std::forward<OtherIndices>(other_indices));
    }
  };

#ifdef USE_XSIMD
  // specialize for SIMD only if available
  // otherwise use the std::copy fallback
  template <class Op, class vectorizer, class E, class F>
  void vbroadcast_update(E &&self, F const &other)
  {
    using T = typename F::dtype;
    using vT = typename xsimd::batch<T>;
    long other_size = other.size();

    static const std::size_t vN = vT::size;
    auto oiter = vectorizer::vbegin(other);
    auto iter = vectorizer::vbegin(self);
    const long bound = other.size() / vN * vN;

    for (auto end = vectorizer::vend(self); iter != end; ++iter, ++oiter) {
      iter.store(Op{}(*iter, *oiter));
    }
    // tail
    {
      auto siter = self.begin();
      auto oiter = other.begin();
      for (long i = bound; i < other_size; ++i)
        Op{}(*(siter + i), *(oiter + i));
    }
  }

  template <class Op, class vectorizer, class E, class F0, class F1>
  void vbroadcast_update(E &&self, types::broadcast<F0, F1> const &other)
  {
    auto value = *other.begin();
    for (auto siter = self.begin(), send = self.end(); siter != send; ++siter)
      Op{}(*siter, value);
  }

  template <class Op, class vectorizer, class E, class F>
  void vbroadcast_update(E &&self, types::broadcasted<F> const &other)
  {
    auto value = *other.end();
    for (auto siter = self.begin(), send = self.end(); siter != send; ++siter)
      Op{}(*siter, value);
  }

  template <class Op>
  struct _broadcast_update<Op, types::vectorizer, 1, 0> {
    template <class... Args>
    void operator()(Args &&...args)
    {
      vbroadcast_update<Op, types::vectorizer>(std::forward<Args>(args)...);
    }
  };
  template <class Op>
  struct _broadcast_update<Op, types::vectorizer_nobroadcast, 1, 0> {
    template <class... Args>
    void operator()(Args &&...args)
    {
      vbroadcast_update<Op, types::vectorizer_nobroadcast>(std::forward<Args>(args)...);
    }
  };

#endif
  template <class Op, bool vector_form, class E, class F, size_t N, size_t D>
  struct broadcast_update_dispatcher;

  template <class Op, class E, class F, size_t N, size_t D>
  struct broadcast_update_dispatcher<Op, false, E, F, N, D> {
    void operator()(E &self, F const &other)
    {
      if (utils::no_broadcast_ex(other))
        _broadcast_update<Op, fast_novectorize, N, D>{}(self, other);
      else
        _broadcast_update<Op, types::novectorize, N, D>{}(self, other);
    }
  };
  template <class Op, class E, class F, size_t N, size_t D>
  struct broadcast_update_dispatcher<Op, true, E, F, N, D> {
    void operator()(E &self, F const &other)
    {
      if (utils::no_broadcast_vectorize(other))
        _broadcast_update<Op, types::vectorizer_nobroadcast, N, D>{}(self, other);
      else
        _broadcast_update<Op, types::vectorizer, N, D>{}(self, other);
    }
  };

  template <class Op, class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_update(E &self, F const &other)
  {
    if (self.size())
#ifdef USE_XSIMD
      broadcast_update_dispatcher<Op, vector_form, E, F, N, D>{}(self, other);
#else
      broadcast_update_dispatcher<Op, false, E, F, N, D>{}(self, other);
#endif
    return self;
  }

} // namespace utils
PYTHONIC_NS_END

#endif
