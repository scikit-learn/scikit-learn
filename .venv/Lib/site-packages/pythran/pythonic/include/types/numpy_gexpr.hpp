#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_GEXPR_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_GEXPR_HPP

#include "pythonic/include/types/numpy_iexpr.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/include/utils/meta.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  /* helper to count new axis
   */
  template <class... T>
  struct count_new_axis;
  template <>
  struct count_new_axis<> {
    static constexpr size_t value = 0;
  };

  template <>
  struct count_new_axis<types::none_type> {
    static constexpr size_t value = 1;
  };

  template <class T>
  struct count_new_axis<T> {
    static constexpr size_t value = 0;
  };

  template <class T0, class... T>
  struct count_new_axis<T0, T...> {
    static constexpr size_t value = count_new_axis<T0>::value + count_new_axis<T...>::value;
  };

  /* helper to cast slices
   */
  template <size_t T, size_t Tp>
  long recast_slice(long n)
  {
    return (n * (long)T) / (long)Tp;
  }
  template <size_t T, size_t Tp>
  cstride_slice<1> recast_slice(normalized_slice s)
  {
    assert(s.step == 1 && "To change to a dtype of a different size, the last "
                          "axis must be contiguous");
    return {s.lower * (long)T / (long)Tp, s.upper * (long)T / (long)Tp};
  }
  template <size_t T, size_t Tp, long stride>
  cstride_slice<stride> recast_slice(cstride_normalized_slice<stride> s)
  {
    return {s.lower * (long)T / (long)Tp, s.upper * (long)T / (long)Tp};
  }

  /* helper to turn a new axis into a slice
   */
  template <class T>
  struct to_slice {
    using type = T;
    static constexpr bool is_new_axis = false;
    T operator()(T value);
  };

  template <>
  struct to_slice<none_type> {
    using type = fast_contiguous_slice;
    static constexpr bool is_new_axis = true;
    fast_contiguous_slice operator()(none_type);
  };

  template <class T>
  struct to_normalized_slice {
    using type = T;
    T operator()(T value);
  };

  template <>
  struct to_normalized_slice<none_type> {
    using type = cstride_normalized_slice<1>;
    type operator()(none_type);
  };

  /* helper to build a new shape out of a shape and a slice with new axis
   */
  template <size_t N, class pS, class IsNewAxis>
  auto make_reshape(pS const &shape, IsNewAxis is_new_axis)
      -> decltype(sutils::copy_new_axis<pS::value + N>(shape, is_new_axis));

  /* helper to build an extended slice aka numpy_gexpr out of a subscript
   */
  template <size_t C>
  struct extended_slice {
    template <class E, class... S>
    auto operator()(E &&expr, S const &...s)
        -> decltype(std::forward<E>(expr).reshape(make_reshape<C>(
            expr, std::tuple<std::integral_constant<bool, to_slice<S>::is_new_axis>...>()))(
            to_slice<S>{}(s)...))

    {
      return std::forward<E>(expr).reshape(make_reshape<C>(
          expr, std::tuple<std::integral_constant<bool, to_slice<S>::is_new_axis>...>()))(
          to_slice<S>{}(s)...);
    }
  };

  template <>
  struct extended_slice<0> {
    template <class E, class... S>
    auto operator()(E &&expr, long const &s0, S const &...s)
        -> std::enable_if_t<utils::all_of<std::is_integral<S>::value...>::value,
                            decltype(std::forward<E>(expr)[types::make_tuple(s0, s...)])>
    {
      return std::forward<E>(expr)[types::make_tuple(s0, s...)];
    }
    template <class E, class... S>
    auto operator()(E &&expr, long const &s0, S const &...s)
        -> std::enable_if_t<!utils::all_of<std::is_integral<S>::value...>::value,
                            decltype(std::forward<E>(expr)[s0](s...))>
    {
      return std::forward<E>(expr)[s0](s...);
    }

    template <class E, class... S, size_t... Is>
    numpy_gexpr<std::decay_t<E>, normalize_t<S>...> fwd(E &&expr, std::tuple<S...> const &s,
                                                        std::index_sequence<Is...>)
    {
      return {std::forward<E>(expr), std::get<Is>(s).normalize(expr.template shape<Is>())...};
    }

    template <class E, class Sp, class... S>
    std::enable_if_t<is_slice<Sp>::value, numpy_gexpr<E, normalize_t<Sp>, normalize_t<S>...>>
    operator()(E &&expr, Sp const &s0, S const &...s)
    {
      return make_gexpr(std::forward<E>(expr), s0, s...);
    }

    template <class E, class F, class... S>
    std::enable_if_t<!is_slice<F>::value,
                     numpy_gexpr<ndarray<typename std::decay_t<E>::dtype,
                                         array_tuple<long, std::decay_t<E>::value>>,
                                 cstride_normalized_slice<1>, normalize_t<S>...>>
    operator()(E &&expr, F const &s0, S const &...s)
    {
      return numpy_vexpr<
          ndarray<typename std::decay_t<E>::dtype, array_tuple<long, std::decay_t<E>::value>>, F>{
          std::forward<E>(expr), s0}(fast_contiguous_slice(none_type{}, none_type{}), s...);
    }
  };

  /* Meta-Function to count the number of slices in a type list
   */
  template <class... Types>
  struct count_long;

  template <>
  struct count_long<long> {
    static constexpr size_t value = 1;
  };

  template <>
  struct count_long<normalized_slice> {
    static constexpr size_t value = 0;
  };

  template <long stride>
  struct count_long<cstride_normalized_slice<stride>> {
    static constexpr size_t value = 0;
  };

  template <class T, class... Types>
  struct count_long<T, Types...> {
    static constexpr size_t value = count_long<T>::value + count_long<Types...>::value;
  };

  template <>
  struct count_long<> {
    static constexpr size_t value = 0;
  };

  /* helper to get the type of the nth element of an array
   */
  template <class T, size_t N>
  struct nth_value_type {
    using type = typename nth_value_type<typename T::value_type, N - 1>::type;
  };

  template <class T>
  struct nth_value_type<T, 0> {
    using type = T;
  };

  /* helper that yields true if the first slice of a pack is a contiguous
   * slice
   */
  template <class... S>
  struct is_contiguous {
    static const bool value = false;
  };

  template <class... S>
  struct is_contiguous<cstride_normalized_slice<1>, S...> {
    static const bool value = true;
  };

  /* numpy_gexpr factory
   *
   * replaces the constructor, in order to properly merge gexpr composition
   *into a single gexpr
   */
  namespace details
  {

    template <class T, class Ts, size_t... Is>
    std::tuple<T, std::tuple_element_t<Is, Ts>...> tuple_push_head(T const &val, Ts const &vals,
                                                                   std::index_sequence<Is...>)
    {
      return std::tuple<T, std::tuple_element_t<Is, Ts>...>{val, std::get<Is>(vals)...};
    }

    template <class T, class Ts>
    auto tuple_push_head(T const &val, Ts const &vals)
        -> decltype(tuple_push_head(val, vals,
                                    std::make_index_sequence<std::tuple_size<Ts>::value>()))
    {
      return tuple_push_head(val, vals, std::make_index_sequence<std::tuple_size<Ts>::value>());
    }

    // this struct is specialized for every type combination && takes care of
    // the slice merge
    template <class T, class Tp>
    struct merge_gexpr;

    template <>
    struct merge_gexpr<std::tuple<>, std::tuple<>> {
      template <size_t I, class S>
      std::tuple<> run(S const &, std::tuple<> const &t0, std::tuple<> const &);
    };

    template <class... T0>
    struct merge_gexpr<std::tuple<T0...>, std::tuple<>> {
      template <size_t I, class S>
      std::tuple<T0...> run(S const &, std::tuple<T0...> const &t0, std::tuple<>);
      static_assert(utils::all_of<std::is_same<T0, normalize_t<T0>>::value...>::value,
                    "all slices are normalized");
    };

    template <class... T1>
    struct merge_gexpr<std::tuple<>, std::tuple<T1...>> {
      template <size_t I, class S>
      std::tuple<normalize_t<T1>...> run(S const &, std::tuple<>, std::tuple<T1...> const &t1);
    };

    template <class S0, class... T0, class S1, class... T1>
    struct merge_gexpr<std::tuple<S0, T0...>, std::tuple<S1, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<S0, T0...> const &t0, std::tuple<S1, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t0) * std::get<0>(t1),
                                      merge_gexpr<std::tuple<T0...>, std::tuple<T1...>>{}
                                          .template run<I + 1>(s, tuple_tail(t0), tuple_tail(t1))))
      {
        return tuple_push_head(
            std::get<0>(t0) * std::get<0>(t1),
            merge_gexpr<std::tuple<T0...>, std::tuple<T1...>>{}.template run<I + 1>(
                s, tuple_tail(t0), tuple_tail(t1)));
      }
      static_assert(
          std::is_same<decltype(std::declval<S0>() * std::declval<S1>()),
                       normalize_t<decltype(std::declval<S0>() * std::declval<S1>())>>::value,
          "all slices are normalized");
    };

    template <class... T1>
    struct merge_gexpr<std::tuple<>, std::tuple<none_type, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<> const &t0, std::tuple<none_type, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t1),
                                      merge_gexpr<std::tuple<>, std::tuple<T1...>>{}
                                          .template run<I + 1>(s, t0, tuple_tail(t1))))
      {
        return tuple_push_head(std::get<0>(t1),
                               merge_gexpr<std::tuple<>, std::tuple<T1...>>{}.template run<I + 1>(
                                   s, t0, tuple_tail(t1)));
      }
    };

    template <class S0, class... T0, class... T1>
    struct merge_gexpr<std::tuple<S0, T0...>, std::tuple<none_type, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<S0, T0...> const &t0, std::tuple<none_type, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t1),
                                      merge_gexpr<std::tuple<S0, T0...>, std::tuple<T1...>>{}
                                          .template run<I + 1>(s, t0, tuple_tail(t1))))
      {
        return tuple_push_head(
            std::get<0>(t1),
            merge_gexpr<std::tuple<S0, T0...>, std::tuple<T1...>>{}.template run<I + 1>(
                s, t0, tuple_tail(t1)));
      }
    };

    template <class... T0, class S1, class... T1>
    struct merge_gexpr<std::tuple<long, T0...>, std::tuple<S1, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<long, T0...> const &t0, std::tuple<S1, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t0),
                                      merge_gexpr<std::tuple<T0...>, std::tuple<S1, T1...>>{}
                                          .template run<I>(s, tuple_tail(t0), t1)))
      {
        return tuple_push_head(
            std::get<0>(t0),
            merge_gexpr<std::tuple<T0...>, std::tuple<S1, T1...>>{}.template run<I>(
                s, tuple_tail(t0), t1));
      }
    };
    template <class... T0, class... T1>
    struct merge_gexpr<std::tuple<long, T0...>, std::tuple<none_type, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<long, T0...> const &t0,
               std::tuple<none_type, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t0),
                                      merge_gexpr<std::tuple<T0...>, std::tuple<none_type, T1...>>{}
                                          .template run<I>(s, tuple_tail(t0), t1)))
      {
        return tuple_push_head(
            std::get<0>(t0),
            merge_gexpr<std::tuple<T0...>, std::tuple<none_type, T1...>>{}.template run<I>(
                s, tuple_tail(t0), t1));
      }
    };

    template <class S0, class... T0, class... T1>
    struct merge_gexpr<std::tuple<S0, T0...>, std::tuple<long, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<S0, T0...> const &t0, std::tuple<long, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t1),
                                      merge_gexpr<std::tuple<T0...>, std::tuple<T1...>>{}
                                          .template run<I + 1>(s, tuple_tail(t0), tuple_tail(t1))))
      {
        return tuple_push_head(
            std::get<0>(t1) * std::get<0>(t0).step + std::get<0>(t0).lower,
            merge_gexpr<std::tuple<T0...>, std::tuple<T1...>>{}.template run<I + 1>(
                s, tuple_tail(t0), tuple_tail(t1)));
      }
    };

    template <class... T0, class... T1>
    struct merge_gexpr<std::tuple<long, T0...>, std::tuple<long, T1...>> {
      template <size_t I, class S>
      auto run(S const &s, std::tuple<long, T0...> const &t0, std::tuple<long, T1...> const &t1)
          -> decltype(tuple_push_head(std::get<0>(t0),
                                      merge_gexpr<std::tuple<T0...>, std::tuple<long, T1...>>{}
                                          .template run<I>(s, tuple_tail(t0), t1)))
      {
        return tuple_push_head(
            std::get<0>(t0),
            merge_gexpr<std::tuple<T0...>, std::tuple<long, T1...>>{}.template run<I>(
                s, tuple_tail(t0), t1));
      }
    };

    template <class Arg, class... Sp>
    std::enable_if_t<count_new_axis<Sp...>::value == 0, numpy_gexpr<Arg, Sp...>>
    _make_gexpr(Arg arg, std::tuple<Sp...> const &t);

    template <class Arg, class S, size_t... Is>
    numpy_gexpr<Arg, typename to_normalized_slice<std::tuple_element_t<Is, S>>::type...>
    _make_gexpr_helper(Arg arg, S const &s, std::index_sequence<Is...>);

    template <class Arg, class... Sp>
    auto _make_gexpr(Arg arg, std::tuple<Sp...> const &s) -> std::enable_if_t<
        count_new_axis<Sp...>::value != 0,
        decltype(_make_gexpr_helper(
            arg.reshape(make_reshape<count_new_axis<Sp...>::value>(
                arg, std::tuple<std::integral_constant<bool, to_slice<Sp>::is_new_axis>...>())),
            s, std::make_index_sequence<sizeof...(Sp)>()))>;

    template <class Arg, class... S>
    struct make_gexpr {
      template <size_t... Is>
      numpy_gexpr<Arg, normalize_t<S>...> operator()(Arg arg, std::tuple<S...>,
                                                     std::index_sequence<Is...>);
      numpy_gexpr<Arg, normalize_t<S>...> operator()(Arg arg, S const &...s);
    };

    // this specialization is in charge of merging gexpr
    template <class Arg, class... S, class... Sp>
    struct make_gexpr<numpy_gexpr<Arg, S...> const &, Sp...> {
      auto operator()(numpy_gexpr<Arg, S...> const &arg, Sp const &...s) -> decltype(_make_gexpr(
          std::declval<Arg>(), merge_gexpr<std::tuple<S...>, std::tuple<Sp...>>{}.template run<0>(
                                   arg, std::tuple<S...>(), std::tuple<Sp...>())))
      {
        return _make_gexpr(arg.arg,
                           merge_gexpr<std::tuple<S...>, std::tuple<Sp...>>{}.template run<0>(
                               arg, arg.slices, std::make_tuple(s...)));
      }
    };
  } // namespace details

  template <class Arg, class... S>
  auto make_gexpr(Arg &&arg, S const &...s)
      -> decltype(details::make_gexpr<Arg, S...>{}(std::forward<Arg>(arg), s...));

  /* type-based compile time overlapping detection: detect if a type may
   *overlap with another
   * the goal is to detect whether the following operation
   *
   * a[...] = b
   *
   * requires a copy.
   *
   * It requires a copy if b = a[...], as in
   *
   * a[1:] = a[:-1]
   *
   * because this is *!* equivalent to for i in range(0, n-1): a[i+1] = a[i]
   *
   * to avoid the copy, we rely on the lhs type
   */

  template <class E>
  struct may_overlap_gexpr : std::integral_constant<bool, !is_dtype<E>::value> {
  };

  template <class T0, class T1>
  struct may_overlap_gexpr<broadcast<T0, T1>> : std::false_type {
  };

  template <class E>
  struct may_overlap_gexpr<broadcasted<E>> : std::false_type {
  };

  template <class E>
  struct may_overlap_gexpr<E &> : may_overlap_gexpr<E> {
  };

  template <class E>
  struct may_overlap_gexpr<E const &> : may_overlap_gexpr<E> {
  };

  template <class T, class pS>
  struct may_overlap_gexpr<ndarray<T, pS>> : std::false_type {
  };

  template <class E>
  struct may_overlap_gexpr<numpy_iexpr<E>> : may_overlap_gexpr<E> {
  };

  template <class E>
  struct may_overlap_gexpr<numpy_texpr<E>> : may_overlap_gexpr<E> {
  };

  template <class E>
  struct may_overlap_gexpr<list<E>> : std::integral_constant<bool, !is_dtype<E>::value> {
  };

  template <class E, size_t N, class V>
  struct may_overlap_gexpr<array_base<E, N, V>> : may_overlap_gexpr<E> {
  };

  template <class Op, class... Args>
  struct may_overlap_gexpr<numpy_expr<Op, Args...>>
      : utils::any_of<may_overlap_gexpr<Args>::value...> {
  };

  template <class OpS, class pS, class... S>
  struct gexpr_shape;

  template <class... Tys, class... oTys>
  struct gexpr_shape<pshape<Tys...>, pshape<oTys...>> {
    using type = pshape<Tys..., oTys...>;
  };

  template <class... Tys>
  struct gexpr_shape<pshape<Tys...>, array_tuple<long, 0>> {
    using type = pshape<Tys...>;
  };

  template <class... Tys, size_t N>
  struct gexpr_shape<pshape<Tys...>, array_tuple<long, N>>
      : gexpr_shape<pshape<Tys..., long>, array_tuple<long, N - 1>> {
  };

  template <class... Tys, class... oTys, class... S, long stride>
  struct gexpr_shape<pshape<Tys...>, pshape<std::integral_constant<long, 1>, oTys...>,
                     cstride_normalized_slice<stride>, S...>
      : gexpr_shape<pshape<Tys..., std::integral_constant<long, 1>>, pshape<oTys...>, S...> {
  };
  template <class... Tys, class... oTys, class... S>
  struct gexpr_shape<pshape<Tys...>, pshape<std::integral_constant<long, 1>, oTys...>,
                     normalized_slice, S...>
      : gexpr_shape<pshape<Tys..., std::integral_constant<long, 1>>, pshape<oTys...>, S...> {
  };

  template <class... Tys, class oT, class... oTys, class... S>
  struct gexpr_shape<pshape<Tys...>, pshape<oT, oTys...>, long, S...>
      : gexpr_shape<pshape<Tys...>, pshape<oTys...>, S...> {
  };
  template <class... Tys, class oT, class... oTys, class cS, class... S>
  struct gexpr_shape<pshape<Tys...>, pshape<oT, oTys...>, cS, S...>
      : gexpr_shape<pshape<Tys..., long>, pshape<oTys...>, S...> {
  };
  template <class... Tys, size_t N, class... S>
  struct gexpr_shape<pshape<Tys...>, array_tuple<long, N>, long, S...>
      : gexpr_shape<pshape<Tys...>, array_tuple<long, N - 1>, S...> {
  };
  template <class... Tys, size_t N, class cS, class... S>
  struct gexpr_shape<pshape<Tys...>, array_tuple<long, N>, cS, S...>
      : gexpr_shape<pshape<Tys..., long>, array_tuple<long, N - 1>, S...> {
  };

  template <class pS, class... S>
  using gexpr_shape_t = typename gexpr_shape<pshape<>, pS, S...>::type;

  /* Expression template for numpy expressions - extended slicing operators
   */
  template <class Arg, class... S>
  struct numpy_gexpr {
    static_assert(utils::all_of<std::is_same<S, normalize_t<S>>::value...>::value,
                  "all slices are normalized");
    static_assert(
        utils::all_of<(std::is_same<S, long>::value || is_normalized_slice<S>::value)...>::value,
        "all slices are valid");
    static_assert(std::decay_t<Arg>::value >= sizeof...(S), "slicing respects array shape");

    // numpy_gexpr is a wrapper for extended sliced array around a numpy
    // expression.
    // It contains compacted sorted slices value in lower, step && upper is
    // the same as shape.
    // indices for long index are store in the indices array.
    // position for slice and long value in the extended slice can be found
    // through the S... template
    // && compacted values as we know that first S is a slice.

    static_assert(utils::all_of<std::is_same<S, std::decay_t<S>>::value...>::value,
                  "no modifiers on slices");

    using dtype = typename std::remove_reference_t<Arg>::dtype;
    static constexpr size_t value = std::remove_reference_t<Arg>::value - count_long<S...>::value;

    using last_arg_stride_t = decltype(std::declval<Arg>().template strides<sizeof...(S) - 1>());
    using last_slice_t = std::tuple_element_t<sizeof...(S) - 1, std::tuple<S...>>;

    // It is not possible to vectorize everything. We only vectorize if the
    // last dimension is contiguous, which happens if
    // 1. Arg is an ndarray (this is too strict)
    // 2. the size of the gexpr is lower than the dim of arg, || it's the
    // same, but the last slice is contiguous
    static const bool is_vectorizable =
        std::remove_reference_t<Arg>::is_vectorizable &&
        (sizeof...(S) < std::remove_reference_t<Arg>::value ||
         std::is_same<cstride_normalized_slice<1>, last_slice_t>::value);
    static const bool is_flat =
        std::remove_reference_t<Arg>::is_flat && value == 1 &&
        utils::all_of<std::is_same<cstride_normalized_slice<1>, S>::value...>::value;
    static const bool is_strided =
        std::remove_reference_t<Arg>::is_strided ||
        (((sizeof...(S) - count_long<S...>::value) == value) &&
         !std::is_same<cstride_normalized_slice<1>, last_slice_t>::value);

    using value_type =
        std::decay_t<decltype(numpy_iexpr_helper<value>::get(std::declval<numpy_gexpr>(), 1))>;

    using iterator = std::conditional_t<is_strided || value != 1, nditerator<numpy_gexpr>, dtype *>;
    using const_iterator =
        std::conditional_t<is_strided || value != 1, const_nditerator<numpy_gexpr>, dtype const *>;

    std::remove_cv_t<Arg> arg;

    std::tuple<S...> slices;

    using shape_t = gexpr_shape_t<typename std::remove_reference_t<Arg>::shape_t, S...>;

    shape_t _shape;
    dtype *buffer;

    static constexpr types::pshape<std::integral_constant<long, 1>>
        last_stride(cstride_normalized_slice<1>);
    template <long stride>
    static constexpr types::pshape<std::integral_constant<long, stride>>
        last_stride(cstride_normalized_slice<stride>);
    static constexpr types::array_tuple<long, 1> last_stride(...);

    sutils::concat_t<types::array_tuple<long, value - 1>,
                     std::conditional_t<sizeof...(S) == std::decay_t<Arg>::value,
                                        decltype(last_stride(std::declval<last_slice_t>())),
                                        types::array_tuple<long, 1>>>
        _strides; // strides

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

    numpy_gexpr();
    numpy_gexpr(numpy_gexpr const &) = default;
    numpy_gexpr(numpy_gexpr &&) = default;

    template <class Argp> // ! using the default one, to make it possible to
    // accept reference && non reference version of
    // Argp
    numpy_gexpr(numpy_gexpr<Argp, S...> const &other);

    template <size_t J, class Slice>
    std::enable_if_t<is_normalized_slice<Slice>::value, void>
    init_shape(Slice const &s, utils::int_<1>, utils::int_<J>);

    template <size_t I, size_t J, class Slice>
    std::enable_if_t<is_normalized_slice<Slice>::value, void>
    init_shape(Slice const &s, utils::int_<I>, utils::int_<J>);

    template <size_t J>
    void init_shape(long cs, utils::int_<1>, utils::int_<J>);

    template <size_t I, size_t J>
    void init_shape(long cs, utils::int_<I>, utils::int_<J>);

    // private because we must use the make_gexpr factory to create a gexpr
  private:
    template <class _Arg, class... _other_classes>
    friend struct details::make_gexpr;

    friend struct array_base_slicer;
    template <class _Arg, class... _other_classes>
    friend std::enable_if_t<count_new_axis<_other_classes...>::value == 0,
                            numpy_gexpr<_Arg, _other_classes...>>
    details::_make_gexpr(_Arg arg, std::tuple<_other_classes...> const &t);
    template <class _Arg, class _other_classes, size_t... Is>
    friend numpy_gexpr<
        _Arg, typename to_normalized_slice<std::tuple_element_t<Is, _other_classes>>::type...>
    details::_make_gexpr_helper(_Arg arg, _other_classes const &s, std::index_sequence<Is...>);

    template <size_t C>
    friend struct extended_slice;

#ifdef ENABLE_PYTHON_MODULE
    template <typename T>
    friend struct pythonic::from_python;
#endif

    // When we create a new numpy_gexpr, we deduce step, lower && shape from
    // slices
    // && indices from long value.
    // Also, last shape information are set from origin array like in :
    // >>> a = numpy.arange(2*3*4).reshape(2,3,4)
    // >>> a[:, 1]
    // the last dimension (4) is missing from slice information
    // Finally, if origin expression was already sliced, lower bound && step
    // have to
    // be increased
    numpy_gexpr(Arg const &arg, std::tuple<S const &...> const &values);
    numpy_gexpr(Arg const &arg, S const &...s);

  public:
    template <class Argp, class... Sp>
    numpy_gexpr(numpy_gexpr<Argp, Sp...> const &expr, Arg arg);

    template <class G>
    numpy_gexpr(G const &expr, Arg &&arg);
    template <class pS>
    ndarray<dtype, pS> reshape(pS const &shape) const
    {
      return copy().reshape(shape);
    }

    template <class E>
    std::enable_if_t<may_overlap_gexpr<E>::value, numpy_gexpr &> _copy(E const &expr);

    template <class E>
    std::enable_if_t<!may_overlap_gexpr<E>::value, numpy_gexpr &> _copy(E const &expr);

    template <class E>
    numpy_gexpr &_copy_restrict(E const &expr);

    template <class E>
    numpy_gexpr &operator=(E const &expr);

    numpy_gexpr &operator=(numpy_gexpr const &expr);

    template <class Argp>
    numpy_gexpr &operator=(numpy_gexpr<Argp, S...> const &expr);

    template <class Op, class E>
    std::enable_if_t<may_overlap_gexpr<E>::value, numpy_gexpr &> update_(E const &expr);

    template <class Op, class E>
    std::enable_if_t<!may_overlap_gexpr<E>::value, numpy_gexpr &> update_(E const &expr);

    template <class E>
    numpy_gexpr &operator+=(E const &expr);

    numpy_gexpr &operator+=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator-=(E const &expr);

    numpy_gexpr &operator-=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator*=(E const &expr);

    numpy_gexpr &operator*=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator/=(E const &expr);

    numpy_gexpr &operator/=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator|=(E const &expr);

    numpy_gexpr &operator|=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator&=(E const &expr);

    numpy_gexpr &operator&=(numpy_gexpr const &expr);

    template <class E>
    numpy_gexpr &operator^=(E const &expr);

    numpy_gexpr &operator^=(numpy_gexpr const &expr);

    const_iterator begin() const;
    const_iterator end() const;

    iterator begin();
    iterator end();

    auto fast(long i) const & -> decltype(numpy_iexpr_helper<value>::get(*this, i))
    {
      return numpy_iexpr_helper<value>::get(*this, i);
    }

    auto fast(long i) & -> decltype(numpy_iexpr_helper<value>::get(*this, i))
    {
      return numpy_iexpr_helper<value>::get(*this, i);
    }

    template <class E, class... Indices>
    void store(E elt, Indices... indices)
    {
      static_assert(is_dtype<E>::value, "valid store");
      *(buffer + noffset<value>{}(*this, array_tuple<long, value>{{indices...}})) =
          static_cast<E>(elt);
    }
    template <class... Indices>
    dtype load(Indices... indices) const
    {
      return *(buffer + noffset<value>{}(*this, array_tuple<long, value>{{indices...}}));
    }
    template <class Op, class E, class... Indices>
    void update(E elt, Indices... indices) const
    {
      static_assert(is_dtype<E>::value, "valid store");
      Op{}(*(buffer + noffset<value>{}(*this, array_tuple<long, value>{{indices...}})),
           static_cast<E>(elt));
    }

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<numpy_gexpr>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    template <class... Sp>
    auto operator()(Sp const &...s) const -> decltype(make_gexpr(*this, s...));

    template <class Sp>
    auto operator[](Sp const &s) const
        -> std::enable_if_t<is_slice<Sp>::value, decltype(make_gexpr(*this, (s.lower, s)))>;

    template <size_t M>
    auto fast(array_tuple<long, M> const &indices) const & -> decltype(nget<M - 1>().fast(*this,
                                                                                          indices));

    template <size_t M>
    auto
    fast(array_tuple<long, M> const &indices) && -> decltype(nget<M - 1>().fast(std::move(*this),
                                                                                indices));

    template <size_t M>
    auto operator[](
        array_tuple<long, M> const &indices) const & -> decltype(nget<M - 1>()(*this, indices));

    template <size_t M>
    auto
    operator[](array_tuple<long, M> const &indices) && -> decltype(nget<M - 1>()(std::move(*this),
                                                                                 indices));

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value,
                     numpy_vexpr<numpy_gexpr, F>>
    operator[](F const &filter) const
    {
      return {*this, filter};
    }
    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value,
                     numpy_vexpr<numpy_gexpr, F>>
    fast(F const &filter) const
    {
      return {*this, filter};
    }

    template <class F>
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value,
                     numpy_vexpr<numpy_gexpr, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F>
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value,
                     numpy_vexpr<numpy_gexpr, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;
    auto operator[](long i) const -> decltype(this->fast(i));
    auto operator[](long i) -> decltype(this->fast(i));

    // template <class... Sp>
    // auto operator()(long i, Sp const &... s) const
    //    -> decltype((*this)[i](s...));

    explicit operator bool() const;

    dtype *data()
    {
      return buffer;
    }
    const dtype *data() const
    {
      return buffer;
    }

    long flat_size() const;
    long size() const;
    ndarray<dtype, shape_t> copy() const
    {
      return {*this};
    }

    intptr_t baseid() const
    {
      return arg.baseid();
    }

    template <class Tp, size_t... Is>
    auto recast(std::index_sequence<Is...>)
    {
      return make_gexpr(arg.template recast<Tp>(),
                        recast_slice<sizeof(dtype), sizeof(Tp)>(std::get<Is>(slices))...);
    }

    template <class Tp>
    auto recast()
    {
      return recast<Tp>(std::make_index_sequence<sizeof...(S)>());
    }
  };
} // namespace types

template <class E, class... S>
struct assignable_noescape<types::numpy_gexpr<E, S...>> {
  using type = types::numpy_gexpr<E, S...>;
};

template <class T, class pS, class... S>
struct assignable<types::numpy_gexpr<types::ndarray<T, pS> const &, S...>> {
  using type = types::numpy_gexpr<types::ndarray<T, pS>, S...>;
};

template <class T, class pS, class... S>
struct assignable<types::numpy_gexpr<types::ndarray<T, pS> &, S...>> {
  using type = types::numpy_gexpr<types::ndarray<T, pS>, S...>;
};

template <class Arg, class... S>
struct assignable<types::numpy_gexpr<Arg, S...>> {
  using type = types::numpy_gexpr<typename assignable<Arg>::type, S...>;
};

template <class Arg, class... S>
struct lazy<types::numpy_gexpr<Arg, S...>> : assignable<types::numpy_gexpr<Arg, S...>> {
};

PYTHONIC_NS_END
/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class Arg, class... S>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>,
                  pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::numpy_gexpr<Arg, S...>;
};

template <class Arg, class... S, class Argp, class... Sp>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>,
                  pythonic::types::numpy_gexpr<Argp, Sp...>> {
  using t0 = pythonic::types::numpy_gexpr<Arg, S...>;
  using t1 = pythonic::types::numpy_gexpr<Argp, Sp...>;
  using type =
      pythonic::types::ndarray < typename __combined<typename t0::dtype, typename t1::dtype>::type,
        pythonic::types::array_tuple<long, t0::value<t1::value ? t1::value : t0::value>>;
};

template <class Arg, class... S, class O>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>, O> {
  using type = pythonic::types::numpy_gexpr<Arg, S...>;
};
template <class Arg, class... S, class T>
struct __combined<pythonic::types::list<T>, pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::list<
      typename __combined<typename pythonic::types::numpy_gexpr<Arg, S...>::value_type, T>::type>;
};
template <class Arg, class... S, class T>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>, pythonic::types::list<T>> {
  using type = pythonic::types::list<
      typename __combined<typename pythonic::types::numpy_gexpr<Arg, S...>::value_type, T>::type>;
};
template <class Arg, class... S>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>, pythonic::types::none_type> {
  using type = pythonic::types::none<pythonic::types::numpy_gexpr<Arg, S...>>;
};

template <class Arg, class... S, class O>
struct __combined<O, pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::numpy_gexpr<Arg, S...>;
};

template <class Arg, class... S, class O>
struct __combined<O &, pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::numpy_gexpr<Arg, S...>;
};

template <class Arg, class... S>
struct __combined<pythonic::types::none_type, pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::none<pythonic::types::numpy_gexpr<Arg, S...>>;
};

/* combined are sorted such that the assigned type comes first */
template <class Arg, class... S, class T, class pS>
struct __combined<pythonic::types::numpy_gexpr<Arg, S...>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class Arg, class... S, class T, class pS>
struct __combined<pythonic::types::ndarray<T, pS>, pythonic::types::numpy_gexpr<Arg, S...>> {
  using type = pythonic::types::ndarray<T, pS>;
};

#endif
