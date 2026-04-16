#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_EXPR_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_EXPR_HPP

#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/utils/meta.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  template <size_t I, class Args>
  bool is_trivial_broadcast()
  {
    return std::is_same<
        std::tuple_element_t<0, typename std::decay_t<std::tuple_element_t<I, Args>>::shape_t>,
        std::integral_constant<long, 1>>::value;
  }

  template <class... Tys>
  struct count_non_integral;
  template <>
  struct count_non_integral<long> : std::integral_constant<long, 1> {
  };
  template <class T>
  struct count_non_integral<T> : std::integral_constant<long, 0> {
  };
  template <class Ty0, class Ty1, class... Tys>
  struct count_non_integral<Ty0, Ty1, Tys...>
      : std::integral_constant<long, count_non_integral<Ty0>::value +
                                         count_non_integral<Ty1, Tys...>::value> {
  };
  template <class P>
  struct is_perfect_stepping;
  template <class... Tys>
  struct is_perfect_stepping<pshape<Tys...>>
      : std::integral_constant<bool, count_non_integral<Tys...>::value == 1> {
  };
  template <size_t value, class Args, size_t N, size_t... Is>
  struct all_valid_indices;

  template <size_t value, class Args, size_t... Is>
  struct all_valid_indices<value, Args, 0, Is...> {
    using type = std::index_sequence<Is...>;
  };
  template <size_t value, class Args, size_t N, size_t... Is>
  struct all_valid_indices
      : std::conditional_t<(value <=
                            std::remove_reference_t<std::tuple_element_t<N - 1, Args>>::value),
                           all_valid_indices<value, Args, N - 1, Is..., N - 1>,
                           all_valid_indices<value, Args, N - 1, Is...>> {
  };

  template <size_t value, class Args>
  using valid_indices = typename all_valid_indices<value, Args, std::tuple_size<Args>::value>::type;

  template <class Expr>
  struct is_numexpr_arg;

  template <class Expr, class... Slice>
  struct numpy_gexpr;

  template <class Op>
  struct Dereferencer {
    template <class Ts, size_t... I>
    auto operator()(Ts const &iters, std::index_sequence<I...>)
        -> decltype(Op{}(*std::get<I>(iters)...))
    {
      return Op{}(*std::get<I>(iters)...);
    }
  };
  template <class T>
  struct step {
    using type = typename T::step_type;
  };
  namespace details
  {
    template <size_t I, class Args, size_t... Is>
    long init_shape_element(Args const &args, std::index_sequence<Is...>);
  }

  template <class Op, class Steps, class... Iters>
  struct numpy_expr_iterator
      : std::iterator<
            std::random_access_iterator_tag,
            std::remove_reference_t<decltype(std::declval<Op>()(*std::declval<Iters>()...))>> {
    Steps steps_;
    std::tuple<Iters...> iters_;

    numpy_expr_iterator(Steps steps, Iters... iters) : steps_(steps), iters_(iters...)
    {
    }

    numpy_expr_iterator(numpy_expr_iterator const &other)
        : steps_(other.steps_), iters_(other.iters_)
    {
    }

    numpy_expr_iterator &operator=(numpy_expr_iterator const &other)
    {
      iters_ = other.iters_;
      return *this;
    }

    template <size_t... I>
    auto _dereference(std::index_sequence<I...> s) const -> decltype(Dereferencer<Op>{}(iters_, s))
    {
      return Dereferencer<Op>{}(iters_, s);
    }

    auto operator*() const
        -> decltype(this->_dereference(std::make_index_sequence<sizeof...(Iters)>{}))
    {
      return _dereference(std::make_index_sequence<sizeof...(Iters)>{});
    }

    template <size_t I>
    bool _incr_opt(std::integral_constant<bool, true> long_step)
    {
      if (is_perfect_stepping<Steps>::value)
        ++std::get<I>(iters_);
      else
        std::get<I>(iters_) += std::get<I>(steps_);
      return true;
    }

    template <size_t I>
    bool _incr_opt(std::integral_constant<bool, false> long_step)
    {
      if (std::tuple_element_t<I, Steps>::value)
        ++std::get<I>(iters_);
      return true;
    }

    template <size_t... I>
    void _incr(std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{_incr_opt<I>(
          std::integral_constant<bool,
                                 std::is_same<long, std::tuple_element_t<I, Steps>>::value>{})...};
    }
    numpy_expr_iterator &operator++()
    {
      _incr(std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    numpy_expr_iterator operator+(long i) const
    {
      numpy_expr_iterator other(*this);
      return other += i;
    }

    template <size_t... I>
    void _update(long i, std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{(std::get<I>(iters_) += std::get<I>(steps_) * i, true)...};
    }
    numpy_expr_iterator &operator+=(long i)
    {
      _update(i, std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    template <size_t... I>
    long _difference(numpy_expr_iterator const &other, std::index_sequence<I...>) const
    {
      std::initializer_list<long> distances{
          (static_cast<long>(std::get<I>(iters_) - std::get<I>(other.iters_)))...};
      return *std::max_element(distances.begin(), distances.end());
    }

    long operator-(numpy_expr_iterator const &other) const
    {
      return _difference(other, std::make_index_sequence<sizeof...(Iters)>{});
    }

    bool _neq(numpy_expr_iterator const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _neq(numpy_expr_iterator const &other, utils::int_<I>) const
    {
      return (std::get<I - 1>(steps_) &&
              (std::get<I - 1>(iters_) != std::get<I - 1>(other.iters_))) ||
             _neq(other, utils::int_<I - 1>{});
    }

    bool operator!=(numpy_expr_iterator const &other) const
    {
      return _neq(other, utils::int_<sizeof...(Iters)>{});
    }
    bool _eq(numpy_expr_iterator const &other, utils::int_<0u>) const
    {
      return true;
    }

    template <size_t I>
    bool _eq(numpy_expr_iterator const &other, utils::int_<I>) const
    {
      return (!std::get<I - 1>(steps_) ||
              (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_))) &&
             _eq(other, utils::int_<I - 1>{});
    }

    bool operator==(numpy_expr_iterator const &other) const
    {
      return _eq(other, utils::int_<sizeof...(Iters)>{});
    }

    bool _lt(numpy_expr_iterator const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _lt(numpy_expr_iterator const &other, utils::int_<I>) const
    {
      if (!std::get<I - 1>(steps_) || (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_)))
        return _lt(other, utils::int_<I - 1>{});
      else
        return std::get<I - 1>(steps_) && (std::get<I - 1>(iters_) < std::get<I - 1>(other.iters_));
    }

    bool operator<(numpy_expr_iterator const &other) const
    {
      return _lt(other, utils::int_<sizeof...(Iters)>{});
    }

    bool operator<=(numpy_expr_iterator const &other) const
    {
      return *this < other || *this == other;
    }
  };
#ifdef USE_XSIMD
  template <class E, class Op, class Steps, class SIters, class... Iters>
  struct numpy_expr_simd_iterator
      : std::iterator<
            std::random_access_iterator_tag,
            std::remove_reference_t<decltype(std::declval<Op>()(*std::declval<Iters>()...))>> {
    Steps steps_;
    std::tuple<Iters...> iters_;
    SIters siters_;

    numpy_expr_simd_iterator(array_tuple<long, sizeof...(Iters)> steps, SIters const &siters,
                             Iters... iters)
        : steps_(steps), iters_(iters...), siters_(siters)
    {
    }

    numpy_expr_simd_iterator(numpy_expr_simd_iterator const &other)
        : steps_(other.steps_), iters_(other.iters_), siters_(other.siters_)
    {
    }

    numpy_expr_simd_iterator &operator=(numpy_expr_simd_iterator const &other)
    {
      iters_ = other.iters_;
      siters_ = other.siters_;
      return *this;
    }

    template <size_t... I>
    auto _dereference(std::index_sequence<I...>) const -> decltype(Op{}(*std::get<I>(iters_)...))
    {
      return Op{}(((std::get<I>(steps_)) ? (*std::get<I>(iters_)) : (std::get<I>(siters_)))...);
    }

    auto operator*() const
        -> decltype(this->_dereference(std::make_index_sequence<sizeof...(Iters)>{}))
    {
      return _dereference(std::make_index_sequence<sizeof...(Iters)>{});
    }

    template <size_t I>
    bool _incr_opt(std::integral_constant<bool, true> long_step)
    {
      if (is_perfect_stepping<Steps>::value)
        ++std::get<I>(iters_);
      else
        std::get<I>(iters_) += std::get<I>(steps_);
      return true;
    }

    template <size_t I>
    bool _incr_opt(std::integral_constant<bool, false> long_step)
    {
      if (std::tuple_element_t<I, Steps>::value)
        ++std::get<I>(iters_);
      return true;
    }

    template <size_t... I>
    void _incr(std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{_incr_opt<I>(
          std::integral_constant<bool,
                                 std::is_same<long, std::tuple_element_t<I, Steps>>::value>{})...};
    }
    numpy_expr_simd_iterator &operator++()
    {
      _incr(std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    numpy_expr_simd_iterator operator+(long i) const
    {
      numpy_expr_simd_iterator other(*this);
      return other += i;
    }

    template <size_t... I>
    void _update(long i, std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{(std::get<I>(iters_) += std::get<I>(steps_) * i, true)...};
    }
    numpy_expr_simd_iterator &operator+=(long i)
    {
      _update(i, std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    template <size_t... I>
    long _difference(numpy_expr_simd_iterator const &other, std::index_sequence<I...>) const
    {
      std::initializer_list<long> distances{(std::get<I>(iters_) - std::get<I>(other.iters_))...};
      return *std::max_element(distances.begin(), distances.end());
    }

    long operator-(numpy_expr_simd_iterator const &other) const
    {
      return _difference(other, std::make_index_sequence<sizeof...(Iters)>{});
    }

    bool _neq(numpy_expr_simd_iterator const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _neq(numpy_expr_simd_iterator const &other, utils::int_<I>) const
    {
      return (std::get<I - 1>(steps_) &&
              (std::get<I - 1>(iters_) != std::get<I - 1>(other.iters_))) ||
             _neq(other, utils::int_<I - 1>{});
    }

    bool operator!=(numpy_expr_simd_iterator const &other) const
    {
      return _neq(other, utils::int_<sizeof...(Iters)>{});
    }

    bool _eq(numpy_expr_simd_iterator const &other, utils::int_<0u>) const
    {
      return true;
    }

    template <size_t I>
    bool _eq(numpy_expr_simd_iterator const &other, utils::int_<I>) const
    {
      return (std::get<I - 1>(steps_) &&
              (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_))) &&
             _eq(other, utils::int_<I - 1>{});
    }

    bool operator==(numpy_expr_simd_iterator const &other) const
    {
      return _eq(other, utils::int_<sizeof...(Iters)>{});
    }

    bool _lt(numpy_expr_simd_iterator const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _lt(numpy_expr_simd_iterator const &other, utils::int_<I>) const
    {
      if (std::get<I - 1>(steps_) && (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_)))
        return _lt(other, utils::int_<I - 1>{});
      else
        return std::get<I - 1>(steps_) && (std::get<I - 1>(iters_) < std::get<I - 1>(other.iters_));
    }

    bool operator<(numpy_expr_simd_iterator const &other) const
    {
      return _lt(other, utils::int_<sizeof...(Iters)>{});
    }
  };

  template <class E, class Op, class... Iters>
  struct numpy_expr_simd_iterator_nobroadcast
      : std::iterator<
            std::random_access_iterator_tag,
            std::remove_reference_t<decltype(std::declval<Op>()(*std::declval<Iters>()...))>> {
    std::tuple<Iters...> iters_;

    numpy_expr_simd_iterator_nobroadcast(Iters... iters) : iters_(iters...)
    {
    }

    numpy_expr_simd_iterator_nobroadcast(numpy_expr_simd_iterator_nobroadcast const &other)
        : iters_(other.iters_)
    {
    }

    numpy_expr_simd_iterator_nobroadcast &
    operator=(numpy_expr_simd_iterator_nobroadcast const &other)
    {
      iters_ = other.iters_;
      return *this;
    }

    template <size_t... I>
    auto _dereference(std::index_sequence<I...>) const -> decltype(Op{}(*std::get<I>(iters_)...))
    {
      return Op{}((*std::get<I>(iters_))...);
    }

    auto operator*() const
        -> decltype(this->_dereference(std::make_index_sequence<sizeof...(Iters)>{}))
    {
      return _dereference(std::make_index_sequence<sizeof...(Iters)>{});
    }

    template <size_t... I>
    void _incr(std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{(++std::get<I>(iters_), true)...};
    }
    numpy_expr_simd_iterator_nobroadcast &operator++()
    {
      _incr(std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    template <size_t... I>
    long _difference(numpy_expr_simd_iterator_nobroadcast const &other,
                     std::index_sequence<I...>) const
    {
      std::initializer_list<long> distances{(std::get<I>(iters_) - std::get<I>(other.iters_))...};
      return *std::max_element(distances.begin(), distances.end());
    }

    long operator-(numpy_expr_simd_iterator_nobroadcast const &other) const
    {
      return _difference(other, std::make_index_sequence<sizeof...(Iters)>{});
    }

    numpy_expr_simd_iterator_nobroadcast operator+(long i) const
    {
      numpy_expr_simd_iterator_nobroadcast other(*this);
      return other += i;
    }

    template <size_t... I>
    void _update(long i, std::index_sequence<I...>)
    {
      (void)std::initializer_list<bool>{(std::get<I>(iters_) += i, true)...};
    }
    numpy_expr_simd_iterator_nobroadcast &operator+=(long i)
    {
      _update(i, std::make_index_sequence<sizeof...(Iters)>{});
      return *this;
    }

    bool _neq(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _neq(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<I>) const
    {
      return (std::get<I - 1>(iters_) != std::get<I - 1>(other.iters_)) ||
             _neq(other, utils::int_<I - 1>{});
    }

    bool operator!=(numpy_expr_simd_iterator_nobroadcast const &other) const
    {
      return _neq(other, utils::int_<sizeof...(Iters)>{});
    }

    bool _eq(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<0u>) const
    {
      return true;
    }

    template <size_t I>
    bool _eq(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<I>) const
    {
      return (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_)) &&
             _eq(other, utils::int_<I - 1>{});
    }

    bool operator==(numpy_expr_simd_iterator_nobroadcast const &other) const
    {
      return _eq(other, utils::int_<sizeof...(Iters)>{});
    }

    bool _lt(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<0u>) const
    {
      return false;
    }

    template <size_t I>
    bool _lt(numpy_expr_simd_iterator_nobroadcast const &other, utils::int_<I>) const
    {
      if (std::get<I - 1>(iters_) == std::get<I - 1>(other.iters_))
        return _lt(other, utils::int_<I - 1>{});
      else
        return std::get<I - 1>(iters_) < std::get<I - 1>(other.iters_);
    }

    bool operator<(numpy_expr_simd_iterator_nobroadcast const &other) const
    {
      return _lt(other, utils::int_<sizeof...(Iters)>{});
    }
  };
#endif

  template <long N0, long N1>
  std::integral_constant<long, N0 == N1> make_step(std::integral_constant<long, N0>,
                                                   std::integral_constant<long, N1>)
  {
    return {};
  }
  template <class T0, class T1>
  long make_step(T0 n0, T1 n1)
  {
    return (long)n0 == (long)n1;
  }

  template <class S>
  constexpr size_t count_none(size_t I)
  {
    return I == 0 ? 0 : std::is_same<S, none_type>::value;
  }

  template <class S, class Sp, class... Ss>
  constexpr size_t count_none(size_t I)
  {
    return I == 0 ? 0 : (std::is_same<S, none_type>::value + count_none<Sp, Ss...>(I - 1));
  }

  template <class BT, class T>
  using step_type_t =
      decltype(make_step(std::get<0>(std::declval<BT>()), std::get<0>(std::declval<T>())));

  constexpr size_t clamp(size_t i, size_t j)
  {
    return i > j ? j : i;
  }

  template <size_t... J, class Arg, class Shp, class... S>
  auto make_subslice(std::index_sequence<J...>, Arg const &arg, Shp const &shp,
                     std::tuple<S...> const &ss) -> decltype(arg(std::get<J>(ss)...))
  {
    // we need to adapt_slice to take broadcasting into account
    return arg(adapt_slice(
        std::get<J>(ss), shp.template shape<clamp(J - count_none<S...>(J), Shp::value - 1)>(),
        arg.template shape<clamp(J - count_none<S...>(J), Arg::value - 1)>())...);
  }

  /* Expression template for numpy expressions - binary operators
   */
  template <class Op, class... Args>
  struct numpy_expr {
    using first_arg = typename utils::front<Args...>::type;
    static const bool is_vectorizable =
        utils::all_of<std::remove_reference_t<Args>::is_vectorizable...>::value &&
        utils::all_of<std::is_same<
            typename std::remove_cv_t<std::remove_reference_t<first_arg>>::dtype,
            typename std::remove_cv_t<std::remove_reference_t<Args>>::dtype>::value...>::value &&
        types::is_vector_op<Op, typename std::remove_reference_t<Args>::dtype...>::value;
    static const bool is_flat = false;
    static const bool is_strided =
        utils::any_of<std::remove_reference_t<Args>::is_strided...>::value;

    static constexpr size_t value =
        utils::max_element<std::remove_reference_t<Args>::value...>::value;
    using value_type =
        decltype(Op()(std::declval<typename std::remove_reference_t<Args>::value_type>()...));
    using dtype = decltype(Op()(std::declval<typename std::remove_reference_t<Args>::dtype>()...));

#ifdef CYTHON_ABI
    std::tuple<std::remove_reference_t<Args>...> args;
#else
    std::tuple<Args...> args;
#endif
    using shape_t =
        sutils::merged_shapes_t<value, typename std::remove_reference_t<Args>::shape_t...>;
    using steps_t =
        pshape<step_type_t<shape_t, typename std::remove_reference_t<Args>::shape_t>...>;
    static_assert(value == std::tuple_size<shape_t>::value, "consistent shape and size");
    using const_iterator =
        numpy_expr_iterator<Op, steps_t, typename std::remove_reference_t<Args>::const_iterator...>;
    using iterator =
        numpy_expr_iterator<Op, steps_t, typename std::remove_reference_t<Args>::iterator...>;
    using const_fast_iterator = const_nditerator<numpy_expr>;

    numpy_expr() = default;
    numpy_expr(numpy_expr const &) = default;
    numpy_expr(numpy_expr &&) = default;

    template <class... Argp>
    numpy_expr(numpy_expr<Op, Argp...> const &other) : args(other.args)
    {
    }

    numpy_expr(Args const &...args);

    template <size_t... I>
    const_iterator _begin(std::index_sequence<I...>) const;
    const_iterator begin() const;

    template <size_t... I>
    const_iterator _end(std::index_sequence<I...>) const;
    const_iterator end() const;

    const_fast_iterator begin(types::fast) const;
    const_fast_iterator end(types::fast) const;

    template <size_t... I>
    iterator _begin(std::index_sequence<I...>);
    iterator begin();

    template <size_t... I>
    iterator _end(std::index_sequence<I...>);
    iterator end();

    template <size_t... I>
    auto _fast(long i, std::index_sequence<I...>) const
        -> decltype(Op()(std::get<I>(args).fast(i)...))
    {
      return Op()(std::get<I>(args).fast(i)...);
    }

    auto fast(long i) const
        -> decltype(this->_fast(i, std::make_index_sequence<sizeof...(Args)>{}));

    template <class... Indices, size_t... I>
    auto _load(std::index_sequence<I...>, Indices... indices) const
        -> decltype(Op()(std::get<I>(args).load(indices...)...))
    {
      return Op()(std::get<I>(args).load(indices...)...);
    }

    template <class... Indices>
    auto load(Indices... indices) const
        -> decltype(this->_load(std::make_index_sequence<sizeof...(Args)>{}, indices...))
    {
      return this->_load(std::make_index_sequence<sizeof...(Args)>{}, indices...);
    }

    template <size_t... I>
    auto _map_fast(array_tuple<long, sizeof...(I)> const &indices, std::index_sequence<I...>) const
        -> decltype(Op()(std::get<I>(args).fast(std::get<I>(indices))...))
    {
      return Op()(std::get<I>(args).fast(std::get<I>(indices))...);
    }

    template <class... Indices>
    auto map_fast(Indices... indices) const
        -> decltype(this->_map_fast(array_tuple<long, sizeof...(Indices)>{{indices...}},
                                    std::make_index_sequence<sizeof...(Args)>{}));

  public:
    template <size_t I>
    auto shape() const
        -> decltype(details::init_shape_element<I>(args,
                                                   valid_indices<value, std::tuple<Args...>>{}))
    {
      return details::init_shape_element<I>(args, valid_indices<value, std::tuple<Args...>>{});
    }
    template <size_t... I>
    bool _no_broadcast(std::index_sequence<I...>) const;
    bool no_broadcast() const;
    template <size_t... I>
    bool _no_broadcast_vectorize(std::index_sequence<I...>) const;
    bool no_broadcast_vectorize() const;
    template <size_t... I>
    bool _no_broadcast_ex(std::index_sequence<I...>) const;
    bool no_broadcast_ex() const;

#ifdef USE_XSIMD
    using simd_iterator = numpy_expr_simd_iterator<
        numpy_expr, Op,
        pshape<step_type_t<shape_t, typename std::remove_reference_t<Args>::shape_t>...>,
        std::tuple<xsimd::batch<typename std::remove_reference_t<Args>::value_type>...>,
        typename std::remove_reference_t<Args>::simd_iterator...>;
    using simd_iterator_nobroadcast = numpy_expr_simd_iterator_nobroadcast<
        numpy_expr, Op, typename std::remove_reference_t<Args>::simd_iterator_nobroadcast...>;
    template <size_t... I>
    simd_iterator _vbegin(types::vectorize, std::index_sequence<I...>) const;
    simd_iterator vbegin(types::vectorize) const;
    template <size_t... I>
    simd_iterator _vend(types::vectorize, std::index_sequence<I...>) const;
    simd_iterator vend(types::vectorize) const;

    template <size_t... I>
    simd_iterator_nobroadcast _vbegin(types::vectorize_nobroadcast,
                                      std::index_sequence<I...>) const;
    simd_iterator_nobroadcast vbegin(types::vectorize_nobroadcast) const;
    template <size_t... I>
    simd_iterator_nobroadcast _vend(types::vectorize_nobroadcast, std::index_sequence<I...>) const;
    simd_iterator_nobroadcast vend(types::vectorize_nobroadcast) const;

#endif

    template <size_t... I, class... S>
    auto _get(std::index_sequence<I...> is, S const &...s) const
        -> decltype(Op{}(make_subslice(std::make_index_sequence<sizeof...(S)>{}, std::get<I>(args),
                                       *this, std::make_tuple(s...))...))
    {
      return Op{}(make_subslice(std::make_index_sequence<sizeof...(S)>{}, std::get<I>(args), *this,
                                std::make_tuple(s...))...);
    }

    template <class... S>
    auto operator()(S const &...s) const
        -> decltype(this->_get(std::make_index_sequence<sizeof...(Args)>{}, s...));

    template <class F>
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         !is_pod_array<F>::value,
                     numpy_vexpr<numpy_expr, ndarray<long, pshape<long>>>>
    fast(F const &filter) const;

    template <class F>
    std::enable_if_t<is_numexpr_arg<F>::value && std::is_same<bool, typename F::dtype>::value &&
                         !is_pod_array<F>::value,
                     numpy_vexpr<numpy_expr, ndarray<long, pshape<long>>>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_expr, F>>
    operator[](F const &filter) const;

    template <class F> // indexing through an array of indices -- a view
    std::enable_if_t<is_numexpr_arg<F>::value && !is_array_index<F>::value &&
                         !std::is_same<bool, typename F::dtype>::value && !is_pod_array<F>::value,
                     numpy_vexpr<numpy_expr, F>>
    fast(F const &filter) const;

    // FIXME: this does not take into account bounds and broadcasting
    auto operator[](long i) const -> decltype(this->fast(i));

    template <size_t... I, class S>
    auto _index(S s, std::index_sequence<I...>) const -> decltype(Op{}(std::get<I>(args)[s]...))
    {
      return Op{}(std::get<I>(args)[s]...);
    }
    template <class S>
    auto operator[](S s) const
        -> decltype((*this)._index((s.lower, s), std::make_index_sequence<sizeof...(Args)>{}))
    {
      return _index(s, std::make_index_sequence<sizeof...(Args)>{});
    }

    dtype operator[](array_tuple<long, value> const &indices) const
    {
      return _index(indices, std::make_index_sequence<sizeof...(Args)>{});
    }

    explicit operator bool() const;

    long flat_size() const;

    long size() const;
  };
} // namespace types

template <class Op, class... Args>
struct assignable<types::numpy_expr<Op, Args...>> {
  using type = types::ndarray<typename pythonic::types::numpy_expr<Op, Args...>::dtype,
                              typename pythonic::types::numpy_expr<Op, Args...>::shape_t>;
};

template <class Op, class... Arg>
struct lazy<types::numpy_expr<Op, Arg...>> {
  using type = types::numpy_expr<Op, typename lazy<Arg>::type...>;
};
PYTHONIC_NS_END
/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"
template <class Op, class K, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, indexable<K>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class K, class... Args>
struct __combined<indexable<K>, pythonic::types::numpy_expr<Op, Args...>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class K, class V, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, indexable_container<K, V>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class K, class V, class... Args>
struct __combined<indexable_container<K, V>, pythonic::types::numpy_expr<Op, Args...>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class K, class... Args>
struct __combined<container<K>, pythonic::types::numpy_expr<Op, Args...>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class K, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, container<K>> {
  using type = pythonic::types::numpy_expr<Op, Args...>;
};

template <class Op, class Op2, class... Args, class... Args2>
struct __combined<pythonic::types::numpy_expr<Op, Args...>,
                  pythonic::types::numpy_expr<Op2, Args2...>> {
  using type = pythonic::types::ndarray<
      typename pythonic::types::numpy_expr<Op, Args...>::dtype,
      pythonic::types::array_tuple<long, pythonic::types::numpy_expr<Op, Args...>::value>>;
};
template <class E, class Op, class... Args>
struct __combined<pythonic::types::numpy_iexpr<E>, pythonic::types::numpy_expr<Op, Args...>> {
  using type = pythonic::types::numpy_iexpr<E>;
};
template <class E, class Op, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, pythonic::types::numpy_iexpr<E>> {
  using type = pythonic::types::numpy_iexpr<E>;
};

template <class T, class pS, class Op, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, pythonic::types::ndarray<T, pS>> {
  using type = pythonic::types::ndarray<T, pS>;
};

template <class T, class Op, class... Args>
struct __combined<pythonic::types::numpy_expr<Op, Args...>, pythonic::types::numpy_texpr<T>> {
  using type = pythonic::types::ndarray<
      typename pythonic::types::numpy_expr<Op, Args...>::dtype,
      pythonic::types::array_tuple<long, pythonic::types::numpy_expr<Op, Args...>::value>>;
};

template <class T, class Op, class... Args>
struct __combined<pythonic::types::numpy_texpr<T>, pythonic::types::numpy_expr<Op, Args...>> {
  using type = pythonic::types::ndarray<
      typename pythonic::types::numpy_expr<Op, Args...>::dtype,
      pythonic::types::array_tuple<long, pythonic::types::numpy_expr<Op, Args...>::value>>;
};

/*}*/
#endif
