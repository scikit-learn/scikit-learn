#ifndef NUMPY_NARY_FUNC_NAME
#error NUMPY_NARY_FUNC_NAME undefined
#endif

#ifndef NUMPY_NARY_FUNC_SYM
#error NUMPY_NARY_FUNC_SYM undefined
#endif

#ifndef NUMPY_NARY_RESHAPE_MODE
#define NUMPY_NARY_RESHAPE_MODE adapt_type
#endif

#ifndef NUMPY_NARY_EXTRA_METHOD
#define NUMPY_NARY_EXTRA_METHOD
#endif

#define STR(a) STR_(a)

namespace functor
{

  struct NUMPY_NARY_FUNC_NAME {
    using callable = void;

    // We accept implementation here
    NUMPY_NARY_EXTRA_METHOD

    template <typename... T>
    auto operator()(T &&...args) const
        -> std::enable_if_t<!types::valid_numexpr_parameters<std::decay_t<T>...>::value,
                            decltype(NUMPY_NARY_FUNC_SYM(std::forward<T>(args)...))>;

    template <class... E>
    std::enable_if_t<types::valid_numexpr_parameters<std::decay_t<E>...>::value,
                     types::numpy_expr<NUMPY_NARY_FUNC_NAME,
                                       typename types::NUMPY_NARY_RESHAPE_MODE<E, E...>::type...>>
    operator()(E &&...args) const;

    friend std::ostream &operator<<(std::ostream &os, NUMPY_NARY_FUNC_NAME)
    {
      return os << STR(NUMPY_NARY_FUNC_NAME);
    }
  };
} // namespace functor

#undef NUMPY_NARY_FUNC_NAME
#undef NUMPY_NARY_FUNC_SYM
#undef NUMPY_NARY_RESHAPE_MODE
#undef NUMPY_NARY_EXTRA_METHOD
#undef STR
