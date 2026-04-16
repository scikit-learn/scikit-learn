#ifndef NUMPY_BINARY_FUNC_NAME
#error NUMPY_BINARY_FUNC_NAME undefined
#endif
#ifndef NUMPY_BINARY_FUNC_SYM
#error NUMPY_BINARY_FUNC_SYM undefined
#endif

template <class E0, class E1>
std::enable_if_t<types::valid_numop_parameters<std::decay_t<E0>, std::decay_t<E1>>::value,
                 types::numpy_expr<NUMPY_BINARY_FUNC_SYM, typename types::adapt_type<E0, E1>::type,
                                   typename types::adapt_type<E1, E0>::type>>
NUMPY_BINARY_FUNC_NAME(E0 &&self, E1 &&other)
{
  return {std::forward<E0>(self), std::forward<E1>(other)};
}

#undef NUMPY_BINARY_FUNC_NAME
#undef NUMPY_BINARY_FUNC_SYM
