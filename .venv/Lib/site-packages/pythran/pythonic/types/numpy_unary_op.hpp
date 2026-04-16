#ifndef NUMPY_UNARY_FUNC_NAME
#error NUMPY_UNARY_FUNC_NAME undefined
#endif
#ifndef NUMPY_UNARY_FUNC_SYM
#error NUMPY_UNARY_FUNC_SYM undefined
#endif

template <class E>
std::enable_if_t<types::valid_numop_parameters<std::decay_t<E>>::value,
                 types::numpy_expr<NUMPY_UNARY_FUNC_SYM, E>>
NUMPY_UNARY_FUNC_NAME(E &&self)
{
  return {std::forward<E>(self)};
}

#undef NUMPY_UNARY_FUNC_NAME
#undef NUMPY_UNARY_FUNC_SYM
