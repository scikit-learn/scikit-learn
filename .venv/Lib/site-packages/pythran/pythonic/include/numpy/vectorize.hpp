#ifndef PYTHONIC_INCLUDE_NUMPY_VECTORIZE_HPP
#define PYTHONIC_INCLUDE_NUMPY_VECTORIZE_HPP

#include "pythonic/include/types/numpy_expr.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F>
  struct vectorized {
    using callable = void;
    template <typename... T>
    auto operator()(T &&...args) const
        -> std::enable_if_t<!types::valid_numexpr_parameters<std::decay_t<T>...>::value,
                            decltype(F{}(std::forward<T>(args)...))>;

    template <class... E>
    std::enable_if_t<types::valid_numexpr_parameters<std::decay_t<E>...>::value,
                     types::numpy_expr<F, typename types::adapt_type<E, E...>::type...>>
    operator()(E &&...args) const;
  };

  template <class F>
  vectorized<F> vectorize(F const &);

  DEFINE_FUNCTOR(pythonic::numpy, vectorize);
} // namespace numpy
PYTHONIC_NS_END

#endif
