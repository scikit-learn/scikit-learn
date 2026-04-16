#ifndef PYTHONIC_INCLUDE_NUMPY_MEAN_HPP
#define PYTHONIC_INCLUDE_NUMPY_MEAN_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/expand_dims.hpp"
#include "pythonic/include/numpy/sum.hpp"
#include "pythonic/include/types/immediate.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <size_t N>
    struct make_scalar_pshape : sutils::concat<types::pshape<std::integral_constant<long, 1>>,
                                               typename make_scalar_pshape<N - 1>::type> {
    };

    template <>
    struct make_scalar_pshape<1> {
      using type = types::pshape<std::integral_constant<long, 1>>;
    };

    template <class dtype>
    struct dtype_or_double_helper {
      using type = typename dtype::type;
    };
    template <>
    struct dtype_or_double_helper<types::none_type> {
      using type = double;
    };
    template <class dtype>
    using dtype_or_double = typename dtype_or_double_helper<dtype>::type;
  } // namespace details

  template <class E, class dtype = types::none_type>
  auto mean(E const &expr, types::none_type axis = {}, dtype d = {}, types::none_type out = {},
            types::false_immediate keep_dims = {})
      -> decltype(sum(expr, axis, d) / details::dtype_or_double<dtype>(expr.flat_size()));

  template <class E, class dtype = types::none_type>
  auto mean(E const &expr, long axis, dtype d = {}, types::none_type out = {},
            types::false_immediate keep_dims = {}) -> decltype(sum(expr, axis, d));

  template <class E, class dtype>
  types::ndarray<details::dtype_or_double<dtype>,
                 typename details::make_scalar_pshape<E::value>::type>
  mean(E const &expr, types::none_type axis, dtype d, types::none_type out,
       types::true_immediate keep_dims);

  template <class E, class dtype>
  auto mean(E const &expr, long axis, dtype d, types::none_type out,
            types::true_immediate keep_dims) -> decltype(expand_dims(mean(expr, axis, d), axis));

  DEFINE_FUNCTOR(pythonic::numpy, mean);
} // namespace numpy
PYTHONIC_NS_END

#endif
