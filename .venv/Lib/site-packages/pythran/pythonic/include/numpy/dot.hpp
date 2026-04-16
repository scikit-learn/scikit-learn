#ifndef PYTHONIC_INCLUDE_NUMPY_DOT_HPP
#define PYTHONIC_INCLUDE_NUMPY_DOT_HPP

#include "pythonic/include/numpy/sum.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_expr.hpp"
#include "pythonic/include/types/traits.hpp"

template <class T>
struct is_blas_type : pythonic::types::is_complex<T> {
};

template <>
struct is_blas_type<float> : std::true_type {
};

template <>
struct is_blas_type<double> : std::true_type {
};

template <class E>
struct is_strided {
  template <class T>
  static decltype(T::is_strided, std::true_type{}) get(T *);
  static std::false_type get(...);
  static constexpr bool value = decltype(get((E *)nullptr))::value;
};

template <class E>
struct is_blas_array {
  static constexpr bool value = pythonic::types::has_buffer<E>::value &&
                                is_blas_type<typename pythonic::types::dtype_of<E>::type>::value &&
                                !is_strided<E>::value;
};

template <class E>
struct is_blas_view {
  static constexpr bool value = pythonic::types::has_buffer<E>::value &&
                                is_blas_type<typename pythonic::types::dtype_of<E>::type>::value;
};

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  std::enable_if_t<types::is_dtype<E>::value && types::is_dtype<F>::value,
                   decltype(std::declval<E>() * std::declval<F>())>
  dot(E const &e, F const &f);

  /// Vector / Vector multiplication
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value && types::is_numexpr_arg<F>::value &&
                       E::value == 1 && F::value == 1 &&
                       (!is_blas_view<E>::value || !is_blas_view<F>::value ||
                        !std::is_same<typename E::dtype, typename F::dtype>::value),
                   typename __combined<typename E::dtype, typename F::dtype>::type>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, float>::value &&
                       std::is_same<typename F::dtype, float>::value && is_blas_array<E>::value &&
                       is_blas_array<F>::value,
                   float>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, double>::value &&
                       std::is_same<typename F::dtype, double>::value && is_blas_array<E>::value &&
                       is_blas_array<F>::value,
                   double>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<float>>::value &&
                       std::is_same<typename F::dtype, std::complex<float>>::value &&
                       is_blas_array<E>::value && is_blas_array<F>::value,
                   std::complex<float>>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<double>>::value &&
                       std::is_same<typename F::dtype, std::complex<double>>::value &&
                       is_blas_array<E>::value && is_blas_array<F>::value,
                   std::complex<double>>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, float>::value &&
                       std::is_same<typename F::dtype, float>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   float>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, double>::value &&
                       std::is_same<typename F::dtype, double>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   double>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<float>>::value &&
                       std::is_same<typename F::dtype, std::complex<float>>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   std::complex<float>>
  dot(E const &e, F const &f);

  template <class E, class F>
  std::enable_if_t<E::value == 1 && F::value == 1 &&
                       std::is_same<typename E::dtype, std::complex<double>>::value &&
                       std::is_same<typename F::dtype, std::complex<double>>::value &&
                       (is_blas_view<E>::value && is_blas_view<F>::value &&
                        !(is_blas_array<E>::value && is_blas_array<F>::value)),
                   std::complex<double>>
  dot(E const &e, F const &f);

  /// Matrix / Vector multiplication

  // We transpose the matrix to reflect our C order
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 1,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &f, types::ndarray<E, pS1> const &e);

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 1,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &f, types::ndarray<E, pS1> const &e);

  // The trick is to not transpose the matrix so that MV become VM
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 1 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &e, types::ndarray<E, pS1> const &f);

  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 1 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::pshape<long>>>
  dot(types::ndarray<E, pS0> const &e, types::numpy_texpr<types::ndarray<E, pS1>> const &f);

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 2 && F::value == 1,     // And it is matrix / vect
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f);

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 1 && F::value == 2,     // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f);

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 1 && F::value == 2, // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f);

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 2 && F::value == 1, // And it is vect / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::pshape<long>>>
  dot(E const &e, F const &f);

  /// Matrix / Matrix multiplication

  // The trick is to use the transpose arguments to reflect C order.
  // We want to perform A * B in C order but blas order is F order.
  // So we compute B'A' == (AB)'. As this equality is perform with F order
  // We doesn't have to return a texpr because we want a C order matrice!!
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b);

  template <class E, class pS0, class pS1, class pS2>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2 && std::tuple_size<pS2>::value == 2,
                   types::ndarray<E, pS2>> &
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b, types::ndarray<E, pS2> &c);

  // texpr variants: MT, TM, TT
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a, types::ndarray<E, pS1> const &b);
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::ndarray<E, pS0> const &a, types::numpy_texpr<types::ndarray<E, pS1>> const &b);
  template <class E, class pS0, class pS1>
  std::enable_if_t<is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
                       std::tuple_size<pS1>::value == 2,
                   types::ndarray<E, types::array_tuple<long, 2>>>
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a,
      types::numpy_texpr<types::ndarray<E, pS1>> const &b);

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  std::enable_if_t<types::is_numexpr_arg<E>::value &&
                       types::is_numexpr_arg<F>::value // It is an array_like
                       && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
                           !std::is_same<typename E::dtype, typename F::dtype>::value) &&
                       is_blas_type<typename E::dtype>::value &&
                       is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                              // blas
                       && E::value == 2 && F::value == 2,     // And both are matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, 2>>>
  dot(E const &e, F const &f);

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix multiplication.
  template <class E, class F>
  std::enable_if_t<(!is_blas_type<typename E::dtype>::value ||
                    !is_blas_type<typename F::dtype>::value) &&
                       E::value == 2 && F::value == 2, // And it is matrix / matrix
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, 2>>>
  dot(E const &e, F const &f);

  // N x M where N >= 3 and M == 1
  template <class E, class F>
  std::enable_if_t<(E::value >= 3 && F::value == 1),
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, E::value - 1>>>
  dot(E const &e, F const &f);

  // N x M where N >= 3 and M >= 2
  template <class E, class F>
  std::enable_if_t<(E::value >= 3 && F::value >= 2),
                   types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                                  types::array_tuple<long, E::value - 1>>>
  dot(E const &e, F const &f);

  DEFINE_FUNCTOR(pythonic::numpy, dot);
} // namespace numpy
PYTHONIC_NS_END

#endif
