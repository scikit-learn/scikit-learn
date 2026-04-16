#ifndef PYTHONIC_INCLUDE_NUMPY_REDUCE_HPP
#define PYTHONIC_INCLUDE_NUMPY_REDUCE_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/types/ndarray.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN
namespace operator_
{
  namespace functor
  {
    struct imax;
    struct imin;
  } // namespace functor
} // namespace operator_

namespace numpy
{

  namespace
  {
    template <class Op, class E, class T>
    struct reduce_result_type_helper {
      using type = typename T::type;
    };
    template <class Op, class E>
    struct reduce_result_type_helper<Op, E, types::none_type> {
      using type = std::conditional_t<
          std::is_integral<typename types::dtype_of<E>::type>::value &&
              (sizeof(typename types::dtype_of<E>::type) < sizeof(long)) &&
              !std::is_same<Op, operator_::functor::imin>::value &&
              !std::is_same<Op, operator_::functor::imax>::value,
          std::conditional_t<
              std::is_same<typename types::dtype_of<E>::type, bool>::value, long,
              std::conditional_t<std::is_signed<typename types::dtype_of<E>::type>::value, long,
                                 unsigned long>>,
          typename types::dtype_of<E>::type>;
    };
    template <class Op, class E, class T = types::none_type>
    using reduce_result_type = typename reduce_result_type_helper<Op, E, T>::type;
  } // namespace

  template <class Op, class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, E>
  reduce(E const &expr, types::none_type _ = types::none_type());

  template <class Op, class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, E>
  reduce(E const &array, long axis);

  template <class Op, class E, class dtype = types::none_type>
  std::enable_if_t<types::is_numexpr_arg<E>::value, reduce_result_type<Op, E, dtype>>
  reduce(E const &expr, types::none_type axis = {}, dtype d = {});

  template <class Op, class E, class dtype = types::none_type>
  reduce_result_type<Op, E> reduce(types::numpy_texpr<E> const &expr, types::none_type axis = {},
                                   dtype d = {})
  {
    return reduce<Op>(expr.arg, axis, d);
  }

  template <class Op, class E, class dtype = types::none_type>
  std::enable_if_t<E::value == 1, reduce_result_type<Op, E, dtype>>
  reduce(E const &array, long axis, dtype d = {}, types::none_type out = {});

  template <class Op, class E, class Out>
  std::enable_if_t<E::value == 1, reduce_result_type<Op, E>>
  reduce(E const &array, long axis, types::none_type dtype, Out &&out);

  namespace
  {
    template <class E, class Op, class dtype = types::none_type>
    using reduced_type =
        types::ndarray<reduce_result_type<Op, E, dtype>, types::array_tuple<long, E::value - 1>>;
  }

  template <class Op, class E, class dtype = types::none_type>
  std::enable_if_t<E::value != 1, reduced_type<E, Op, dtype>>
  reduce(E const &array, long axis, dtype d = {}, types::none_type out = types::none_type());

  template <class Op, class E>
  reduced_type<E, Op> reduce(types::numpy_texpr<E> const &array, long axis,
                             types::none_type dtype = types::none_type(),
                             types::none_type out = types::none_type())
  {
    return reduce<Op>(array.arg, (axis + 1) % 2);
  }

  template <class Op, class E, class Out>
  std::enable_if_t<E::value != 1, reduced_type<E, Op>> reduce(E const &array, long axis,
                                                              types::none_type dtype, Out &&out);
} // namespace numpy
PYTHONIC_NS_END

#endif
