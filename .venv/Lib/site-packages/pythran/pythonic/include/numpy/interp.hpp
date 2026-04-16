#ifndef PYTHONIC_INCLUDE_NUMPY_INTERP_HPP
#define PYTHONIC_INCLUDE_NUMPY_INTERP_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  template <class T>
  using interp_out_type =
      std::conditional_t<types::is_complex<typename T::dtype>::value, std::complex<double>, double>;

  // None,None,None
  template <class T1, class T2, class T3>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left = types::none_type{},
         types::none_type right = types::none_type{}, types::none_type period = types::none_type{});

  // left None None
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, t1 left, types::none_type right = types::none_type{},
         types::none_type period = types::none_type{});

  // None right None
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, t1 right,
         types::none_type period = types::none_type{});
  // None None period
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, t1 period);

  // left right None
  template <class T1, class T2, class T3, typename t1, typename t2>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, t1 left, t2 right, types::none_type period = types::none_type{});

  ////////////////////////// NUMERIC TYPES for x.
  template <class T1, class T2, class T3>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left = types::none_type{},
         types::none_type right = types::none_type{}, types::none_type period = types::none_type{});

  // left None None
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, t1 left, types::none_type right = types::none_type{},
         types::none_type period = types::none_type{});

  // None right None
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, t1 right,
         types::none_type period = types::none_type{});

  // None None period
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, t1 period);

  // left right None
  template <class T1, class T2, class T3, typename t1, typename t2>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, t1 left, t2 right, types::none_type period = types::none_type{});

  NUMPY_EXPR_TO_NDARRAY0_DECL(interp);
  DEFINE_FUNCTOR(pythonic::numpy, interp);
} // namespace numpy
PYTHONIC_NS_END

#endif
