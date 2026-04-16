#ifndef PYTHONIC_NUMPY_INTERP_HPP
#define PYTHONIC_NUMPY_INTERP_HPP

#include "pythonic/include/numpy/argsort.hpp"
#include "pythonic/include/numpy/interp.hpp"
#include "pythonic/include/numpy/remainder.hpp"
#include <pythonic/include/numpy/concatenate.hpp>

#include "pythonic/numpy/argsort.hpp"
#include <pythonic/numpy/concatenate.hpp>
#include <pythonic/numpy/remainder.hpp>

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/interp_core.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T1, class T2, class T3, typename t1, typename t2, typename t3>
  types::ndarray<interp_out_type<T3>, types::pshape<long>> interp(T1 x, T2 xp, T3 fp, t1 _left,
                                                                  t2 _right, t3 _period)
  {
    interp_out_type<T3> left = _left;
    interp_out_type<T3> right = _right;
    double period = _period;
    // Todo: what to do if this condition isn't satisfied? Can't use a statis
    // assert because the size isn't known at compile time.
    assert(xp.template shape<0>() == fp.template shape<0>());
    interp_out_type<T3> outVal(0);

    types::ndarray<interp_out_type<T3>, types::pshape<long>> out = {(long)(x.template shape<0>()),
                                                                    outVal};

    if (period) {
      auto x_rem = pythonic::numpy::functor::remainder{}(x, period);
      auto xp_rem = pythonic::numpy::functor::remainder{}(xp, period);
      auto idx = pythonic::numpy::functor::argsort{}(xp_rem);
      auto xp_sorted = xp_rem[idx];
      auto fp_sorted = fp[idx];

      auto left_pad_xp = types::ndarray<typename T2::dtype, types::pshape<long>>(
          types::pshape<long>(1), xp_sorted[-1] - period);
      auto right_pad_xp = types::ndarray<typename T2::dtype, types::pshape<long>>(
          types::pshape<long>(1), xp_sorted[0] + period);
      auto new_xp = pythonic::numpy::functor::concatenate{}(
          pythonic::types::make_tuple(left_pad_xp, xp_sorted, right_pad_xp));

      auto left_pad_fp = types::ndarray<interp_out_type<T3>, types::pshape<long>>(
          types::pshape<long>(1), fp_sorted[-1]);
      auto right_pad_fp = types::ndarray<interp_out_type<T3>, types::pshape<long>>(
          types::pshape<long>(1), fp_sorted[0]);
      auto new_fp = pythonic::numpy::functor::concatenate{}(
          pythonic::types::make_tuple(left_pad_fp, fp_sorted, right_pad_fp));

      auto lenxp = new_xp.size();
      auto lenx = x_rem.size();
      do_interp(x_rem, new_xp, new_fp, out, lenxp, lenx, 0., 0.);
    } else {
      auto lenxp = xp.size();
      auto lenx = x.size();
      do_interp(x, xp, fp, out, lenxp, lenx, left, right);
    }

    return out;
  }

  // No parameter specified
  template <class T1, class T2, class T3>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, types::none_type period)
  {
    auto _left = fp[0];
    auto _right = fp[-1];
    return interp(x, xp, fp, _left, _right, 0.);
  }

  // left specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, t1 left, types::none_type right, types::none_type period)
  {
    auto _right = fp[-1];
    return interp(x, xp, fp, left, _right, 0.);
  }
  // right specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, t1 right, types::none_type period)
  {
    auto _left = fp[0];
    return interp(x, xp, fp, _left, right, 0.);
  }
  // period specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, t1 period)
  {
    assert(period != 0);
    return interp(x, xp, fp, 0., 0., period);
  }

  // left and right specified
  template <class T1, class T2, class T3, typename t1, typename t2>
  std::enable_if_t<!std::is_arithmetic<T1>::value,
                   types::ndarray<interp_out_type<T3>, types::pshape<long>>>
  interp(T1 x, T2 xp, T3 fp, t1 left, t2 right, types::none_type period)
  {
    return interp(x, xp, fp, left, right, 0.);
  }

  // No parameter specified
  template <class T1, class T2, class T3>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, types::none_type period)
  {
    auto _left = fp[0];
    auto _right = fp[-1];
    auto temp_array = types::ndarray<double, types::pshape<long>>(types::pshape<long>(1), x);
    return interp(temp_array, xp, fp, _left, _right, 0.)[0];
  }

  // left specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, t1 left, types::none_type right, types::none_type period)
  {
    auto _right = fp[-1];
    auto temp_array = types::ndarray<double, types::pshape<long>>(types::pshape<long>(1), x);
    return interp(temp_array, xp, fp, left, _right, 0.)[0];
  }
  // right specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, t1 right, types::none_type period)
  {
    auto _left = fp[0];
    auto temp_array = types::ndarray<double, types::pshape<long>>(types::pshape<long>(1), x);
    return interp(temp_array, xp, fp, _left, right, 0.)[0];
  }
  // period specified
  template <class T1, class T2, class T3, typename t1>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, types::none_type left, types::none_type right, t1 period)
  {
    assert(period != 0);
    auto temp_array = types::ndarray<double, types::pshape<long>>(types::pshape<long>(1), x);
    return interp(temp_array, xp, fp, 0., 0., period)[0];
  }

  // left and right specified,
  template <class T1, class T2, class T3, typename t1, typename t2>
  std::enable_if_t<std::is_arithmetic<T1>::value, interp_out_type<T3>>
  interp(T1 x, T2 xp, T3 fp, t1 left, t2 right, types::none_type period)
  {
    auto temp_array = types::ndarray<double, types::pshape<long>>(types::pshape<long>(1), x);
    return interp(temp_array, xp, fp, left, right, 0.)[0];
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(interp);
} // namespace numpy
PYTHONIC_NS_END

#endif
