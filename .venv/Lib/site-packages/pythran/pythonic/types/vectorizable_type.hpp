#ifndef PYTHONIC_TYPES_VECTORIZABLE_TYPE_HPP
#define PYTHONIC_TYPES_VECTORIZABLE_TYPE_HPP

#include "pythonic/include/types/vectorizable_type.hpp"

#include "pythonic/include/numpy/bool_.hpp"
#include "pythonic/include/numpy/float32.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/numpy/int16.hpp"
#include "pythonic/include/numpy/int32.hpp"
#include "pythonic/include/numpy/int64.hpp"
#include "pythonic/include/numpy/int8.hpp"
#include "pythonic/include/numpy/uint16.hpp"
#include "pythonic/include/numpy/uint32.hpp"
#include "pythonic/include/numpy/uint64.hpp"
#include "pythonic/include/numpy/uint8.hpp"

PYTHONIC_NS_BEGIN
namespace operator_
{
  namespace functor
  {
    struct mod;
    struct div;
  } // namespace functor
} // namespace operator_

namespace builtins
{
  namespace pythran
  {
    namespace functor
    {
      struct abssqr;
    }
  } // namespace pythran
} // namespace builtins

namespace numpy
{
  namespace functor
  {
    struct arctan2;
    struct angle_in_rad;
    struct asarray_chkfinite;
    struct clip;
    struct copysign;
    struct divide;
    struct fix;
    struct floor_divide;
    struct fmod;
    struct heaviside;
    struct hypot;
    struct isfinite;
    struct isinf;
    struct isnan;
    struct isposinf;
    struct ldexp;
    struct logaddexp;
    struct logaddexp2;
    struct maximum;
    struct minimum;
    struct nan_to_num;
    struct nextafter;
    struct power;
    struct remainder;
    struct rint;
    struct signbit;
    struct spacing;
    struct true_divide;
    struct where;
  } // namespace functor
} // namespace numpy
namespace scipy
{
  namespace special
  {
    namespace functor
    {
      struct binom;
      struct gammaincinv;
      struct hankel1;
      struct hankel2;
      struct jv;
      struct i0;
      struct i0e;
      struct iv;
      struct kv;
      struct yv;
      struct jvp;
      struct ivp;
      struct kvp;
      struct yvp;
      struct ndtr;
      struct ndtri;
      struct spherical_jn;
      struct spherical_yn;
    } // namespace functor
  } // namespace special
} // namespace scipy
namespace types
{
  template <class O, class... Args>
  struct is_vector_op {

    // vectorize everything but these ops. They require special handling for
    // vectorization, and SG did not invest enough time in those
    static const bool value =
        !std::is_same<O, operator_::functor::mod>::value &&
        (!std::is_same<O, operator_::functor::div>::value ||
         utils::all_of<std::is_same<Args, decltype(std::declval<O>()(
                                              std::declval<Args>()...))>::value...>::value) &&
        !std::is_same<O, numpy::functor::logaddexp2>::value &&
        // Return type for generic function should be generic
        !std::is_same<O, numpy::functor::angle_in_rad>::value &&
        !std::is_same<O, numpy::functor::ldexp>::value &&
        !std::is_same<O, numpy::functor::isfinite>::value &&
        !std::is_same<O, numpy::functor::fix>::value &&
        !std::is_same<O, numpy::functor::isinf>::value &&
        !std::is_same<O, numpy::functor::isnan>::value &&
        !std::is_same<O, numpy::functor::isposinf>::value &&
        !std::is_same<O, numpy::functor::rint>::value &&
        !std::is_same<O, numpy::functor::signbit>::value &&
        // conditional processing doesn't permit SIMD
        !std::is_same<O, numpy::functor::nan_to_num>::value &&
        !std::is_same<O, numpy::functor::asarray_chkfinite>::value &&
        !std::is_same<O, numpy::functor::clip>::value &&
        !std::is_same<O, numpy::functor::where>::value &&
        // not supported by xsimd
        !std::is_same<O, numpy::functor::nextafter>::value &&
        !std::is_same<O, numpy::functor::spacing>::value &&
        // not supported for complex numbers
        !(utils::any_of<is_complex<typename dtype_of<Args>::type>::value...>::value &&
          (std::is_same<O, numpy::functor::floor_divide>::value ||
           std::is_same<O, numpy::functor::maximum>::value ||
           std::is_same<O, builtins::pythran::functor::abssqr>::value ||
           std::is_same<O, numpy::functor::minimum>::value)) &&
        // transtyping
        !std::is_same<O, numpy::functor::bool_>::value &&
        !std::is_same<O, numpy::functor::int8>::value &&
        !std::is_same<O, numpy::functor::int16>::value &&
        !std::is_same<O, numpy::functor::int32>::value &&
        !std::is_same<O, numpy::functor::int64>::value &&
        !std::is_same<O, numpy::functor::uint8>::value &&
        !std::is_same<O, numpy::functor::uint16>::value &&
        !std::is_same<O, numpy::functor::uint32>::value &&
        !std::is_same<O, numpy::functor::uint64>::value &&
        !std::is_same<O, numpy::functor::float32>::value &&
        !std::is_same<O, numpy::functor::float64>::value &&
        // not supported for integral numbers
        !(utils::any_of<std::is_integral<typename dtype_of<Args>::type>::value...>::value &&
          (std::is_same<O, numpy::functor::floor_divide>::value ||
           std::is_same<O, numpy::functor::true_divide>::value ||
           std::is_same<O, numpy::functor::divide>::value ||
           std::is_same<O, numpy::functor::arctan2>::value ||
           std::is_same<O, numpy::functor::copysign>::value ||
           std::is_same<O, numpy::functor::logaddexp>::value ||
           std::is_same<O, numpy::functor::power>::value ||
           std::is_same<O, numpy::functor::remainder>::value ||
           std::is_same<O, numpy::functor::hypot>::value ||
           std::is_same<O, numpy::functor::fmod>::value)) &&
        // special functions not in the scope of xsimd
        !std::is_same<O, numpy::functor::heaviside>::value &&
        !std::is_same<O, scipy::special::functor::binom>::value &&
        !std::is_same<O, scipy::special::functor::gammaincinv>::value &&
        !std::is_same<O, scipy::special::functor::hankel1>::value &&
        !std::is_same<O, scipy::special::functor::hankel2>::value &&
        !std::is_same<O, scipy::special::functor::jv>::value &&
        !std::is_same<O, scipy::special::functor::i0>::value &&
        !std::is_same<O, scipy::special::functor::i0e>::value &&
        !std::is_same<O, scipy::special::functor::iv>::value &&
        !std::is_same<O, scipy::special::functor::kv>::value &&
        !std::is_same<O, scipy::special::functor::yv>::value &&
        !std::is_same<O, scipy::special::functor::jvp>::value &&
        !std::is_same<O, scipy::special::functor::ivp>::value &&
        !std::is_same<O, scipy::special::functor::kvp>::value &&
        !std::is_same<O, scipy::special::functor::ndtr>::value &&
        !std::is_same<O, scipy::special::functor::ndtri>::value &&
        !std::is_same<O, scipy::special::functor::yvp>::value &&
        !std::is_same<O, scipy::special::functor::spherical_jn>::value &&
        !std::is_same<O, scipy::special::functor::spherical_yn>::value &&
        //
        true;
  };
} // namespace types
PYTHONIC_NS_END

#endif
