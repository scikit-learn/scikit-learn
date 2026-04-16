#ifndef PYTHONIC_INCLUDE_NUMPY_FLOORDIVIDE_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOORDIVIDE_HPP

#include "pythonic/include//numpy/floor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class Arg0, class Arg1>
    std::complex<std::common_type_t<Arg0, Arg1>> divfloor(std::complex<Arg0> const &arg0,
                                                          std::complex<Arg1> const &arg1)
    {
      return {functor::floor{}(std::real(arg0 / arg1)), 0};
    }

    template <class Arg0, class Arg1>
    auto divfloor(Arg0 const &arg0, Arg1 const &arg1)
        -> std::enable_if_t<(std::is_integral<Arg0>::value && std::is_integral<Arg1>::value),
                            decltype(arg0 / arg1)>
    {
      bool opposite_sign = (arg0 >= 0 && arg1 < 0) || (arg0 < 0 && arg1 >= 0);
      return (arg0 + opposite_sign * (-arg1 + 1)) / arg1;
    }

    template <class Arg0, class Arg1>
    auto divfloor(Arg0 const &arg0, Arg1 const &arg1)
        -> std::enable_if_t<!std::is_integral<Arg0>::value || !std::is_integral<Arg1>::value,
                            decltype(functor::floor{}(arg0 / arg1))>
    {
      return functor::floor{}(arg0 / arg1);
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME floor_divide
#define NUMPY_NARY_FUNC_SYM wrapper::divfloor
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
