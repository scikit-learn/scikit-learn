#ifndef UFUNC_NAME
#error missing UFUNC_NAME
#endif
#ifndef UFUNC_INAME
#error missing UFUNC_INAME
#endif

// clang-format off
#include INCLUDE_FILE(pythonic/include/operator_,UFUNC_INAME)
// clang-format on
#include "pythonic/include/numpy/reduce.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace UFUNC_NAME
  {

    template <class Arg>
    auto reduce(Arg &&arg)
        -> decltype(numpy::reduce<operator_::functor::UFUNC_INAME>(std::forward<Arg>(arg), 0L))
    {
      return numpy::reduce<operator_::functor::UFUNC_INAME>(std::forward<Arg>(arg), 0L);
    }
    template <class... Args>
    auto reduce(Args &&...args) -> std::enable_if_t<
        sizeof...(Args) != 1,
        decltype(numpy::reduce<operator_::functor::UFUNC_INAME>(std::forward<Args>(args)...))>
    {
      return numpy::reduce<operator_::functor::UFUNC_INAME>(std::forward<Args>(args)...);
    }

    DEFINE_FUNCTOR(pythonic::numpy::UFUNC_NAME, reduce);
  } // namespace UFUNC_NAME
} // namespace numpy
PYTHONIC_NS_END
