#ifndef PYTHONIC_NUMPY_WHERE_HPP
#define PYTHONIC_NUMPY_WHERE_HPP

#include "pythonic/include/numpy/where.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/copy.hpp"
#include "pythonic/numpy/nonzero.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace impl
  {
    template <class E, class F, class G>
    typename __combined<F, G>::type where(E const &cond, F const &true_, G const &false_)
    {
      if (cond)
        return true_;
      else
        return false_;
    }
  } // namespace impl

#define NUMPY_NARY_FUNC_NAME where
#define NUMPY_NARY_FUNC_SYM impl::where
#define NUMPY_NARY_RESHAPE_MODE reshape_type
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy

namespace types
{

  template <>
  struct Dereferencer<numpy::functor::where> {

    template <class Ts>
    auto operator()(Ts const &iters, std::index_sequence<0, 1, 2>) -> std::enable_if_t<
        types::is_dtype<
            std::remove_cv_t<std::remove_reference_t<decltype(*std::get<0>(iters))>>>::value &&
            types::is_dtype<
                std::remove_cv_t<std::remove_reference_t<decltype(*std::get<1>(iters))>>>::value &&
            types::is_dtype<
                std::remove_cv_t<std::remove_reference_t<decltype(*std::get<2>(iters))>>>::value,
        decltype(numpy::impl::where(*std::get<0>(iters), *std::get<1>(iters), *std::get<2>(iters)))>
    {
      if (*std::get<0>(iters))
        return *std::get<1>(iters);
      else
        return *std::get<2>(iters);
    }

    template <class Ts, size_t... I>
    auto operator()(Ts const &iters, std::index_sequence<I...>, ...)
        -> decltype(numpy::functor::where{}(*std::get<I>(iters)...))
    {
      return numpy::functor::where{}(*std::get<I>(iters)...);
    }
  };
} // namespace types
PYTHONIC_NS_END

#endif
