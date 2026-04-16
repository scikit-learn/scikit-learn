#ifndef PYTHONIC_NUMPY_LINALG_NORM_HPP
#define PYTHONIC_NUMPY_LINALG_NORM_HPP
#include "pythonic/builtins/NotImplementedError.hpp"
#include "pythonic/builtins/pythran/abssqr.hpp"
#include "pythonic/include/numpy/linalg/norm.hpp"
#include "pythonic/numpy/abs.hpp"
#include "pythonic/numpy/asfarray.hpp"
#include "pythonic/numpy/conj.hpp"
#include "pythonic/numpy/inf.hpp"
#include "pythonic/numpy/max.hpp"
#include "pythonic/numpy/min.hpp"
#include "pythonic/numpy/power.hpp"
#include "pythonic/numpy/real.hpp"
#include "pythonic/numpy/sqrt.hpp"
#include "pythonic/numpy/sum.hpp"
PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace linalg
  {
    template <class Array>
    auto norm(Array &&array, types::none_type ord, types::none_type axis)
        -> decltype(pythonic::numpy::functor::sqrt{}(pythonic::numpy::functor::sum{}(
            pythonic::builtins::pythran::functor::abssqr{}(std::forward<Array>(array)))))
    {
      return pythonic::numpy::functor::sqrt{}(pythonic::numpy::functor::sum{}(
          pythonic::builtins::pythran::functor::abssqr{}(std::forward<Array>(array))));
    }

    template <class Array>
    norm_t<Array> norm(Array &&x, double ord, types::none_type)
    {
      switch (std::decay_t<Array>::value) {
      case 1:
        return norm(std::forward<Array>(x), ord, 0L);
      case 2:
        return norm(std::forward<Array>(x), ord, types::array_tuple<long, 2>{{0L, 1L}});
      default:
        throw pythonic::builtins::NotImplementedError("Invalid norm order for matrices.");
      }
    }

    template <class Array>
    norm_t<Array> norm(Array &&x, double ord, long axis)
    {
      auto &&y = pythonic::numpy::functor::asfarray{}(x);
      if (ord == inf)
        return pythonic::numpy::functor::max{}(pythonic::numpy::functor::abs{}(y), axis);
      else if (ord == -inf)
        return pythonic::numpy::functor::min{}(pythonic::numpy::functor::abs{}(y), axis);
      else if (ord == 0.)
        return pythonic::numpy::functor::sum{}(y != 0., axis);
      else if (ord == 1.)
        return pythonic::numpy::functor::sum{}(pythonic::numpy::functor::abs{}(y), axis);
      else if (ord == 2.)
        return pythonic::numpy::functor::sqrt{}(pythonic::numpy::functor::sum{}(
            pythonic::numpy::functor::real{}(pythonic::numpy::functor::conj{}(y)*y), axis));
      else {
        return pythonic::numpy::functor::power{}(
            pythonic::numpy::functor::sum{}(
                pythonic::numpy::functor::power{}(pythonic::numpy::functor::abs{}(y), ord), axis),
            1. / ord);
      }
    }
    template <class Array>
    norm_t<Array> norm(Array &&x, types::none_type ord, double axis)
    {
      return norm(std::forward<Array>(x), 2., axis);
    }
    template <class Array>
    norm_t<Array> norm(Array &&x, double ord, types::array_tuple<long, 1> axis)
    {
      return norm(std::forward<Array>(x), ord, axis[0]);
    }
    template <class Array>
    norm_t<Array> norm(Array &&array, double ord, types::array_tuple<long, 2> axis)
    {
      throw pythonic::builtins::NotImplementedError("We need more dev!");
    }
  } // namespace linalg
} // namespace numpy
PYTHONIC_NS_END

#endif
