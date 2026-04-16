#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_CHOICE_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_CHOICE_HPP

#include "pythonic/include/numpy/random/randint.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS, class P>
    types::ndarray<long, pS> choice(long max, pS const &shape, bool replace, P const &p);

    template <class P>
    types::ndarray<long, types::pshape<long>> choice(long max, long size, bool replace, P &&p);

    template <class T>
    auto choice(long max, T &&size) -> decltype(randint(0, max, std::forward<T>(size)));

    long choice(long max);

    template <class T>
    typename T::dtype choice(T const &a);

    template <class T, class pS>
    types::ndarray<typename T::dtype, pS> choice(T const &a, pS const &shape);

    template <class T>
    types::ndarray<typename T::dtype, types::pshape<long>> choice(T &&a, long size);

    template <class T, class pS, class P>
    types::ndarray<typename T::dtype, pS> choice(T const &a, pS const &shape, bool replace,
                                                 P const &p);

    template <class T, class P>
    types::ndarray<typename T::dtype, types::pshape<long>> choice(T &&a, long size, bool replace,
                                                                  P &&p);

    DEFINE_FUNCTOR(pythonic::numpy::random, choice);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
