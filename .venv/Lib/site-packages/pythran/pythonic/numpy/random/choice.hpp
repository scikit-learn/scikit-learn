#ifndef PYTHONIC_NUMPY_RANDOM_CHOICE_HPP
#define PYTHONIC_NUMPY_RANDOM_CHOICE_HPP

#include "pythonic/include/numpy/random/choice.hpp"
#include "pythonic/include/numpy/random/generator.hpp"

#include "pythonic/builtins/NotImplementedError.hpp"
#include "pythonic/numpy/random/randint.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>
#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    /***********************************************************
     * Implementation with long as first argument
     **********************************************************/
    template <class pS, class P>
    types::ndarray<long, pS> choice(long max, pS const &shape, bool replace, P const &p)
    {
      if (!replace)
        throw pythonic::builtins::NotImplementedError(
            "Choice without replacement is ! implemented, ask if you want "
            "it");

      types::ndarray<long, pS> result{shape, types::none_type()};
      std::discrete_distribution<long> distribution{p.begin(), p.end()};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    template <class P>
    types::ndarray<long, types::pshape<long>> choice(long max, long size, bool replace, P &&p)
    {
      return choice(max, types::pshape<long>{size}, replace, std::forward<P>(p));
    }

    template <class T>
    auto choice(long max, T &&size) -> decltype(randint(0, max, std::forward<T>(size)))
    {
      return randint(0, max, std::forward<T>(size));
    }

    inline long choice(long max)
    {
      return randint(max);
    }

    /***********************************************************
     * Implementation with array as first argument
     **********************************************************/

    template <class T>
    typename T::dtype choice(T const &a)
    {
      // This is a numpy constraint
      static_assert(T::value == 1, "ValueError: a must be 1-dimensional");

      return a.fast(randint(a.size()));
    }

    template <class T, class pS>
    types::ndarray<typename T::dtype, pS> choice(T const &a, pS const &shape)
    {
      // This is a numpy constraint
      static_assert(T::value == 1, "ValueError: a must be 1-dimensional");

      types::ndarray<typename T::dtype, pS> result{shape, types::none_type()};
      std::uniform_int_distribution<long> distribution{0, a.size() - 1};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return a[distribution(details::generator)]; });
      return result;
    }

    template <class T>
    types::ndarray<typename T::dtype, types::pshape<long>> choice(T &&a, long size)
    {
      return choice(std::forward<T>(a), types::pshape<long>{size});
    }

    template <class T, class pS, class P>
    types::ndarray<typename T::dtype, pS> choice(T const &a, pS const &shape, bool replace,
                                                 P const &p)
    {
      // This is a numpy constraint
      static_assert(T::value == 1, "ValueError: a must be 1-dimensional");

      if (!replace)
        throw pythonic::builtins::NotImplementedError(
            "Choice without replacement is ! implemented, ask if you want "
            "it");

      types::ndarray<typename T::dtype, pS> result{shape, types::none_type()};
      std::discrete_distribution<long> distribution{p.begin(), p.end()};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return a[distribution(details::generator)]; });
      return result;
    }

    template <class T, class P>
    types::ndarray<typename T::dtype, types::pshape<long>> choice(T &&a, long size, bool replace,
                                                                  P &&p)
    {
      return choice(std::forward<T>(a), types::pshape<long>{size}, replace, std::forward<P>(p));
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
