#ifndef PYTHONIC_NUMPY_RANDOM_SHUFFLE_HPP
#define PYTHONIC_NUMPY_RANDOM_SHUFFLE_HPP

#include "pythonic/include/numpy/random/shuffle.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace random
  {
    template <class T>
    types::none_type shuffle(T &seq)
    {
      std::shuffle(seq.begin(), seq.end(), details::generator);
      return builtins::None;
    }
  } // namespace random
} // namespace numpy

PYTHONIC_NS_END

#endif
