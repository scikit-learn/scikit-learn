#ifndef PYTHONIC_NUMPY_RANDOM_BYTES_HPP
#define PYTHONIC_NUMPY_RANDOM_BYTES_HPP

#include "pythonic/include/numpy/random/bytes.hpp"
#include "pythonic/include/numpy/random/generator.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <random>
#include <string>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    inline types::str bytes(long length)
    {
      // dummy init + rewrite is faster than reserve && push_back
      types::str result(std::string(length, 0));
      std::uniform_int_distribution<long> distribution{0, 255};
      std::generate(result.chars().begin(), result.chars().end(),
                    [&]() { return static_cast<char>(distribution(details::generator)); });
      return result;
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
