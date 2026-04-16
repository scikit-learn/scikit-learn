#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_GENERATOR_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_GENERATOR_HPP

#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    namespace details
    {

      /*
       * PCG Random Number Generation for C.
       *
       * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
       *
       * Licensed under the Apache License, Version 2.0 (the "License");
       * you may not use this file except in compliance with the License.
       * You may obtain a copy of the License at
       *
       *     http://www.apache.org/licenses/LICENSE-2.0
       *
       * Unless required by applicable law or agreed to in writing, software
       * distributed under the License is distributed on an "AS IS" BASIS,
       * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
       *implied.
       * See the License for the specific language governing permissions and
       * limitations under the License.
       *
       * For additional information about the PCG random number generation
       *scheme,
       * including its license and other licensing options, visit
       *
       *       http://www.pcg-random.org
       */

      class pcg
      {
        uint64_t state;
        static constexpr uint64_t inc = 0xda3e39cb94b95bdbULL;

      public:
        using result_type = uint32_t;
        static constexpr result_type min()
        {
          return 0;
        }
        static constexpr result_type max()
        {
          return std::numeric_limits<uint32_t>::max();
        }
        friend bool operator==(pcg const &self, pcg const &other)
        {
          return self.state == other.state;
        }
        friend bool operator!=(pcg const &self, pcg const &other)
        {
          return self.state != other.state;
        }

        pcg() : state(0)
        {
        }
        explicit pcg(std::random_device &rd)
        {
          seed(rd());
        }

        void seed(uint64_t value = 0)
        {
          state = value;
          (void)operator()();
        }

        result_type operator()()
        {
          uint64_t oldstate = state;
          state = oldstate * 6364136223846793005ULL + inc;
          uint32_t xorshifted = uint32_t(((oldstate >> 18u) ^ oldstate) >> 27u);
          int rot = oldstate >> 59u;
          return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        }

        void discard(std::size_t n)
        {
          for (std::size_t i = 0; i < n; ++i)
            operator()();
        }

      private:
      };

      std::random_device rd;
      pcg generator(rd);
    } // namespace details
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
