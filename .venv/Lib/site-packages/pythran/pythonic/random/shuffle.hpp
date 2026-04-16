#ifndef PYTHONIC_RANDOM_SHUFFLE_HPP
#define PYTHONIC_RANDOM_SHUFFLE_HPP

#include "pythonic/include/random/shuffle.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/random/random.hpp"
#include "pythonic/utils/functor.hpp"

#include <limits>

PYTHONIC_NS_BEGIN

namespace random
{

  template <class T>
  types::none_type shuffle(T &seq)
  {
    std::shuffle(seq.begin(), seq.end(), __random_generator);
    return builtins::None;
  }

  namespace details
  {
    template <class function>
    struct URG {
      URG(function &&f) : randf(f)
      {
      }

      typedef unsigned result_type;
      static constexpr result_type min()
      {
        return 0;
      }
      /* -1 because of the floor() operation performed by the float->unsigned
       * conversion */
      static constexpr result_type max()
      {
        return std::numeric_limits<result_type>::max() - 1;
      }
      result_type operator()()
      {
        return randf() * std::numeric_limits<result_type>::max();
      }

      function randf;
    };
  } // namespace details

  template <class T, class function>
  types::none_type shuffle(T &seq, function &&randf)
  {
    std::shuffle(seq.begin(), seq.end(), details::URG<function>(std::forward<function>(randf)));
    return builtins::None;
  }
} // namespace random

PYTHONIC_NS_END

#endif
