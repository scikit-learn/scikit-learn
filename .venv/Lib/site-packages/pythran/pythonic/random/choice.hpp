#ifndef PYTHONIC_RANDOM_CHOICE_HPP
#define PYTHONIC_RANDOM_CHOICE_HPP

#include "pythonic/include/random/choice.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  namespace details
  {

    template <class Seq>
    std::enable_if_t<types::has_size<Seq>::value, typename Seq::value_type> choice(Seq const &seq)
    {
      auto tmp = seq.begin();
      // std::advance not usable because it requires operator--
      for (long n = random() * seq.size(); n; --n)
        ++tmp;
      return *tmp;
    }

    template <class Seq>
    std::enable_if_t<!types::has_size<Seq>::value, typename Seq::value_type> choice(Seq const &seq)
    {
      using dtype = std::decay_t<typename Seq::value_type>;
      std::vector<dtype, utils::allocator<dtype>> tmp(seq.begin(), seq.end());
      return tmp[long(random() * tmp.size())];
    }
  } // namespace details

  template <class Seq>
  typename Seq::value_type choice(Seq const &seq)
  {
    return details::choice(seq);
  }
} // namespace random
PYTHONIC_NS_END

#endif
