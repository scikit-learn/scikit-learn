#ifndef PYTHONIC_UTILS_FWD_HPP
#define PYTHONIC_UTILS_FWD_HPP

#include "pythonic/include/utils/fwd.hpp"

PYTHONIC_NS_BEGIN

namespace utils
{

  template <typename... Types>
  void fwd(Types const &...types)
  {
  }
} // namespace utils
PYTHONIC_NS_END

#endif
