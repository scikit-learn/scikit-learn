#ifndef PYTHONIC_INCLUDE_UTILS_FWD_HPP
#define PYTHONIC_INCLUDE_UTILS_FWD_HPP

PYTHONIC_NS_BEGIN

namespace utils
{

  template <typename... Types>
  void fwd(Types const &...types);
}
PYTHONIC_NS_END

#endif
