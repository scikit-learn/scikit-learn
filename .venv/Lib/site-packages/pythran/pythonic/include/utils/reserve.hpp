#ifndef PYTHONIC_INCLUDE_UTILS_RESERVE_HPP
#define PYTHONIC_INCLUDE_UTILS_RESERVE_HPP

PYTHONIC_NS_BEGIN

namespace utils
{

  template <class Container, class From>
  void reserve(Container &, From &&); // do nothing unless specialized
}
PYTHONIC_NS_END

#endif
