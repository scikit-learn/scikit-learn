#ifndef PYTHONIC_TYPES_FINFO_HPP
#define PYTHONIC_TYPES_FINFO_HPP

#include "pythonic/include/types/finfo.hpp"

#include "pythonic/types/attr.hpp"

#include <limits>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  T finfo<std::complex<T>>::eps() const
  {
    return std::numeric_limits<T>::epsilon();
  }

  template <class T>
  T finfo<T>::eps() const
  {
    return std::numeric_limits<T>::epsilon();
  }
} // namespace types
PYTHONIC_NS_END

/* pythran attribute system { */
PYTHONIC_NS_BEGIN
namespace builtins
{
  template <class T>
  auto getattr(types::attr::EPS, pythonic::types::finfo<T> const &f) -> decltype(f.eps())
  {
    return f.eps();
  }
} // namespace builtins
PYTHONIC_NS_END
/* } */
#endif
