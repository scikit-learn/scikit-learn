#ifndef PYTHONIC_INCLUDE_TYPES_FINFO_HPP
#define PYTHONIC_INCLUDE_TYPES_FINFO_HPP

#include "pythonic/include/types/attr.hpp"

#include <limits>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct finfo {
    T eps() const;
  };

  template <class T>
  struct finfo<std::complex<T>> {
    T eps() const;
  };
} // namespace types
PYTHONIC_NS_END

/* pythran attribute system { */
PYTHONIC_NS_BEGIN
namespace builtins
{
  template <class T>
  auto getattr(types::attr::EPS, pythonic::types::finfo<T> const &f) -> decltype(f.eps());
}
PYTHONIC_NS_END
/* } */
#endif
