#ifndef PYTHONIC_TYPES_RAW_ARRAY_HPP
#define PYTHONIC_TYPES_RAW_ARRAY_HPP

#include "pythonic/builtins/MemoryError.hpp"
#include "pythonic/include/types/raw_array.hpp"
#include "pythonic/utils/allocate.hpp"

#include <sstream>

PYTHONIC_NS_BEGIN

namespace types
{
  /* Wrapper class to store an array pointer
   *
   * for internal use only, meant to be stored in a shared_ptr
   */
  template <class T>
  raw_array<T>::raw_array() : data(nullptr), external(false)
  {
  }

  template <class T>
  raw_array<T>::raw_array(size_t n) : data(utils::allocate<T>(n)), external(false)
  {
    if (!data) {
      std::ostringstream oss;
      oss << "unable to allocate " << n << " bytes";
      throw types::MemoryError(oss.str());
    }
  }

  template <class T>
  raw_array<T>::raw_array(T *d, ownership o) : data(d), external(o == ownership::external)
  {
  }

  template <class T>
  raw_array<T>::raw_array(raw_array<T> &&d) : data(d.data), external(d.external)
  {
    d.data = nullptr;
  }

  template <class T>
  raw_array<T>::~raw_array()
  {
    if (data && !external)
      utils::deallocate(data);
  }

  template <class T>
  void raw_array<T>::forget()
  {
    external = true;
  }
} // namespace types
PYTHONIC_NS_END

#endif
