#ifndef PYTHONIC_TYPES_DYNAMIC_TUPLE_HPP
#define PYTHONIC_TYPES_DYNAMIC_TUPLE_HPP

#include "pythonic/include/types/dynamic_tuple.hpp"

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/nditerator.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/nested_container.hpp"
#include "pythonic/utils/seq.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include <algorithm>
#include <functional>

PYTHONIC_NS_BEGIN
namespace types
{
  template <typename T>
  template <class E>
  long dynamic_tuple<T>::_flat_size(E const &e, utils::int_<1>) const
  {
    return e.size();
  }
  template <class T>
  intptr_t dynamic_tuple<T>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*data));
  }

  template <typename T>
  template <class E, size_t L>
  long dynamic_tuple<T>::_flat_size(E const &e, utils::int_<L>) const
  {
    return e.size() * _flat_size(e.fast(0), utils::int_<L - 1>{});
  }

  template <typename T>
  long dynamic_tuple<T>::flat_size() const
  {
    return _flat_size(*this, utils::int_<value>{});
  }
  template <typename T>
  bool dynamic_tuple<T>::operator==(dynamic_tuple<T> const &other) const
  {
    return size() == other.size() && std::equal(begin(), end(), other.begin());
  }

  template <typename T>
  bool dynamic_tuple<T>::operator!=(dynamic_tuple<T> const &other) const
  {
    return !(*this == other);
  }

  template <typename T>
  bool dynamic_tuple<T>::operator<(dynamic_tuple<T> const &other) const
  {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end(), std::less<T>());
  }

  template <typename T>
  bool dynamic_tuple<T>::operator<=(dynamic_tuple<T> const &other) const
  {
    if (size() == other.size() && std::equal(begin(), end(), other.begin()))
      return true;
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end(), std::less<T>());
  }

  template <typename T>
  bool dynamic_tuple<T>::operator>(dynamic_tuple<T> const &other) const
  {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end(),
                                        std::greater<T>());
  }

  template <typename T>
  bool dynamic_tuple<T>::operator>=(dynamic_tuple<T> const &other) const
  {
    if (size() == other.size() && std::equal(begin(), end(), other.begin()))
      return true;
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end(),
                                        std::greater<T>());
  }

  template <typename T>
  dynamic_tuple<T> dynamic_tuple<T>::operator+(dynamic_tuple<T> const &other) const
  {
    dynamic_tuple<T> result(begin(), end());
    result.data->resize(size() + other.size());
    std::copy(other.begin(), other.end(), result.data->begin() + size());
    return result;
  }
} // namespace types
PYTHONIC_NS_END

namespace std
{
  template <class T>
  size_t hash<pythonic::types::dynamic_tuple<T>>::operator()(
      pythonic::types::dynamic_tuple<T> const &l) const
  {
    std::hash<T> hasher;
    size_t seed = 0x9e3779b9;
    for (auto &&v : l)
      seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
} // namespace std

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/include/utils/fwd.hpp"
#include "pythonic/include/utils/seq.hpp"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

#ifdef Py_LIMITED_API
#define PyTuple_SET_ITEM PyTuple_SetItem
#endif

template <typename T>
PyObject *to_python<types::dynamic_tuple<T>>::convert(types::dynamic_tuple<T> const &t)
{
  size_t N = t.size();
  PyObject *out = PyTuple_New(N);
  for (size_t i = 0; i < N; ++i)
    PyTuple_SET_ITEM(out, i, ::to_python(t.fast(i)));
  return out;
}

#ifdef Py_LIMITED_API
#undef PyTuple_SET_ITEM
#endif
PYTHONIC_NS_END
#endif

#endif
