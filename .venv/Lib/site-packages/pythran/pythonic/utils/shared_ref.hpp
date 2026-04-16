#ifndef PYTHONIC_UTILS_SHARED_REF_HPP
#define PYTHONIC_UTILS_SHARED_REF_HPP

#include "pythonic/include/utils/shared_ref.hpp"
#include "pythonic/utils/allocate.hpp"

#include <unordered_map>
#include <utility>
#ifdef _OPENMP
#include <atomic>
#endif

PYTHONIC_NS_BEGIN

namespace utils
{

  /** Light-weight shared_ptr like-class
   *
   *  Unlike std::shared_ptr, it allocates the memory itself using new.
   */
  template <class T>
  template <class... Types>
  shared_ref<T>::memory::memory(Types &&...args)
      : ptr(std::forward<Types>(args)...), count(1), foreign(nullptr)
  {
  }

  template <class T>
  shared_ref<T>::shared_ref(no_memory const &) noexcept : mem(nullptr)
  {
  }

  template <class T>
  shared_ref<T>::shared_ref(no_memory &&) noexcept : mem(nullptr)
  {
  }

  template <class T>
  template <class... Types>
  shared_ref<T>::shared_ref(Types &&...args)
      : mem(new(utils::allocate<memory>(1)) memory(std::forward<Types>(args)...))
  {
  }

  template <class T>
  shared_ref<T>::shared_ref(shared_ref<T> &&p) noexcept : mem(p.mem)
  {
    p.mem = nullptr;
  }

  template <class T>
  shared_ref<T>::shared_ref(shared_ref<T> const &p) noexcept : mem(p.mem)
  {
    if (mem)
      acquire();
  }

  template <class T>
  template <class Tp>
  shared_ref<Tp> shared_ref<T>::recast() noexcept
  {
    shared_ref<Tp> res = no_memory{};
    if (mem)
      acquire();
    res.mem = reinterpret_cast<typename shared_ref<Tp>::memory *>(mem);
    return res;
  }

  template <class T>
  shared_ref<T>::shared_ref(shared_ref<T> &p) noexcept : mem(p.mem)
  {
    if (mem)
      acquire();
  }

  template <class T>
  shared_ref<T>::~shared_ref() noexcept
  {
    dispose();
  }

  template <class T>
  void shared_ref<T>::swap(shared_ref<T> &rhs) noexcept
  {
    using std::swap;
    swap(mem, rhs.mem);
  }

  template <class T>
  shared_ref<T> &shared_ref<T>::operator=(shared_ref<T> p) noexcept
  {
    swap(p);
    return *this;
  }

  template <class T>
  T &shared_ref<T>::operator*() const noexcept
  {
    assert(mem);
    return mem->ptr;
  }

  template <class T>
  T *shared_ref<T>::operator->() const noexcept
  {
    assert(mem);
    return &mem->ptr;
  }

  template <class T>
  bool shared_ref<T>::operator!=(shared_ref<T> const &other) const noexcept
  {
    return mem != other.mem;
  }

  template <class T>
  bool shared_ref<T>::operator==(shared_ref<T> const &other) const noexcept
  {
    return mem == other.mem;
  }

  template <class T>
  void shared_ref<T>::external(extern_type obj_ptr)
  {
    assert(mem);
    mem->foreign = obj_ptr;
  }

  template <class T>
  inline extern_type shared_ref<T>::get_foreign()
  {
    assert(mem);
    return mem->foreign;
  }

  template <class T>
  void shared_ref<T>::dispose()
  {
    if (mem && --mem->count == 0) {
#ifdef ENABLE_PYTHON_MODULE
      if (mem->foreign) {
        Py_DECREF(mem->foreign);
      }
#endif
      mem->~memory();
      utils::deallocate(mem);
      mem = nullptr;
    }
  }

  template <class T>
  void shared_ref<T>::acquire()
  {
    assert(mem);
    ++mem->count;
  }
} // namespace utils
PYTHONIC_NS_END

#endif
