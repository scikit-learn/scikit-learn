#ifndef PYTHONIC_INCLUDE_UTILS_SHARED_REF_HPP
#define PYTHONIC_INCLUDE_UTILS_SHARED_REF_HPP

#include <unordered_map>
#include <utility>
#ifdef _OPENMP
#define THREAD_SAFE_REF_COUNT
#endif

#ifdef THREAD_SAFE_REF_COUNT
#include <atomic>
#endif
#ifdef ENABLE_PYTHON_MODULE
#include <Python.h>
#endif

PYTHONIC_NS_BEGIN
#ifdef ENABLE_PYTHON_MODULE
using extern_type = PyObject *;
#else
using extern_type = void *;
#endif

#ifdef THREAD_SAFE_REF_COUNT
using atomic_size_t = std::atomic_size_t;
#else
using atomic_size_t = size_t;
#endif

namespace utils
{

  // Force construction of an uninitialized shared_ref
  struct no_memory {
  };

  /** Light-weight shared_ptr like-class
   *
   *  Unlike std::shared_ptr, it allocates the memory itself using new.
   */
  template <class T>
  class shared_ref
  {
  private:
    struct memory {
      T ptr;
      atomic_size_t count;
      extern_type foreign;
      template <class... Types>
      memory(Types &&...args);
    } *mem;

    template <class Tp>
    friend class shared_ref;

  public:
    // Uninitialized ctor
    shared_ref(no_memory const &) noexcept;

    // Uninitialized ctor (rvalue ref)
    shared_ref(no_memory &&) noexcept;

    // Ctor allocate T and forward all arguments to T ctor
    template <class... Types>
    shared_ref(Types &&...args);

    // Move Ctor
    shared_ref(shared_ref<T> &&p) noexcept;

    // Copy Ctor
    shared_ref(shared_ref<T> const &p) noexcept;

    // Copy Ctor, again
    // Without a non-const copy-ctor here, the greedy variadic template ctor
    // takes over
    shared_ref(shared_ref<T> &p) noexcept;

    ~shared_ref() noexcept;

    // Magic swapperator, help for assignment operators
    void swap(shared_ref<T> &rhs) noexcept;

    // Takes by copy so that acquire/release is handle by ctor
    shared_ref<T> &operator=(shared_ref<T> p) noexcept;

    template <class Tp>
    shared_ref<Tp> recast() noexcept;

    T &operator*() const noexcept;

    T *operator->() const noexcept;

    bool operator!=(shared_ref<T> const &other) const noexcept;
    bool operator==(shared_ref<T> const &other) const noexcept;

    // Save pointer to the external object to decref once we doesn't
    // use it anymore
    void external(extern_type obj_ptr);

    extern_type get_foreign();
    bool is_foreign() const;

  private:
    void dispose();
    void acquire();
  };
} // namespace utils
PYTHONIC_NS_END

#endif
