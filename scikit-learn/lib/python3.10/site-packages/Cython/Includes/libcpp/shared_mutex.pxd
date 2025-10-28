from libcpp cimport bool
from libcpp.mutex cimport defer_lock_t, defer_lock, try_to_lock_t, try_to_lock, adopt_lock_t, adopt_lock

cdef extern from "<shared_mutex>" namespace "std" nogil:
    cppclass shared_mutex:
        # may not be present, and we know nothing about it
        cppclass native_handle_type:
            pass

        # We strongly recommend not calling lock with the GIL held (to avoid deadlocks)
        void lock() except+
        bool try_lock()
        void unlock()

        # We strongly recommend not calling lock_shared with the GIL held (to avoid deadlocks)
        void lock_shared() except+
        bool try_lock_shared()
        void unlock_shared()

        native_handle_type native_handle() except+

    cppclass shared_timed_mutex:
        # may not be present, and we know nothing about it.
        # For shared_timed_mutex cppreference doesn't mention this
        cppclass native_handle_type:
            pass

        # We strongly recommend not calling lock with the GIL held (to avoid deadlocks)
        # and moderately recommend not calling the timed lock functions with the GIL either.
        void lock() except+
        bool try_lock()
        bool try_lock_for[T](const T& duration) except+
        bool try_lock_until[T](const T& time_point) except+
        void unlock()

        void lock_shared() except+
        bool try_lock_shared()
        bool try_lock_shared_for[T](const T& duration) except+
        bool try_lock_shared_until[T](const T& time_point) except+
        void unlock_shared()

        native_handle_type native_handle() except+

    cppclass shared_lock[T]:
        ctypedef T mutex_type
        # This covers both the plain regular constructor, the 3 versions with tags
        # and two std::chrono constructors.  The two templated chrono versions stop
        # us from declaring the overloads explicitly.
        shared_lock()
        shared_lock(mutex_type&, ...) except+
        #shared_lock(mutex_type&, defer_lock_t)
        #shared_lock(mutex_type&, try_to_lock_t) except+
        ## this feels like it should be noexcept, but cppreference implies it isn't
        #shared_lock(mutex_type&, adopt_lock_t) except+

        # We strongly recommend not calling lock with the GIL held (to avoid deadlocks)
        void lock() except+
        bool try_lock() except+
        bool try_lock_for[T](const T& duration) except+
        bool try_lock_until[T](const T& time_point) except+
        void unlock() except+

        # We strongly recommend not calling lock_shared with the GIL held (to avoid deadlocks)
        void swap(shared_lock& other)
        # "release" is definitely not the same as unlock. Noted here because
        # DW always makes this mistake and regrets it and wants to save you
        # from the same fate.
        mutex_type* release()

        mutex_type* mutex()
        bool owns_lock()
        bool operator bool()
