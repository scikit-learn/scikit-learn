
cdef extern from "<atomic>" namespace "std" nogil:

    cdef enum memory_order:
        memory_order_relaxed
        memory_order_consume
        memory_order_acquire
        memory_order_release
        memory_order_acq_rel
        memory_order_seq_cst

    cdef cppclass atomic[T]:
        atomic()
        atomic(T)

        bint is_lock_free()
        void store(T)
        void store(T, memory_order)
        T load()
        T load(memory_order)
        T exchange(T)
        T exchange(T, memory_order)

        bint compare_exchange_weak(T&, T, memory_order, memory_order)
        bint compare_exchange_weak(T&, T, memory_order)
        bint compare_exchange_weak(T&, T)
        bint compare_exchange_strong(T&, T, memory_order, memory_order)
        bint compare_exchange_strong(T&, T, memory_order)
        bint compare_exchange_strong(T&, T)

        T fetch_add(T, memory_order)
        T fetch_add(T)
        T fetch_sub(T, memory_order)
        T fetch_sub(T)
        T fetch_and(T, memory_order)
        T fetch_and(T)
        T fetch_or(T, memory_order)
        T fetch_or(T)
        T fetch_xor(T, memory_order)
        T fetch_xor(T)

        T operator++()
        T operator++(int)
        T operator--()
        T operator--(int)

        # modify-in-place operators not yet supported by Cython:
        # T operator+=(T)
        # T operator-=(T)
        # T operator&=(T)
        # T operator|=(T)
        # T operator^=(T)

        bint operator==(atomic[T]&, atomic[T]&)
        bint operator==(atomic[T]&, T&)
        bint operator==(T&, atomic[T]&)
        bint operator!=(atomic[T]&, atomic[T]&)
        bint operator!=(atomic[T]&, T&)
        bint operator!=(T&, atomic[T]&)
