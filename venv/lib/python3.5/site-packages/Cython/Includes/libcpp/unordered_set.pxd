from .utility cimport pair

cdef extern from "<unordered_set>" namespace "std" nogil:
    cdef cppclass unordered_set[T,HASH=*,PRED=*,ALLOCATOR=*]:
        ctypedef T value_type
        cppclass iterator:
            T& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        cppclass reverse_iterator:
            T& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(reverse_iterator)
            bint operator!=(reverse_iterator)
        cppclass const_iterator(iterator):
            pass
        cppclass const_reverse_iterator(reverse_iterator):
            pass
        unordered_set() except +
        unordered_set(unordered_set&) except +
        #unordered_set(key_compare&)
        #unordered_set& operator=(unordered_set&)
        bint operator==(unordered_set&, unordered_set&)
        bint operator!=(unordered_set&, unordered_set&)
        bint operator<(unordered_set&, unordered_set&)
        bint operator>(unordered_set&, unordered_set&)
        bint operator<=(unordered_set&, unordered_set&)
        bint operator>=(unordered_set&, unordered_set&)
        iterator begin()
        const_iterator const_begin "begin"()
        void clear()
        size_t count(T&)
        bint empty()
        iterator end()
        const_iterator const_end "end"()
        pair[iterator, iterator] equal_range(T&)
        #pair[const_iterator, const_iterator] equal_range(T&)
        void erase(iterator)
        void erase(iterator, iterator)
        size_t erase(T&)
        iterator find(T&)
        const_iterator const_find "find"(T&)
        pair[iterator, bint] insert(T&)
        iterator insert(iterator, T&)
        #void insert(input_iterator, input_iterator)
        #key_compare key_comp()
        iterator lower_bound(T&)
        const_iterator const_lower_bound "lower_bound"(T&)
        size_t max_size()
        reverse_iterator rbegin()
        const_reverse_iterator const_rbegin "rbegin"()
        reverse_iterator rend()
        const_reverse_iterator const_rend "rend"()
        size_t size()
        void swap(unordered_set&)
        iterator upper_bound(T&)
        const_iterator const_upper_bound "upper_bound"(T&)
        #value_compare value_comp()
