cdef extern from "<list>" namespace "std" nogil:
    cdef cppclass list[T,ALLOCATOR=*]:
        ctypedef T value_type
        ctypedef ALLOCATOR allocator_type

        # these should really be allocator_type.size_type and
        # allocator_type.difference_type to be true to the C++ definition
        # but cython doesn't support deferred access on template arguments
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        cppclass const_iterator
        cppclass iterator:
            iterator() except +
            iterator(iterator&) except +
            value_type& operator*()
            iterator operator++()
            iterator operator--()
            iterator operator++(int)
            iterator operator--(int)
            bint operator==(iterator)
            bint operator==(const_iterator)
            bint operator!=(iterator)
            bint operator!=(const_iterator)
        cppclass const_iterator:
            const_iterator() except +
            const_iterator(iterator&) except +
            const_iterator(const_iterator&) except +
            operator=(iterator&) except +
            const value_type& operator*()
            const_iterator operator++()
            const_iterator operator--()
            const_iterator operator++(int)
            const_iterator operator--(int)
            bint operator==(iterator)
            bint operator==(const_iterator)
            bint operator!=(iterator)
            bint operator!=(const_iterator)

        cppclass const_reverse_iterator
        cppclass reverse_iterator:
            reverse_iterator() except +
            reverse_iterator(reverse_iterator&) except +
            value_type& operator*()
            reverse_iterator operator++()
            reverse_iterator operator--()
            reverse_iterator operator++(int)
            reverse_iterator operator--(int)
            bint operator==(reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(reverse_iterator)
            bint operator!=(const_reverse_iterator)
        cppclass const_reverse_iterator:
            const_reverse_iterator() except +
            const_reverse_iterator(reverse_iterator&) except +
            operator=(reverse_iterator&) except +
            const value_type& operator*()
            const_reverse_iterator operator++()
            const_reverse_iterator operator--()
            const_reverse_iterator operator++(int)
            const_reverse_iterator operator--(int)
            bint operator==(reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(reverse_iterator)
            bint operator!=(const_reverse_iterator)

        list() except +
        list(list&) except +
        list(size_t, T&) except +
        #list operator=(list&)
        bint operator==(list&, list&)
        bint operator!=(list&, list&)
        bint operator<(list&, list&)
        bint operator>(list&, list&)
        bint operator<=(list&, list&)
        bint operator>=(list&, list&)
        void assign(size_t, T&) except +
        T& back()
        iterator begin()
        const_iterator const_begin "begin"()
        const_iterator cbegin()
        void clear()
        bint empty()
        iterator end()
        const_iterator const_end "end"()
        const_iterator cend()
        iterator erase(iterator)
        iterator erase(iterator, iterator)
        T& front()
        iterator insert(iterator, T&)
        void insert(iterator, size_t, T&)
        size_t max_size()
        void merge(list&) except +
        #void merge(list&, BinPred)
        void pop_back()
        void pop_front()
        void push_back(T&) except +
        void push_front(T&) except +
        reverse_iterator rbegin()
        const_reverse_iterator const_rbegin "rbegin"()
        const_reverse_iterator crbegin()
        void remove(T&) except +
        #void remove_if(UnPred)
        reverse_iterator rend()
        const_reverse_iterator const_rend "rend"()
        const_reverse_iterator crend()
        void resize(size_t, T&) except +
        void reverse()
        size_t size()
        void sort() except +
        #void sort(BinPred)
        void splice(iterator, list&)
        void splice(iterator, list&, iterator)
        void splice(iterator, list&, iterator, iterator)
        void swap(list&)
        void unique()
        #void unique(BinPred)
