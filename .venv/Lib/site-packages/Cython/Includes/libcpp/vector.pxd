cdef extern from "<vector>" namespace "std" nogil:
    cdef cppclass vector[T,ALLOCATOR=*]:
        ctypedef T value_type
        ctypedef ALLOCATOR allocator_type

        # these should really be allocator_type.size_type and
        # allocator_type.difference_type to be true to the C++ definition
        # but cython doesn't support deferred access on template arguments
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        bint operator==(const vector&, const vector&)
        bint operator!=(const vector&, const vector&)
        bint operator<(const vector&, const vector&)
        bint operator>(const vector&, const vector&)
        bint operator<=(const vector&, const vector&)
        bint operator>=(const vector&, const vector&)

        cppclass const_iterator
        cppclass iterator:
            iterator() except +
            iterator(iterator&) except +
            T& operator*()
            iterator operator++()
            iterator operator--()
            iterator operator++(int)
            iterator operator--(int)
            iterator operator+(size_type)
            iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(iterator)
            bint operator==(const_iterator)
            bint operator!=(iterator)
            bint operator!=(const_iterator)
            bint operator<(iterator)
            bint operator<(const_iterator)
            bint operator>(iterator)
            bint operator>(const_iterator)
            bint operator<=(iterator)
            bint operator<=(const_iterator)
            bint operator>=(iterator)
            bint operator>=(const_iterator)
        cppclass const_iterator:
            const_iterator() except +
            const_iterator(iterator&) except +
            const_iterator(const_iterator&) except +
            operator=(iterator&) except +
            const T& operator*()
            const_iterator operator++()
            const_iterator operator--()
            const_iterator operator++(int)
            const_iterator operator--(int)
            const_iterator operator+(size_type)
            const_iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(iterator)
            bint operator==(const_iterator)
            bint operator!=(iterator)
            bint operator!=(const_iterator)
            bint operator<(iterator)
            bint operator<(const_iterator)
            bint operator>(iterator)
            bint operator>(const_iterator)
            bint operator<=(iterator)
            bint operator<=(const_iterator)
            bint operator>=(iterator)
            bint operator>=(const_iterator)

        cppclass const_reverse_iterator
        cppclass reverse_iterator:
            reverse_iterator() except +
            reverse_iterator(reverse_iterator&) except +
            T& operator*()
            reverse_iterator operator++()
            reverse_iterator operator--()
            reverse_iterator operator++(int)
            reverse_iterator operator--(int)
            reverse_iterator operator+(size_type)
            reverse_iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(reverse_iterator)
            bint operator!=(const_reverse_iterator)
            bint operator<(reverse_iterator)
            bint operator<(const_reverse_iterator)
            bint operator>(reverse_iterator)
            bint operator>(const_reverse_iterator)
            bint operator<=(reverse_iterator)
            bint operator<=(const_reverse_iterator)
            bint operator>=(reverse_iterator)
            bint operator>=(const_reverse_iterator)
        cppclass const_reverse_iterator:
            const_reverse_iterator() except +
            const_reverse_iterator(reverse_iterator&) except +
            operator=(reverse_iterator&) except +
            const T& operator*()
            const_reverse_iterator operator++()
            const_reverse_iterator operator--()
            const_reverse_iterator operator++(int)
            const_reverse_iterator operator--(int)
            const_reverse_iterator operator+(size_type)
            const_reverse_iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(reverse_iterator)
            bint operator!=(const_reverse_iterator)
            bint operator<(reverse_iterator)
            bint operator<(const_reverse_iterator)
            bint operator>(reverse_iterator)
            bint operator>(const_reverse_iterator)
            bint operator<=(reverse_iterator)
            bint operator<=(const_reverse_iterator)
            bint operator>=(reverse_iterator)
            bint operator>=(const_reverse_iterator)

        # 22.3.11.2, construct/copy/destroy
        vector() except +  # (1)
        vector(size_type) except +  # (3)
        vector(size_type, const T&) except +  # (4)
        # type (vector&) is needed here so that Cython can compile templated definition
        vector& vector[InputIt](InputIt, InputIt) except +  # (5)
        vector(const vector&) except +  # (6)
        # vector& operator=(const vector&)
        void assign[InputIt](InputIt, InputIt) except +
        void assign(size_type, const T&)
        allocator_type get_allocator()

        # iterators
        iterator begin()
        const_iterator const_begin "begin"()
        iterator end()
        const_iterator const_end "end"()
        reverse_iterator rbegin()
        const_reverse_iterator const_rbegin "rbegin"()
        reverse_iterator rend()
        const_reverse_iterator const_rend "rend"()

        const_iterator cbegin()
        const_iterator cend()
        const_reverse_iterator crbegin()
        const_reverse_iterator crend()

        # 22.3.11.3, capacity
        bint empty()
        size_type size()
        size_type max_size()
        size_type capacity()
        void resize(size_type) except +
        void resize(size_type, const T&) except +
        void reserve(size_type) except +
        void shrink_to_fit() except +  # C++11

        # element access
        T& operator[](size_type)
        T& at(size_type) except +
        T& front()
        T& back()

        # 22.3.11.4, data access - C++11
        T* data()
        const T* const_data "data"()

        # 22.3.11.5, modifiers
        T& emplace_back(...) except +  # C++11
        void push_back(const T&) except +
        void pop_back()

        iterator emplace(const_iterator, ...) except +  # C++11
        iterator insert(iterator, const T&) except +
        iterator insert(const_iterator, const T&) except +
        iterator insert(iterator, size_type, const T&) except +
        iterator insert(const_iterator, size_type, const T&) except +
        iterator insert[InputIt](iterator, InputIt, InputIt) except +
        iterator insert[InputIt](const_iterator, InputIt, InputIt) except +
        iterator erase(iterator) except +
        const_iterator erase(const_iterator) except +  # C++11
        iterator erase(iterator, iterator) except +
        iterator erase(const_iterator, const_iterator) except +  # C++11
        void swap(vector&)
        void clear()
