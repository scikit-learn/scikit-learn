from libcpp.vector cimport vector

cdef extern from "<span>" namespace "std" nogil:
    # Only Extent = std::dynamic_extent is supported until Cython can also
    # support integer templates. See https://github.com/cython/cython/pull/426
    cdef cppclass span[T]:
        ctypedef T value_type
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        size_t extent

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

        # const_iterator was added in C++23, so leaving it out for now

        span()
        span(T*, size_type) except +  # span[It](It, size_type)
        span(T*, T*) except +  # span[It, End](It, End)
        span(vector&)  # span[U, N](array[T, N]& arr)
        span(span&)

        T& operator[](ssize_t)

        T& back()
        iterator begin()
        T* data()
        bint empty()
        iterator end()
        span[T] first(size_type)
        T& front()
        span[T] last(size_type)
        reverse_iterator rbegin()
        reverse_iterator rend()
        size_type size()
        span[T] subspan(size_type)
        span[T] subspan(size_type, size_type)

    cdef size_t dynamic_extent
