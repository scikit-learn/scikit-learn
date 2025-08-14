from libcpp.string cimport string

cdef extern from "<string_view>" namespace "std::string_view" nogil:
    const size_t npos

cdef extern from "<string_view>" namespace "std" nogil:

    cdef cppclass string_view:
        ctypedef char value_type
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        bint operator==(string_view)
        bint operator!= (string_view)
        bint operator< (string_view)
        bint operator> (string_view)
        bint operator<= (string_view)
        bint operator>= (string_view)

        cppclass const_iterator:
            const_iterator()
            const_iterator(const_iterator&)
            operator=(const_iterator&)
            const value_type& operator*()
            const_iterator operator++()
            const_iterator operator--()
            const_iterator operator++(int)
            const_iterator operator--(int)
            const_iterator operator+(size_type)
            const_iterator operator-(size_type)
            difference_type operator-(const_iterator)
            bint operator==(const_iterator)
            bint operator!=(const_iterator)
            bint operator<(const_iterator)
            bint operator>(const_iterator)
            bint operator<=(const_iterator)
            bint operator>=(const_iterator)

        ctypedef const_iterator iterator

        cppclass const_reverse_iterator:
            const_reverse_iterator()
            const_reverse_iterator(const_reverse_iterator&)
            operator=(const_reverse_iterator&)
            const value_type& operator*()
            const_reverse_iterator operator++()
            const_reverse_iterator operator--()
            const_reverse_iterator operator++(int)
            const_reverse_iterator operator--(int)
            const_reverse_iterator operator+(size_type)
            const_reverse_iterator operator-(size_type)
            difference_type operator-(const_reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(const_reverse_iterator)
            bint operator<(const_reverse_iterator)
            bint operator>(const_reverse_iterator)
            bint operator<=(const_reverse_iterator)
            bint operator>=(const_reverse_iterator)

        ctypedef const_reverse_iterator reverse_iterator

        # Constructors
        string_view()  # (1)
        string_view(string_view s)  # (2)
        string_view(const char* s, size_t n)  # (3)
        string_view(const char* s)  # (4)
        string_view& string_view[InputIt](InputIt, InputIt)  # (5)

        # This should be `operator string_view()` on std::string, but
        # Cython doesn't yet allow declaring it.
        string_view(const string&)

        # Assignment
        string_view& operator=(string_view)

        # Iterators
        const_iterator begin()
        const_iterator end()
        const_iterator cbegin()
        const_iterator cend()

        const_reverse_iterator rbegin()
        const_reverse_iterator rend()
        const_reverse_iterator crbegin()
        const_reverse_iterator crend()

        # Capacity
        size_t size()
        size_t length()
        size_t max_size()
        bint empty()

        # Element access
        const char& operator[](size_t pos)
        const char& at(size_t pos) except +
        const char& front()
        const char& back()

        const char* data()

        # Modifiers
        void remove_prefix(size_t n)
        void remove_suffix(size_t n)
        void swap(string_view& other)

        # Operations
        size_t copy(char* dest, size_t count, size_t pos)
        size_t copy(char* dest, size_t count)

        size_t find(string_view s, size_t pos)
        size_t find(string_view s)
        size_t find(const char* s, size_t pos, size_t n)
        size_t find(const char* s, size_t pos)
        size_t find(const char* s)
        size_t find(char c, size_t pos)
        size_t find(char c)

        size_t rfind(string_view, size_t pos)
        size_t rfind(string_view)
        size_t rfind(const char* s, size_t pos, size_t n)
        size_t rfind(const char* s, size_t pos)
        size_t rfind(const char* s)
        size_t rfind(char c, size_t pos)
        size_t rfind(char c)

        size_t find_first_of(string_view, size_t pos)
        size_t find_first_of(string_view)
        size_t find_first_of(const char* s, size_t pos, size_t n)
        size_t find_first_of(const char* s, size_t pos)
        size_t find_first_of(const char* s)
        size_t find_first_of(char c, size_t pos)
        size_t find_first_of(char c)

        size_t find_first_not_of(string_view s, size_t pos)
        size_t find_first_not_of(string_view s)
        size_t find_first_not_of(const char* s, size_t pos, size_t n)
        size_t find_first_not_of(const char* s, size_t pos)
        size_t find_first_not_of(const char*)
        size_t find_first_not_of(char c, size_t pos)
        size_t find_first_not_of(char c)

        size_t find_last_of(string_view s, size_t pos)
        size_t find_last_of(string_view s)
        size_t find_last_of(const char* s, size_t pos, size_t n)
        size_t find_last_of(const char* s, size_t pos)
        size_t find_last_of(const char* s)
        size_t find_last_of(char c, size_t pos)
        size_t find_last_of(char c)

        size_t find_last_not_of(string_view s, size_t pos)
        size_t find_last_not_of(string_view s)
        size_t find_last_not_of(const char* s, size_t pos, size_t n)
        size_t find_last_not_of(const char* s, size_t pos)
        size_t find_last_not_of(const char* s)
        size_t find_last_not_of(char c, size_t pos)
        size_t find_last_not_of(char c)

        string_view substr(size_t pos, size_t len) except +
        string_view substr(size_t pos) except +
        string_view substr()

        int compare(string_view s)
        int compare(size_t pos, size_t len, string_view s)
        int compare(size_t pos, size_t len, string_view s, size_t subpos, size_t sublen)
        int compare(const char* s)
        int compare(size_t pos, size_t len, const char* s)
        int compare(size_t pos, size_t len, const char* s , size_t n)

        # C++20
        bint starts_with(string_view s)
        bint starts_with(char c)
        bint starts_with(const char* s)

        bint ends_with(string_view s)
        bint ends_with(char c)
        bint ends_with(const char* s)

        # C++23
        bint contains(string_view s)
        bint contains(char c)
        bint contains(const char* s)
