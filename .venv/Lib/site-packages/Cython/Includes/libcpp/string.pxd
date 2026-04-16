
# deprecated cimport for backwards compatibility:
from libc.string cimport const_char

from libcpp.string_view cimport string_view

cdef extern from "<string>" namespace "std::string" nogil:
    const size_t npos

cdef extern from "<string>" namespace "std" nogil:
    # 21.3.4, numeric conversions
    int stoi(const string& s, size_t* idx, int base) except +
    long stol(const string& s, size_t* idx, int base) except +
    unsigned long stoul(const string& s, size_t* idx, int base) except +
    long long stoll(const string& s, size_t* idx, int base) except +
    unsigned long long stoull(const string& s, size_t* idx, int base) except +
    float stof(const string& s, size_t* idx) except +
    double stod(const string& s, size_t* idx) except +
    long double stold(const string& s, size_t* idx) except +

    # Duplicate declarations are necessary, because extern functions cannot have default arguments values
    int stoi(const string& s, size_t* idx) except +
    long stol(const string& s, size_t* idx) except +
    unsigned long stoul(const string& s, size_t* idx) except +
    long long stoll(const string& s, size_t* idx) except +
    unsigned long long stoull(const string& s, size_t* idx) except +

    int stoi(const string& s) except +
    long stol(const string& s) except +
    unsigned long stoul(const string& s) except +
    long long stoll(const string& s) except +
    unsigned long long stoull(const string& s) except +
    float stof(const string& s) except +
    double stod(const string& s) except +
    long double stold(const string& s) except +

    string to_string(int val) except +
    string to_string(unsigned val) except +
    string to_string(long val) except +
    string to_string(unsigned long val) except +
    string to_string(long long val) except +
    string to_string(unsigned long long val) except +
    string to_string(float val) except +
    string to_string(double val) except +
    string to_string(long double val) except +

    # Cython doesn't know which size type is equivalent to s/size_t on a given platform
    string to_string(size_t val) except +
    string to_string(ssize_t val) except +

    cdef cppclass string:
        ctypedef char value_type

        # these should really be allocator_type.size_type and
        # allocator_type.difference_type to be true to the C++ definition
        # but cython doesn't support deferred access on template arguments
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        bint operator==(const string&)
        bint operator==(const char*)
        bint operator!= (const string&)
        bint operator!= (const char*)
        bint operator< (const string&)
        bint operator< (const char*)
        bint operator> (const string&)
        bint operator> (const char*)
        bint operator<= (const string&)
        bint operator<= (const char*)
        bint operator>= (const string&)
        bint operator>= (const char*)

        string operator+ (const string&) except +
        string operator+ (const char*) except +

        cppclass const_iterator
        cppclass iterator:
            iterator() except +
            iterator(iterator&) except +
            value_type& operator*()
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
            const value_type& operator*()
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
            value_type& operator*()
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
            const value_type& operator*()
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

        # 21.3.2.2, construct/copy/destroy
        string() except +  # (1)
        string(const string& s) except +  # (3)
        string(const string& s, size_t pos) except +  # (5)
        string(const string& s, size_t pos, size_t len) except +  # (6)
        string(const char* s, size_t n) except +  # (9)
        string(const char* s) except +  # (10)
        string(size_t n, char c) except +  # (11)
        string(const string_view& sv) except +
        # type (string&) is needed here so that Cython can compile templated definition
        string& string[InputIt](InputIt, InputIt) except +  # (12)

        #string& operator= (const string&)
        #string& operator= (const char*)
        #string& operator= (char)

        # 21.3.2.3, iterators
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

        # 21.3.2.4, capacity
        size_t size()
        size_t length()
        size_t max_size()
        void resize(size_t, char) except +
        void resize(size_t) except +
        size_t capacity()
        void reserve(size_t) except +
        void shrink_to_fit() except +
        void clear()
        bint empty()

        # 21.3.2.5, element access
        char& operator[](size_t pos)
        char& at(size_t pos) except +

        char& front()
        char& back()

        # 21.3.2.6, modifiers
        string& append(const string& s) except +  # (1)
        string& append(const string& s, size_t subpos, size_t sublen) except +  # (2)
        string& append(const char* s) except +  # (5)
        string& append(const char* s, size_t n) except +  # (6)
        string& append(size_t n, char c) except +  # (7)
        string& append[InputIt](InputIt, InputIt) except +  # (8)

        void push_back(char c) except +

        string& assign(const string& s) except +  # (1)
        string& assign(const string& s, size_t subpos, size_t sublen) except +  # (3)
        string& assign(const char* s, size_t n) except +  # (6)
        string& assign(const char* s) except +  # (7)
        string& assign(size_t n, char c) except +  # (8)
        string& assign[InputIt](InputIt, InputIt) except +  # (9)

        string& insert(size_t pos, const string& s) except +  # (1)
        string& insert(size_t pos, const string& s, size_t subpos, size_t sublen) except +  # (2)
        string& insert(size_t pos, const char* s, size_t n) except +  # (5)
        string& insert(size_t pos, const char* s) except +  # (6)
        string& insert(size_t pos, size_t n, char c) except +  # (7)
        iterator insert(iterator p, char c) except +  # (8)
        iterator insert(iterator p, size_t n, char c) except +  # (9)
        # This method should be templated in InputIt arguments, but then it is not resolved by Cython
        iterator insert(iterator, iterator, iterator) except +  # (10)

        string& erase(size_t pos, size_t len) except +  # (1)
        string& erase(size_t pos) except +  # (1)
        string& erase() except +  # (1)
        iterator erase(iterator p) except +  # (2)
        iterator erase(const_iterator p) except +
        iterator erase(iterator first, iterator last) except +  # (3)
        iterator erase(const_iterator first, const_iterator last) except +

        void pop_back()

        string& replace(size_t pos, size_t len, const string& str) except +  # (1)
        string& replace(size_t pos, size_t len, const string& str, size_t subpos, size_t sublen) except +  # (2)
        string& replace(size_t pos, size_t len, const char* s, size_t n) except +  # (5)
        string& replace(size_t pos, size_t len, const char* s) except +  # (6)
        string& replace(size_t pos, size_t len, size_t n, char c) except +  # (7)
        string& replace(iterator i1, iterator i2, const string& str) except +  # (8)
        string& replace(iterator i1, iterator i2, const char* s, size_t n) except +  # (10)
        string& replace(iterator i1, iterator i2, const char* s) except +  # (11)
        string& replace(iterator i1, iterator i2, size_t n, char c) except +  # (12)
        string& replace(iterator i1, iterator i2, iterator j1, iterator j2) except +  # (13)

        size_t copy(char* s, size_t len, size_t pos) except +
        size_t copy(char* s, size_t len) except +

        void swap(string& other)

        # 21.3.2.7, string operations
        const char* c_str()
        const char* data()

        size_t find(const string& s, size_t pos)
        size_t find(const string& s)
        size_t find(const char* s, size_t pos, size_t n)
        size_t find(const char* s, size_t pos)
        size_t find(const char* s)
        size_t find(char c, size_t pos)
        size_t find(char c)

        size_t rfind(const string&, size_t pos)
        size_t rfind(const string&)
        size_t rfind(const char* s, size_t pos, size_t n)
        size_t rfind(const char* s, size_t pos)
        size_t rfind(const char* s)
        size_t rfind(char c, size_t pos)
        size_t rfind(char c)

        size_t find_first_of(const string&, size_t pos)
        size_t find_first_of(const string&)
        size_t find_first_of(const char* s, size_t pos, size_t n)
        size_t find_first_of(const char* s, size_t pos)
        size_t find_first_of(const char* s)
        size_t find_first_of(char c, size_t pos)
        size_t find_first_of(char c)

        size_t find_first_not_of(const string& s, size_t pos)
        size_t find_first_not_of(const string& s)
        size_t find_first_not_of(const char* s, size_t pos, size_t n)
        size_t find_first_not_of(const char* s, size_t pos)
        size_t find_first_not_of(const char*)
        size_t find_first_not_of(char c, size_t pos)
        size_t find_first_not_of(char c)

        size_t find_last_of(const string& s, size_t pos)
        size_t find_last_of(const string& s)
        size_t find_last_of(const char* s, size_t pos, size_t n)
        size_t find_last_of(const char* s, size_t pos)
        size_t find_last_of(const char* s)
        size_t find_last_of(char c, size_t pos)
        size_t find_last_of(char c)

        size_t find_last_not_of(const string& s, size_t pos)
        size_t find_last_not_of(const string& s)
        size_t find_last_not_of(const char* s, size_t pos, size_t n)
        size_t find_last_not_of(const char* s, size_t pos)
        size_t find_last_not_of(const char* s)
        size_t find_last_not_of(char c, size_t pos)
        size_t find_last_not_of(char c)

        string substr(size_t pos, size_t len) except +
        string substr(size_t pos) except +
        string substr()

        int compare(const string& s)
        int compare(size_t pos, size_t len, const string& s) except +
        int compare(size_t pos, size_t len, const string& s, size_t subpos, size_t sublen) except +
        int compare(const char* s) except +
        int compare(size_t pos, size_t len, const char* s) except +
        int compare(size_t pos, size_t len, const char* s , size_t n) except +

        # C++20
        bint starts_with(char c) except +
        bint starts_with(const char* s)

        bint ends_with(char c) except +
        bint ends_with(const char* s)

        # C++23
        bint contains(char c) except +
        bint contains(const char* s)
