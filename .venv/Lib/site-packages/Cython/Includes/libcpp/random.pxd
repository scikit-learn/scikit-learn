from libc.stdint cimport uint_fast32_t, uint_fast64_t


cdef extern from "<random>" namespace "std" nogil:
    cdef cppclass random_device:
        ctypedef uint_fast32_t result_type
        random_device() except +
        result_type operator()() except +

    cdef cppclass mt19937:
        ctypedef uint_fast32_t result_type
        mt19937() except +
        mt19937(result_type seed) except +
        result_type operator()() except +
        result_type min() except +
        result_type max() except +
        void discard(size_t z) except +
        void seed(result_type seed) except +

    cdef cppclass mt19937_64:
        ctypedef uint_fast64_t result_type

        mt19937_64() except +
        mt19937_64(result_type seed) except +
        result_type operator()() except +
        result_type min() except +
        result_type max() except +
        void discard(size_t z) except +
        void seed(result_type seed) except +

    cdef cppclass uniform_int_distribution[T]:
        ctypedef T result_type
        uniform_int_distribution() except +
        uniform_int_distribution(T, T) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass uniform_real_distribution[T]:
        ctypedef T result_type
        uniform_real_distribution() except +
        uniform_real_distribution(T, T) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass bernoulli_distribution:
        ctypedef bint result_type
        bernoulli_distribution() except +
        bernoulli_distribution(double) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass binomial_distribution[T]:
        ctypedef T result_type
        binomial_distribution() except +
        binomial_distribution(T, double) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass geometric_distribution[T]:
        ctypedef T result_type
        geometric_distribution() except +
        geometric_distribution(double) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +


    cdef cppclass negative_binomial_distribution[T]:
        ctypedef T result_type
        negative_binomial_distribution() except +
        negative_binomial_distribution(T, double) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass poisson_distribution[T]:
        ctypedef T result_type
        poisson_distribution() except +
        poisson_distribution(double) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass exponential_distribution[T]:
        ctypedef T result_type
        exponential_distribution() except +
        exponential_distribution(result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass gamma_distribution[T]:
        ctypedef T result_type
        gamma_distribution() except +
        gamma_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass weibull_distribution[T]:
        ctypedef T result_type
        weibull_distribution() except +
        weibull_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass extreme_value_distribution[T]:
        ctypedef T result_type
        extreme_value_distribution() except +
        extreme_value_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass normal_distribution[T]:
        ctypedef T result_type
        normal_distribution() except +
        normal_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass lognormal_distribution[T]:
        ctypedef T result_type
        lognormal_distribution() except +
        lognormal_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass chi_squared_distribution[T]:
        ctypedef T result_type
        chi_squared_distribution() except +
        chi_squared_distribution(result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass cauchy_distribution[T]:
        ctypedef T result_type
        cauchy_distribution() except +
        cauchy_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass fisher_f_distribution[T]:
        ctypedef T result_type
        fisher_f_distribution() except +
        fisher_f_distribution(result_type, result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +

    cdef cppclass student_t_distribution[T]:
        ctypedef T result_type
        student_t_distribution() except +
        student_t_distribution(result_type) except +
        result_type operator()[Generator](Generator&) except +
        result_type min() except +
        result_type max() except +
