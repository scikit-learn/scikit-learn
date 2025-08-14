from libcpp cimport bool
from libcpp.utility cimport pair
from libc.stddef import ptrdiff_t


cdef extern from "<algorithm>" namespace "std" nogil:
    # Non-modifying sequence operations
    bool all_of[Iter, Pred](Iter first, Iter last, Pred pred) except +
    bool all_of[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +
    bool any_of[Iter, Pred](Iter first, Iter last, Pred pred) except +
    bool any_of[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +
    bool none_of[Iter, Pred](Iter first, Iter last, Pred pred) except +
    bool none_of[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +

    void for_each[Iter, UnaryFunction](Iter first, Iter last, UnaryFunction f) except +  # actually returns f
    void for_each[ExecutionPolicy, Iter, UnaryFunction](ExecutionPolicy&& policy, Iter first, Iter last, UnaryFunction f) except +  # actually returns f

    ptrdiff_t count[Iter, T](Iter first, Iter last, const T& value) except +
    ptrdiff_t count[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    ptrdiff_t count_if[Iter, Pred](Iter first, Iter last, Pred pred) except +
    ptrdiff_t count_if[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +

    pair[Iter1, Iter2] mismatch[Iter1, Iter2](
        Iter1 first1, Iter1 last1, Iter2 first2) except +  # other overloads are tricky
    pair[Iter1, Iter2] mismatch[ExecutionPolicy, Iter1, Iter2](
        ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2) except +

    Iter find[Iter, T](Iter first, Iter last, const T& value) except +
    Iter find[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +

    Iter find_if[Iter, Pred](Iter first, Iter last, Pred pred) except +
    Iter find_if[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +
    Iter find_if_not[Iter, Pred](Iter first, Iter last, Pred pred) except +
    Iter find_if_not[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred pred) except +

    Iter1 find_end[Iter1, Iter2](Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 find_end[ExecutionPolicy, Iter1, Iter2](ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 find_end[Iter1, Iter2, BinaryPred](
        Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +
    Iter1 find_end[ExecutionPolicy, Iter1, Iter2, BinaryPred](
        ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +


    Iter1 find_first_of[Iter1, Iter2](Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 find_first_of[Iter1, Iter2, BinaryPred](
        Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +
    Iter1 find_first_of[ExecutionPolicy, Iter1, Iter2](ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 find_first_of[ExecutionPolicy, Iter1, Iter2, BinaryPred](
        ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +

    Iter adjacent_find[Iter](Iter first, Iter last) except +
    Iter adjacent_find[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    Iter adjacent_find[Iter, BinaryPred](Iter first, Iter last, BinaryPred pred) except +
    Iter adjacent_find[ExecutionPolicy, Iter, BinaryPred](ExecutionPolicy&& policy, Iter first, Iter last, BinaryPred pred) except +

    Iter1 search[Iter1, Iter2](Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 search[ExecutionPolicy, Iter1, Iter2](ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) except +
    Iter1 search[Iter1, Iter2, BinaryPred](
        Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +
    Iter1 search[ExecutionPolicy, Iter1, Iter2, BinaryPred](
        ExecutionPolicy&& policy, Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, BinaryPred pred) except +
    Iter search_n[Iter, Size, T](Iter first1, Iter last1, Size count, const T& value) except +
    Iter search_n[ExecutionPolicy, Iter, Size, T](ExecutionPolicy&& policy, Iter first1, Iter last1, Size count, const T& value) except +
    Iter search_n[Iter, Size, T, BinaryPred](
        Iter first1, Iter last1, Size count, const T& value, BinaryPred pred) except +
    Iter search_n[ExecutionPolicy, Iter, Size, T, BinaryPred](
        ExecutionPolicy&& policy, Iter first1, Iter last1, Size count, const T& value, BinaryPred pred) except +

    # Modifying sequence operations
    OutputIt copy[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt copy[ExecutionPolicy, InputIt, OutputIt](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt copy_if[InputIt, OutputIt, Pred](InputIt first, InputIt last, OutputIt d_first, Pred pred) except +
    OutputIt copy_if[ExecutionPolicy, InputIt, OutputIt, Pred](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, Pred pred) except +
    OutputIt copy_n[InputIt, Size, OutputIt](InputIt first, Size count, OutputIt result) except +
    OutputIt copy_n[ExecutionPolicy, InputIt, Size, OutputIt](ExecutionPolicy&& policy, InputIt first, Size count, OutputIt result) except +
    Iter2 copy_backward[Iter1, Iter2](Iter1 first, Iter1 last, Iter2 d_last) except +
    Iter2 copy_backward[ExecutionPolicy, Iter1, Iter2](ExecutionPolicy&& policy, Iter1 first, Iter1 last, Iter2 d_last) except +

    OutputIt move[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt move[ExecutionPolicy, InputIt, OutputIt](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first) except +
    Iter2 move_backward[Iter1, Iter2](Iter1 first, Iter1 last, Iter2 d_last) except +
    Iter2 move_backward[ExecutionPolicy, Iter1, Iter2](ExecutionPolicy&& policy, Iter1 first, Iter1 last, Iter2 d_last) except +

    void fill[Iter, T](Iter first, Iter last, const T& value) except +
    void fill[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    Iter fill_n[Iter, Size, T](Iter first, Size count, const T& value) except +
    Iter fill_n[ExecutionPolicy, Iter, Size, T](ExecutionPolicy&& policy, Iter first, Size count, const T& value) except +

    OutputIt transform[InputIt, OutputIt, UnaryOp](
        InputIt first1, InputIt last1, OutputIt d_first, UnaryOp unary_op) except +

    # This overload is ambiguous with the next one. We just let C++ disambiguate from the arguments
    # OutputIt transform[ExecutionPolicy, InputIt, OutputIt, UnaryOp](
    #     ExecutionPolicy&& policy, InputIt first1, InputIt last1, OutputIt d_first, UnaryOp unary_op) except +

    OutputIt transform[InputIt1, InputIt2, OutputIt, BinaryOp](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt d_first, BinaryOp binary_op) except +

    OutputIt transform[ExecutionPolicy, InputIt1, InputIt2, OutputIt, BinaryOp](
        ExecutionPolicy&& policy, InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt d_first, BinaryOp binary_op) except +

    void generate[Iter, Generator](Iter first, Iter last, Generator g) except +
    void generate[ExecutionPolicy, Iter, Generator](ExecutionPolicy&& policy, Iter first, Iter last, Generator g) except +
    void generate_n[Iter, Size, Generator](Iter first, Size count, Generator g) except +
    void generate_n[ExecutionPolicy, Iter, Size, Generator](ExecutionPolicy&& policy, Iter first, Size count, Generator g) except +

    Iter remove[Iter, T](Iter first, Iter last, const T& value) except +
    Iter remove[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    Iter remove_if[Iter, UnaryPred](Iter first, Iter last, UnaryPred pred) except +
    Iter remove_if[ExecutionPolicy, Iter, UnaryPred](ExecutionPolicy&& policy, Iter first, Iter last, UnaryPred pred) except +
    OutputIt remove_copy[InputIt, OutputIt, T](InputIt first, InputIt last, OutputIt d_first, const T& value) except +
    OutputIt remove_copy[ExecutionPolicy, InputIt, OutputIt, T](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, const T& value) except +
    OutputIt remove_copy_if[InputIt, OutputIt, UnaryPred](
        InputIt first, InputIt last, OutputIt d_first, UnaryPred pred) except +
    OutputIt remove_copy_if[ExecutionPolicy, InputIt, OutputIt, UnaryPred](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, UnaryPred pred) except +

    void replace[Iter, T](Iter first, Iter last, const T& old_value, const T& new_value) except +
    void replace[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& old_value, const T& new_value) except +
    void replace_if[Iter, UnaryPred, T](Iter first, Iter last, UnaryPred pred, const T& new_value) except +
    OutputIt replace_copy[InputIt, OutputIt, T](
        InputIt first, InputIt last, OutputIt d_first, const T& old_value, const T& new_value) except +
    void replace_if[ExecutionPolicy, Iter, UnaryPred, T](ExecutionPolicy&& policy, Iter first, Iter last, UnaryPred pred, const T& new_value) except +

    OutputIt replace_copy[ExecutionPolicy, InputIt, OutputIt, T](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, const T& old_value, const T& new_value) except +
    OutputIt replace_copy_if[InputIt, OutputIt, UnaryPred, T](
        InputIt first, InputIt last, OutputIt d_first, UnaryPred pred, const T& new_value) except +
    OutputIt replace_copy_if[ExecutionPolicy, InputIt, OutputIt, UnaryPred, T](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, UnaryPred pred, const T& new_value) except +

    void swap[T](T& a, T& b) except +  # array overload also works
    Iter2 swap_ranges[Iter1, Iter2](Iter1 first1, Iter1 last1, Iter2 first2) except +
    void iter_swap[Iter](Iter a, Iter b) except +

    void reverse[Iter](Iter first, Iter last) except +
    void reverse[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    OutputIt reverse_copy[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt reverse_copy[ExecutionPolicy, InputIt, OutputIt](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first) except +

    Iter rotate[Iter](Iter first, Iter n_first, Iter last) except +
    Iter rotate[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter n_first, Iter last) except +
    OutputIt rotate_copy[InputIt, OutputIt](InputIt first, InputIt n_first, InputIt last, OutputIt d_first) except +
    OutputIt rotate_copy[ExecutionPolicy, InputIt, OutputIt](ExecutionPolicy&& policy, InputIt first, InputIt n_first, InputIt last, OutputIt d_first) except +

    Iter unique[Iter](Iter first, Iter last) except +
    Iter unique[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    Iter unique[Iter, BinaryPred](Iter first, Iter last, BinaryPred p) except +
    Iter unique[ExecutionPolicy, Iter, BinaryPred](ExecutionPolicy&& policy, Iter first, Iter last, BinaryPred p) except +
    OutputIt unique_copy[InputIt, OutputIt](InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt unique_copy[ExecutionPolicy, InputIt, OutputIt](ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first) except +
    OutputIt unique_copy[InputIt, OutputIt, BinaryPred](
        InputIt first, InputIt last, OutputIt d_first, BinaryPred pred) except +
    OutputIt unique_copy[ExecutionPolicy, InputIt, OutputIt, BinaryPred](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, BinaryPred pred) except +

    SampleIt sample[PopulationIt, SampleIt, Distance, URBG](PopulationIt first, PopulationIt last, SampleIt out, Distance n, URBG&& g) except +

    # Partitioning operations
    bool is_partitioned[Iter, Pred](Iter first, Iter last, Pred p) except +
    bool is_partitioned[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred p) except +
    Iter partition[Iter, Pred](Iter first, Iter last, Pred p) except +
    Iter partition[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred p) except +
    pair[OutputIt1, OutputIt2] partition_copy[InputIt, OutputIt1, OutputIt2, Pred](
        InputIt first, InputIt last, OutputIt1 d_first_true, OutputIt2 d_first_false, Pred p) except +
    pair[OutputIt1, OutputIt2] partition_copy[ExecutionPolicy, InputIt, OutputIt1, OutputIt2, Pred](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt1 d_first_true, OutputIt2 d_first_false, Pred p) except +

    Iter stable_partition[Iter, Pred](Iter first, Iter last, Pred p) except +
    Iter stable_partition[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred p) except +
    Iter partition_point[Iter, Pred](Iter first, Iter last, Pred p) except +
    Iter partition_point[ExecutionPolicy, Iter, Pred](ExecutionPolicy&& policy, Iter first, Iter last, Pred p) except +

    # Sorting operations
    bool is_sorted[Iter](Iter first, Iter last) except +
    bool is_sorted[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    bool is_sorted[Iter, Compare](Iter first, Iter last, Compare comp) except +
    bool is_sorted[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter last, Compare comp) except +

    Iter is_sorted_until[Iter](Iter first, Iter last) except +
    Iter is_sorted_until[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    Iter is_sorted_until[Iter, Compare](Iter first, Iter last, Compare comp) except +
    Iter is_sorted_until[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter last, Compare comp) except +

    void sort[Iter](Iter first, Iter last) except +
    void sort[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    void sort[Iter, Compare](Iter first, Iter last, Compare comp) except +
    void sort[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter last, Compare comp) except +

    void partial_sort[Iter](Iter first, Iter middle, Iter last) except +
    void partial_sort[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter middle, Iter last) except +
    void partial_sort[Iter, Compare](Iter first, Iter middle, Iter last, Compare comp) except +
    void partial_sort[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter middle, Iter last, Compare comp) except +

    OutputIt partial_sort_copy[InputIt, OutputIt](
        InputIt first, InputIt last, OutputIt d_first, OutputIt d_last) except +
    OutputIt partial_sort_copy[ExecutionPolicy, InputIt, OutputIt](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last) except +
    OutputIt partial_sort_copy[InputIt, OutputIt, Compare](
        InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, Compare comp) except +
    OutputIt partial_sort_copy[ExecutionPolicy, InputIt, OutputIt, Compare](
        ExecutionPolicy&& policy, InputIt first, InputIt last, OutputIt d_first, OutputIt d_last, Compare comp) except +

    void stable_sort[Iter](Iter first, Iter last) except +
    void stable_sort[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    void stable_sort[Iter, Compare](Iter first, Iter last, Compare comp) except +
    void stable_sort[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter last, Compare comp) except +

    void nth_element[Iter](Iter first, Iter nth, Iter last) except +
    void nth_element[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter nth, Iter last) except +
    void nth_element[Iter, Compare](Iter first, Iter nth, Iter last, Compare comp) except +
    void nth_element[ExecutionPolicy, Iter, Compare](ExecutionPolicy&& policy, Iter first, Iter nth, Iter last, Compare comp) except +

    # Binary search operations (on sorted ranges)
    Iter lower_bound[Iter, T](Iter first, Iter last, const T& value) except +
    Iter lower_bound[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    Iter lower_bound[Iter, T, Compare](Iter first, Iter last, const T& value, Compare comp) except +
    Iter lower_bound[ExecutionPolicy, Iter, T, Compare](ExecutionPolicy&& policy, Iter first, Iter last, const T& value, Compare comp) except +

    Iter upper_bound[Iter, T](Iter first, Iter last, const T& value) except +
    Iter upper_bound[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    Iter upper_bound[Iter, T, Compare](Iter first, Iter last, const T& value, Compare comp) except +
    Iter upper_bound[ExecutionPolicy, Iter, T, Compare](ExecutionPolicy&& policy, Iter first, Iter last, const T& value, Compare comp) except +

    bool binary_search[Iter, T](Iter first, Iter last, const T& value) except +
    bool binary_search[ExecutionPolicy, Iter, T](ExecutionPolicy&& policy, Iter first, Iter last, const T& value) except +
    bool binary_search[Iter, T, Compare](Iter first, Iter last, const T& value, Compare comp) except +
    bool binary_search[ExecutionPolicy, Iter, T, Compare](ExecutionPolicy&& policy, Iter first, Iter last, const T& value, Compare comp) except +

    # Other operations on sorted ranges
    OutputIt merge[InputIt1, InputIt2, OutputIt](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) except +
    OutputIt merge[InputIt1, InputIt2, OutputIt, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out, Compare comp) except +

    void inplace_merge[BidirIt](BidirIt first, BidirIt middle, BidirIt last) except +
    void inplace_merge[BidirIt, Compare](BidirIt first, BidirIt middle, BidirIt last, Compare comp) except +

    # Set operations (on sorted ranges)
    bool includes[InputIt1, InputIt2](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) except +

    bool includes[InputIt1, InputIt2, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, Compare comp) except +

    OutputIt set_difference[InputIt1, InputIt2, OutputIt](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) except +

    OutputIt set_difference[InputIt1, InputIt2, OutputIt, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2,
        OutputIt out, Compare comp) except +

    OutputIt set_intersection[InputIt1, InputIt2, OutputIt](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) except +

    OutputIt set_intersection[InputIt1, InputIt2, OutputIt, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out, Compare comp) except +

    OutputIt set_symmetric_difference[InputIt1, InputIt2, OutputIt](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) except +

    OutputIt set_symmetric_difference[InputIt1, InputIt2, OutputIt, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out, Compare comp) except +

    OutputIt set_union[InputIt1, InputIt2, OutputIt](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out) except +

    OutputIt set_union[InputIt1, InputIt2, OutputIt, Compare](
        InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out, Compare comp) except +

    # Heap operations
    void make_heap[Iter](Iter first, Iter last) except +
    void make_heap[Iter, Compare](Iter first, Iter last, Compare comp) except +

    void push_heap[Iter](Iter first, Iter last) except +
    void push_heap[Iter, Compare](Iter first, Iter last, Compare comp) except +

    void pop_heap[Iter](Iter first, Iter last) except +
    void pop_heap[Iter, Compare](Iter first, Iter last, Compare comp) except +

    void sort_heap[Iter](Iter first, Iter last) except +
    void sort_heap[Iter, Compare](Iter first, Iter last, Compare comp) except +

    # Minimum/maximum operations
    Iter min_element[Iter](Iter first, Iter last) except +
    Iter min_element[Iter, Compare](Iter first, Iter last, Compare comp) except +
    Iter min_element[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    Iter max_element[Iter](Iter first, Iter last) except +
    Iter max_element[Iter, Compare](Iter first, Iter last, Compare comp) except +
    Iter max_element[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    pair[T, T] minmax[T](const T& a, const T& b) except +
    pair[T, T] minmax[T, Compare](const T& a, const T& b, Compare comp) except +
    pair[Iter, Iter] minmax_element[Iter](Iter first, Iter last) except +
    pair[Iter, Iter] minmax_element[Iter, Compare](Iter first, Iter last, Compare comp) except +
    pair[Iter, Iter] minmax_element[ExecutionPolicy, Iter](ExecutionPolicy&& policy, Iter first, Iter last) except +
    const T& clamp[T](const T& v, const T& lo, const T& hi) except +
    const T& clamp[T, Compare](const T& v, const T& lo, const T& hi, Compare comp) except +

    # Comparison operations
    bool equal[InputIt1, InputIt2](InputIt1 first1, InputIt1 last1, InputIt2 first2) except +
    bool equal[InputIt1, InputIt2, BinPred](InputIt1 first1, InputIt1 last1, InputIt2 first2, BinPred pred) except +
    # ambiguous with previous overload
    #bool equal[InputIt1, InputIt2](InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) except +
    bool equal[InputIt1, InputIt2, BinPred](InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinPred pred) except +

    bool lexicographical_compare[InputIt1, InputIt2](InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) except +
    # ambiguous with next overload
    #bool lexicographical_compare[InputIt1, InputIt2, ExecutionPolicy](ExecutionPolicy&& policy, InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) except +
    bool lexicographical_compare[InputIt1, InputIt2, Compare](InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, Compare comp) except +

    # Permutation operations
    bool is_permutation[ForwardIt1, ForwardIt2](ForwardIt1 first1, ForwardIt1 last1, ForwardIt2 first2) except +
    bool is_permutation[ForwardIt1, ForwardIt2, BinaryPred](ForwardIt1 first1, ForwardIt1 last1, ForwardIt2 first2, BinaryPred p) except +
    # ambiguous with previous overload
    #bool is_permutation[ForwardIt1, ForwardIt2](ForwardIt1 first1, ForwardIt1 last1, ForwardIt2 first2, ForwardIt2 last2) except +
    bool is_permutation[ForwardIt1, ForwardIt2, BinaryPred](ForwardIt1 first1, ForwardIt1 last1, ForwardIt2 first2, ForwardIt2 last2, BinaryPred p) except +
    bool next_permutation[BidirIt](BidirIt first, BidirIt last) except +
    bool next_permutation[BidirIt, Compare](BidirIt first, BidirIt last, Compare comp) except +
    bool prev_permutation[BidirIt](BidirIt first, BidirIt last) except +
    bool prev_permutation[BidirIt, Compare](BidirIt first, BidirIt last, Compare comp) except +
