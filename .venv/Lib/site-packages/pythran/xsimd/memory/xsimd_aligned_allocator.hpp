/***************************************************************************
 * Copyright (c) Johan Mabille, Sylvain Corlay, Wolf Vollprecht and         *
 * Martin Renou                                                             *
 * Copyright (c) QuantStack                                                 *
 * Copyright (c) Serge Guelton                                              *
 *                                                                          *
 * Distributed under the terms of the BSD 3-Clause License.                 *
 *                                                                          *
 * The full license is in the file LICENSE, distributed with this software. *
 ****************************************************************************/

#ifndef XSIMD_ALIGNED_ALLOCATOR_HPP
#define XSIMD_ALIGNED_ALLOCATOR_HPP

#include <algorithm>
#include <cstddef>
#include <utility>
#ifdef _WIN32
#include <malloc.h>
#else
#include <cstdlib>
#endif

#include <cassert>
#include <memory>

#include "../config/xsimd_arch.hpp"

namespace xsimd
{

    /**
     * @class aligned_allocator
     * @brief Allocator for aligned memory
     *
     * The aligned_allocator class template is an allocator that
     * performs memory allocation aligned by the specified value.
     *
     * @tparam T type of objects to allocate.
     * @tparam Align alignment in bytes.
     */
    template <class T, size_t Align>
    class aligned_allocator
    {
    public:
        using value_type = T;
        using pointer = T*;
        using const_pointer = const T*;
        using reference = T&;
        using const_reference = const T&;
        using size_type = size_t;
        using difference_type = ptrdiff_t;

        static constexpr size_t alignment = Align;

        template <class U>
        struct rebind
        {
            using other = aligned_allocator<U, Align>;
        };

        XSIMD_INLINE aligned_allocator() noexcept;
        XSIMD_INLINE aligned_allocator(const aligned_allocator& rhs) noexcept;

        template <class U>
        XSIMD_INLINE aligned_allocator(const aligned_allocator<U, Align>& rhs) noexcept;

        XSIMD_INLINE ~aligned_allocator();

        XSIMD_INLINE pointer address(reference) noexcept;
        XSIMD_INLINE const_pointer address(const_reference) const noexcept;

        XSIMD_INLINE pointer allocate(size_type n, const void* hint = 0);
        XSIMD_INLINE void deallocate(pointer p, size_type n);

        XSIMD_INLINE size_type max_size() const noexcept;
        XSIMD_INLINE size_type size_max() const noexcept;

        template <class U, class... Args>
        XSIMD_INLINE void construct(U* p, Args&&... args);

        template <class U>
        XSIMD_INLINE void destroy(U* p);
    };

    template <class T1, size_t Align1, class T2, size_t Align2>
    XSIMD_INLINE bool operator==(const aligned_allocator<T1, Align1>& lhs,
                                 const aligned_allocator<T2, Align2>& rhs) noexcept;

    template <class T1, size_t Align1, class T2, size_t Align2>
    XSIMD_INLINE bool operator!=(const aligned_allocator<T1, Align1>& lhs,
                                 const aligned_allocator<T2, Align2>& rhs) noexcept;

    XSIMD_INLINE void* aligned_malloc(size_t size, size_t alignment);
    XSIMD_INLINE void aligned_free(void* ptr);

    template <class T>
    XSIMD_INLINE size_t get_alignment_offset(const T* p, size_t size, size_t block_size);

    /************************************
     * aligned_allocator implementation *
     ************************************/

    /**
     * Default constructor.
     */
    template <class T, size_t A>
    XSIMD_INLINE aligned_allocator<T, A>::aligned_allocator() noexcept
    {
    }

    /**
     * Copy constructor.
     */
    template <class T, size_t A>
    XSIMD_INLINE aligned_allocator<T, A>::aligned_allocator(const aligned_allocator&) noexcept
    {
    }

    /**
     * Extended copy constructor.
     */
    template <class T, size_t A>
    template <class U>
    XSIMD_INLINE aligned_allocator<T, A>::aligned_allocator(const aligned_allocator<U, A>&) noexcept
    {
    }

    /**
     * Destructor.
     */
    template <class T, size_t A>
    XSIMD_INLINE aligned_allocator<T, A>::~aligned_allocator()
    {
    }

    /**
     * Returns the actual address of \c r even in presence of overloaded \c operator&.
     * @param r the object to acquire address of.
     * @return the actual address of \c r.
     */
    template <class T, size_t A>
    XSIMD_INLINE auto
    aligned_allocator<T, A>::address(reference r) noexcept -> pointer
    {
        return &r;
    }

    /**
     * Returns the actual address of \c r even in presence of overloaded \c operator&.
     * @param r the object to acquire address of.
     * @return the actual address of \c r.
     */
    template <class T, size_t A>
    XSIMD_INLINE auto
    aligned_allocator<T, A>::address(const_reference r) const noexcept -> const_pointer
    {
        return &r;
    }

    /**
     * Allocates <tt>n * sizeof(T)</tt> bytes of uninitialized memory, aligned by \c A.
     * The alignment may require some extra memory allocation.
     * @param n the number of objects to allocate storage for.
     * @param hint unused parameter provided for standard compliance.
     * @return a pointer to the first byte of a memory block suitably aligned and sufficient to
     * hold an array of \c n objects of type \c T.
     */
    template <class T, size_t A>
    XSIMD_INLINE auto
    aligned_allocator<T, A>::allocate(size_type n, const void*) -> pointer
    {
        pointer res = reinterpret_cast<pointer>(aligned_malloc(sizeof(T) * n, A));
#if defined(_CPPUNWIND) || defined(__cpp_exceptions)
        if (res == nullptr)
            throw std::bad_alloc();
#endif
        return res;
    }

    /**
     * Deallocates the storage referenced by the pointer p, which must be a pointer obtained by
     * an earlier call to allocate(). The argument \c n must be equal to the first argument of the call
     * to allocate() that originally produced \c p; otherwise, the behavior is undefined.
     * @param p pointer obtained from allocate().
     * @param n number of objects earlier passed to allocate().
     */
    template <class T, size_t A>
    XSIMD_INLINE void aligned_allocator<T, A>::deallocate(pointer p, size_type)
    {
        aligned_free(p);
    }

    /**
     * Returns the maximum theoretically possible value of \c n, for which the
     * call allocate(n, 0) could succeed.
     * @return the maximum supported allocated size.
     */
    template <class T, size_t A>
    XSIMD_INLINE auto
    aligned_allocator<T, A>::max_size() const noexcept -> size_type
    {
        return size_type(-1) / sizeof(T);
    }

    /**
     * This method is deprecated, use max_size() instead
     */
    template <class T, size_t A>
    XSIMD_INLINE auto
    aligned_allocator<T, A>::size_max() const noexcept -> size_type
    {
        return size_type(-1) / sizeof(T);
    }

    /**
     * Constructs an object of type \c T in allocated uninitialized memory
     * pointed to by \c p, using placement-new.
     * @param p pointer to allocated uninitialized memory.
     * @param args the constructor arguments to use.
     */
    template <class T, size_t A>
    template <class U, class... Args>
    XSIMD_INLINE void aligned_allocator<T, A>::construct(U* p, Args&&... args)
    {
        new ((void*)p) U(std::forward<Args>(args)...);
    }

    /**
     * Calls the destructor of the object pointed to by \c p.
     * @param p pointer to the object that is going to be destroyed.
     */
    template <class T, size_t A>
    template <class U>
    XSIMD_INLINE void aligned_allocator<T, A>::destroy(U* p)
    {
        p->~U();
    }

    /**
     * @defgroup allocator_comparison Comparison operators
     */

    /**
     * @ingroup allocator_comparison
     * Compares two aligned memory allocator for equality. Since allocators
     * are stateless, return \c true iff <tt>A1 == A2</tt>.
     * @param lhs aligned_allocator to compare.
     * @param rhs aligned_allocator to compare.
     * @return true if the allocators have the same alignment.
     */
    template <class T1, size_t A1, class T2, size_t A2>
    XSIMD_INLINE bool operator==(const aligned_allocator<T1, A1>& lhs,
                                 const aligned_allocator<T2, A2>& rhs) noexcept
    {
        return lhs.alignment == rhs.alignment;
    }

    /**
     * @ingroup allocator_comparison
     * Compares two aligned memory allocator for inequality. Since allocators
     * are stateless, return \c true iff <tt>A1 != A2</tt>.
     * @param lhs aligned_allocator to compare.
     * @param rhs aligned_allocator to compare.
     * @return true if the allocators have different alignments.
     */
    template <class T1, size_t A1, class T2, size_t A2>
    XSIMD_INLINE bool operator!=(const aligned_allocator<T1, A1>& lhs,
                                 const aligned_allocator<T2, A2>& rhs) noexcept
    {
        return !(lhs == rhs);
    }

    /****************************************
     * aligned malloc / free implementation *
     ****************************************/

    namespace detail
    {
        XSIMD_INLINE void* xaligned_malloc(size_t size, size_t alignment)
        {
            assert(((alignment & (alignment - 1)) == 0) && "alignment must be a power of two");
            assert((alignment >= sizeof(void*)) && "alignment must be at least the size of a pointer");
            void* res = nullptr;
#ifdef _WIN32
            res = _aligned_malloc(size, alignment);
#else
            if (posix_memalign(&res, alignment, size) != 0)
            {
                res = nullptr;
            }
#endif
            return res;
        }

        XSIMD_INLINE void xaligned_free(void* ptr)
        {
#ifdef _WIN32
            _aligned_free(ptr);
#else
            free(ptr);
#endif
        }
    }

    XSIMD_INLINE void* aligned_malloc(size_t size, size_t alignment)
    {
        return detail::xaligned_malloc(size, alignment);
    }

    XSIMD_INLINE void aligned_free(void* ptr)
    {
        detail::xaligned_free(ptr);
    }

    template <class T>
    XSIMD_INLINE size_t get_alignment_offset(const T* p, size_t size, size_t block_size)
    {
        // size_t block_size = simd_traits<T>::size;
        if (block_size == 1)
        {
            // The simd_block consists of exactly one scalar so that all
            // elements of the array
            // are "well" aligned.
            return 0;
        }
        else if (size_t(p) & (sizeof(T) - 1))
        {
            // The array is not aligned to the size of a single element, so that
            // no element
            // of the array is well aligned
            return size;
        }
        else
        {
            size_t block_mask = block_size - 1;
            return std::min<size_t>(
                (block_size - ((size_t(p) / sizeof(T)) & block_mask)) & block_mask,
                size);
        }
    }

    template <class T, class A = default_arch>
    using default_allocator = typename std::conditional<A::requires_alignment(),
                                                        aligned_allocator<T, A::alignment()>,
                                                        std::allocator<T>>::type;
}

#endif
