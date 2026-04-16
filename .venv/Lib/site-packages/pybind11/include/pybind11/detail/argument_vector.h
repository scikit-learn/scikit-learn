/*
    pybind11/detail/argument_vector.h: small_vector-like containers to
    avoid heap allocation of arguments during function call dispatch.

    Copyright (c) Meta Platforms, Inc. and affiliates.

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <pybind11/pytypes.h>

#include "common.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_WARNING_DISABLE_MSVC(4127)

PYBIND11_NAMESPACE_BEGIN(detail)

// Shared implementation utility for our small_vector-like containers.
// We support C++11 and C++14, so we cannot use
// std::variant. Union with the tag packed next to the inline
// array's size is smaller anyway, allowing 1 extra handle of
// inline storage for free. Compare the layouts (1 line per
// size_t/void*, assuming a 64-bit machine):
// With variant, total is N + 2 for N >= 2:
// - variant tag (cannot be packed with the array size)
// - array size (or first pointer of 3 in std::vector)
// - N pointers of inline storage (or 2 remaining pointers of std::vector)
// Custom union, total is N + 1 for N >= 3:
// - variant tag & array size if applicable
// - N pointers of inline storage (or 3 pointers of std::vector)
//
// NOTE: this is a low-level representational convenience; the two
// use cases of this union are materially different and in particular
// have different semantics for inline_array::size. All that is being
// shared is the memory management behavior.
template <typename ArrayT, std::size_t InlineSize, typename VectorT = ArrayT>
union inline_array_or_vector {
    struct inline_array {
        bool is_inline = true;
        std::uint32_t size = 0;
        std::array<ArrayT, InlineSize> arr;
    };
    struct heap_vector {
        bool is_inline = false;
        std::vector<VectorT> vec;

        heap_vector() = default;
        heap_vector(std::size_t count, VectorT value) : vec(count, value) {}
    };

    inline_array iarray;
    heap_vector hvector;

    inline_array_or_vector() : iarray() {}

    ~inline_array_or_vector() {
        if (is_inline()) {
            iarray.~inline_array();
        } else {
            hvector.~heap_vector();
        }
    }

    // Disable copy ctor and assignment.
    inline_array_or_vector(const inline_array_or_vector &) = delete;
    inline_array_or_vector &operator=(const inline_array_or_vector &) = delete;

    inline_array_or_vector(inline_array_or_vector &&rhs) noexcept {
        if (rhs.is_inline()) {
            new (&iarray) inline_array(std::move(rhs.iarray));
        } else {
            new (&hvector) heap_vector(std::move(rhs.hvector));
        }
        assert(is_inline() == rhs.is_inline());
    }

    inline_array_or_vector &operator=(inline_array_or_vector &&rhs) noexcept {
        if (this == &rhs) {
            return *this;
        }

        if (is_inline()) {
            iarray.~inline_array();
        } else {
            hvector.~heap_vector();
        }

        if (rhs.is_inline()) {
            new (&iarray) inline_array(std::move(rhs.iarray));
        } else {
            new (&hvector) heap_vector(std::move(rhs.hvector));
        }
        return *this;
    }

    bool is_inline() const {
        // It is undefined behavior to access the inactive member of a
        // union directly. However, it is well-defined to reinterpret_cast any
        // pointer into a pointer to char and examine it as an array
        // of bytes. See
        // https://dev-discuss.pytorch.org/t/unionizing-for-profit-how-to-exploit-the-power-of-unions-in-c/444#the-memcpy-loophole-4
        bool result = false;
        static_assert(offsetof(inline_array, is_inline) == 0,
                      "untagged union implementation relies on this");
        static_assert(offsetof(heap_vector, is_inline) == 0,
                      "untagged union implementation relies on this");
        std::memcpy(&result, reinterpret_cast<const char *>(this), sizeof(bool));
        return result;
    }
};

template <typename T, std::size_t InlineSize>
struct small_vector {
public:
    small_vector() = default;

    // Disable copy ctor and assignment.
    small_vector(const small_vector &) = delete;
    small_vector &operator=(const small_vector &) = delete;
    small_vector(small_vector &&) noexcept = default;
    small_vector &operator=(small_vector &&) noexcept = default;

    std::size_t size() const {
        if (is_inline()) {
            return m_repr.iarray.size;
        }
        return m_repr.hvector.vec.size();
    }

    T const *data() const {
        if (is_inline()) {
            return m_repr.iarray.arr.data();
        }
        return m_repr.hvector.vec.data();
    }

    T &operator[](std::size_t idx) {
        assert(idx < size());
        if (is_inline()) {
            return m_repr.iarray.arr[idx];
        }
        return m_repr.hvector.vec[idx];
    }

    T const &operator[](std::size_t idx) const {
        assert(idx < size());
        if (is_inline()) {
            return m_repr.iarray.arr[idx];
        }
        return m_repr.hvector.vec[idx];
    }

    void push_back(const T &x) { emplace_back(x); }

    void push_back(T &&x) { emplace_back(std::move(x)); }

    template <typename... Args>
    void emplace_back(Args &&...x) {
        if (is_inline()) {
            auto &ha = m_repr.iarray;
            if (ha.size == InlineSize) {
                move_to_heap_vector_with_reserved_size(InlineSize + 1);
                m_repr.hvector.vec.emplace_back(std::forward<Args>(x)...);
            } else {
                ha.arr[ha.size++] = T(std::forward<Args>(x)...);
            }
        } else {
            m_repr.hvector.vec.emplace_back(std::forward<Args>(x)...);
        }
    }

    void reserve(std::size_t sz) {
        if (is_inline()) {
            if (sz > InlineSize) {
                move_to_heap_vector_with_reserved_size(sz);
            }
        } else {
            reserve_slow_path(sz);
        }
    }

private:
    using repr_type = inline_array_or_vector<T, InlineSize>;
    repr_type m_repr;

    PYBIND11_NOINLINE void move_to_heap_vector_with_reserved_size(std::size_t reserved_size) {
        assert(is_inline());
        auto &ha = m_repr.iarray;
        using heap_vector = typename repr_type::heap_vector;
        heap_vector hv;
        hv.vec.reserve(reserved_size);
        static_assert(std::is_nothrow_move_constructible<T>::value,
                      "this conversion is not exception safe");
        static_assert(std::is_nothrow_move_constructible<heap_vector>::value,
                      "this conversion is not exception safe");
        std::move(ha.arr.begin(), ha.arr.begin() + ha.size, std::back_inserter(hv.vec));
        new (&m_repr.hvector) heap_vector(std::move(hv));
    }

    PYBIND11_NOINLINE void reserve_slow_path(std::size_t sz) { m_repr.hvector.vec.reserve(sz); }

    bool is_inline() const { return m_repr.is_inline(); }
};

// Container to avoid heap allocation for kRequestedInlineSize or fewer booleans.
template <std::size_t kRequestedInlineSize>
struct small_vector<bool, kRequestedInlineSize> {
private:
public:
    small_vector() = default;

    // Disable copy ctor and assignment.
    small_vector(const small_vector &) = delete;
    small_vector &operator=(const small_vector &) = delete;
    small_vector(small_vector &&) noexcept = default;
    small_vector &operator=(small_vector &&) noexcept = default;

    small_vector(std::size_t count, bool value) {
        if (count > kInlineSize) {
            new (&m_repr.hvector) typename repr_type::heap_vector(count, value);
        } else {
            auto &inline_arr = m_repr.iarray;
            inline_arr.arr.fill(value ? static_cast<std::size_t>(-1) : 0);
            inline_arr.size = static_cast<decltype(inline_arr.size)>(count);
        }
    }

    std::size_t size() const {
        if (is_inline()) {
            return m_repr.iarray.size;
        }
        return m_repr.hvector.vec.size();
    }

    void reserve(std::size_t sz) {
        if (is_inline()) {
            if (sz > kInlineSize) {
                move_to_heap_vector_with_reserved_size(sz);
            }
        } else {
            m_repr.hvector.vec.reserve(sz);
        }
    }

    bool operator[](std::size_t idx) const {
        if (is_inline()) {
            return inline_index(idx);
        }
        assert(idx < m_repr.hvector.vec.size());
        return m_repr.hvector.vec[idx];
    }

    void push_back(bool b) {
        if (is_inline()) {
            auto &ha = m_repr.iarray;
            if (ha.size == kInlineSize) {
                move_to_heap_vector_with_reserved_size(kInlineSize + 1);
                push_back_slow_path(b);
            } else {
                assert(ha.size < kInlineSize);
                const auto wbi = word_and_bit_index(ha.size++);
                assert(wbi.word < kWords);
                assert(wbi.bit < kBitsPerWord);
                if (b) {
                    ha.arr[wbi.word] |= (static_cast<std::size_t>(1) << wbi.bit);
                } else {
                    ha.arr[wbi.word] &= ~(static_cast<std::size_t>(1) << wbi.bit);
                }
                assert(operator[](ha.size - 1) == b);
            }
        } else {
            push_back_slow_path(b);
        }
    }

    void set(std::size_t idx, bool value = true) {
        if (is_inline()) {
            auto &ha = m_repr.iarray;
            assert(ha.size < kInlineSize);
            const auto wbi = word_and_bit_index(idx);
            assert(wbi.word < kWords);
            assert(wbi.bit < kBitsPerWord);
            if (value) {
                ha.arr[wbi.word] |= (static_cast<std::size_t>(1) << wbi.bit);
            } else {
                ha.arr[wbi.word] &= ~(static_cast<std::size_t>(1) << wbi.bit);
            }
        } else {
            m_repr.hvector.vec[idx] = value;
        }
    }

    void swap(small_vector &rhs) noexcept { std::swap(m_repr, rhs.m_repr); }

private:
    struct WordAndBitIndex {
        std::size_t word;
        std::size_t bit;
    };

    static WordAndBitIndex word_and_bit_index(std::size_t idx) {
        return WordAndBitIndex{idx / kBitsPerWord, idx % kBitsPerWord};
    }

    bool inline_index(std::size_t idx) const {
        const auto wbi = word_and_bit_index(idx);
        assert(wbi.word < kWords);
        assert(wbi.bit < kBitsPerWord);
        return m_repr.iarray.arr[wbi.word] & (static_cast<std::size_t>(1) << wbi.bit);
    }

    PYBIND11_NOINLINE void move_to_heap_vector_with_reserved_size(std::size_t reserved_size) {
        auto &inline_arr = m_repr.iarray;
        using heap_vector = typename repr_type::heap_vector;
        heap_vector hv;
        hv.vec.reserve(reserved_size);
        for (std::size_t ii = 0; ii < inline_arr.size; ++ii) {
            hv.vec.push_back(inline_index(ii));
        }
        new (&m_repr.hvector) heap_vector(std::move(hv));
    }

    PYBIND11_NOINLINE void push_back_slow_path(bool b) { m_repr.hvector.vec.push_back(b); }

    static constexpr auto kBitsPerWord = 8 * sizeof(std::size_t);
    static constexpr auto kWords = (kRequestedInlineSize + kBitsPerWord - 1) / kBitsPerWord;
    static constexpr auto kInlineSize = kWords * kBitsPerWord;

    using repr_type = inline_array_or_vector<std::size_t, kWords, bool>;
    repr_type m_repr;

    bool is_inline() const { return m_repr.is_inline(); }
};

// Container to avoid heap allocation for N or fewer arguments.
template <size_t N>
using argument_vector = small_vector<handle, N>;

// Container to avoid heap allocation for N or fewer booleans.
template <size_t N>
using args_convert_vector = small_vector<bool, N>;

/// A small_vector of PyObject* that holds references and releases them on destruction.
/// This provides explicit ownership semantics without relying on py::object's
/// destructor, and avoids the need for reinterpret_cast when passing to vectorcall.
template <std::size_t InlineSize>
class ref_small_vector {
public:
    ref_small_vector() = default;

    ~ref_small_vector() {
        for (std::size_t i = 0; i < m_ptrs.size(); ++i) {
            Py_XDECREF(m_ptrs[i]);
        }
    }

    // Disable copy (prevent accidental double-decref)
    ref_small_vector(const ref_small_vector &) = delete;
    ref_small_vector &operator=(const ref_small_vector &) = delete;

    // Move is allowed
    ref_small_vector(ref_small_vector &&other) noexcept : m_ptrs(std::move(other.m_ptrs)) {
        // other.m_ptrs is now empty, so its destructor won't decref anything
    }

    ref_small_vector &operator=(ref_small_vector &&other) noexcept {
        if (this != &other) {
            // Decref our current contents
            for (std::size_t i = 0; i < m_ptrs.size(); ++i) {
                Py_XDECREF(m_ptrs[i]);
            }
            m_ptrs = std::move(other.m_ptrs);
        }
        return *this;
    }

    /// Add a pointer, taking ownership (no incref, will decref on destruction)
    void push_back_steal(PyObject *p) { m_ptrs.push_back(p); }

    /// Add a pointer, borrowing (increfs now, will decref on destruction)
    void push_back_borrow(PyObject *p) {
        Py_XINCREF(p);
        m_ptrs.push_back(p);
    }

    /// Add a null pointer (for PY_VECTORCALL_ARGUMENTS_OFFSET slot)
    void push_back_null() { m_ptrs.push_back(nullptr); }

    void reserve(std::size_t sz) { m_ptrs.reserve(sz); }

    std::size_t size() const { return m_ptrs.size(); }

    PyObject *operator[](std::size_t idx) const { return m_ptrs[idx]; }

    PyObject *const *data() const { return m_ptrs.data(); }

private:
    small_vector<PyObject *, InlineSize> m_ptrs;
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
