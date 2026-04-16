/*
    pybind11/numpy.h: Basic NumPy support, vectorize() wrapper

    Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"
#include "detail/common.h"
#include "complex.h"
#include "gil_safe_call_once.h"
#include "pytypes.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <typeindex>
#include <utility>
#include <vector>

#ifdef PYBIND11_HAS_SPAN
#    include <span>
#endif

#if defined(PYBIND11_NUMPY_1_ONLY)
#    error "PYBIND11_NUMPY_1_ONLY is no longer supported (see PR #5595)."
#endif

/* This will be true on all flat address space platforms and allows us to reduce the
   whole npy_intp / ssize_t / Py_intptr_t business down to just ssize_t for all size
   and dimension types (e.g. shape, strides, indexing), instead of inflicting this
   upon the library user.
   Note that NumPy 2 now uses ssize_t for `npy_intp` to simplify this. */
static_assert(sizeof(::pybind11::ssize_t) == sizeof(Py_intptr_t), "ssize_t != Py_intptr_t");
static_assert(std::is_signed<Py_intptr_t>::value, "Py_intptr_t must be signed");
// We now can reinterpret_cast between py::ssize_t and Py_intptr_t (MSVC + PyPy cares)

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_WARNING_DISABLE_MSVC(4127)

class dtype; // Forward declaration
class array; // Forward declaration

template <typename>
struct numpy_scalar; // Forward declaration

PYBIND11_NAMESPACE_BEGIN(detail)

template <>
struct handle_type_name<dtype> {
    static constexpr auto name = const_name("numpy.dtype");
};

template <>
struct handle_type_name<array> {
    static constexpr auto name = const_name("numpy.ndarray");
};

template <typename type, typename SFINAE = void>
struct npy_format_descriptor;

/* NumPy 1 proxy (always includes legacy fields) */
struct PyArrayDescr1_Proxy {
    PyObject_HEAD
    PyObject *typeobj;
    char kind;
    char type;
    char byteorder;
    char flags;
    int type_num;
    int elsize;
    int alignment;
    char *subarray;
    PyObject *fields;
    PyObject *names;
};

struct PyArrayDescr_Proxy {
    PyObject_HEAD
    PyObject *typeobj;
    char kind;
    char type;
    char byteorder;
    char _former_flags;
    int type_num;
    /* Additional fields are NumPy version specific. */
};

/* NumPy 2 proxy, including legacy fields */
struct PyArrayDescr2_Proxy {
    PyObject_HEAD
    PyObject *typeobj;
    char kind;
    char type;
    char byteorder;
    char _former_flags;
    int type_num;
    std::uint64_t flags;
    ssize_t elsize;
    ssize_t alignment;
    PyObject *metadata;
    Py_hash_t hash;
    void *reserved_null[2];
    /* The following fields only exist if 0 <= type_num < 2056 */
    char *subarray;
    PyObject *fields;
    PyObject *names;
};

struct PyArray_Proxy {
    PyObject_HEAD
    char *data;
    int nd;
    ssize_t *dimensions;
    ssize_t *strides;
    PyObject *base;
    PyObject *descr;
    int flags;
};

struct PyVoidScalarObject_Proxy {
    PyObject_VAR_HEAD char *obval;
    PyArrayDescr_Proxy *descr;
    int flags;
    PyObject *base;
};

struct numpy_type_info {
    PyObject *dtype_ptr;
    std::string format_str;
};

struct numpy_internals {
    std::unordered_map<std::type_index, numpy_type_info> registered_dtypes;

    numpy_type_info *get_type_info(const std::type_info &tinfo, bool throw_if_missing = true) {
        auto it = registered_dtypes.find(std::type_index(tinfo));
        if (it != registered_dtypes.end()) {
            return &(it->second);
        }
        if (throw_if_missing) {
            pybind11_fail(std::string("NumPy type info missing for ") + tinfo.name());
        }
        return nullptr;
    }

    template <typename T>
    numpy_type_info *get_type_info(bool throw_if_missing = true) {
        return get_type_info(typeid(typename std::remove_cv<T>::type), throw_if_missing);
    }
};

PYBIND11_NOINLINE void load_numpy_internals(numpy_internals *&ptr) {
    ptr = &get_or_create_shared_data<numpy_internals>("_numpy_internals");
}

inline numpy_internals &get_numpy_internals() {
    static numpy_internals *ptr = nullptr;
    if (!ptr) {
        load_numpy_internals(ptr);
    }
    return *ptr;
}

PYBIND11_NOINLINE module_ import_numpy_core_submodule(const char *submodule_name) {
    module_ numpy = module_::import("numpy");
    str version_string = numpy.attr("__version__");
    module_ numpy_lib = module_::import("numpy.lib");
    object numpy_version = numpy_lib.attr("NumpyVersion")(version_string);
    int major_version = numpy_version.attr("major").cast<int>();

    /* `numpy.core` was renamed to `numpy._core` in NumPy 2.0 as it officially
        became a private module. */
    std::string numpy_core_path = major_version >= 2 ? "numpy._core" : "numpy.core";
    return module_::import((numpy_core_path + "." + submodule_name).c_str());
}

template <typename T>
struct same_size {
    template <typename U>
    using as = bool_constant<sizeof(T) == sizeof(U)>;
};

template <typename Concrete>
constexpr int platform_lookup() {
    return -1;
}

// Lookup a type according to its size, and return a value corresponding to the NumPy typenum.
template <typename Concrete, typename T, typename... Ts, typename... Ints>
constexpr int platform_lookup(int I, Ints... Is) {
    return sizeof(Concrete) == sizeof(T) ? I : platform_lookup<Concrete, Ts...>(Is...);
}

struct npy_api {
    // If you change this code, please review `normalized_dtype_num` below.
    enum constants {
        NPY_ARRAY_C_CONTIGUOUS_ = 0x0001,
        NPY_ARRAY_F_CONTIGUOUS_ = 0x0002,
        NPY_ARRAY_OWNDATA_ = 0x0004,
        NPY_ARRAY_FORCECAST_ = 0x0010,
        NPY_ARRAY_ENSUREARRAY_ = 0x0040,
        NPY_ARRAY_ALIGNED_ = 0x0100,
        NPY_ARRAY_WRITEABLE_ = 0x0400,
        NPY_BOOL_ = 0,
        NPY_BYTE_,
        NPY_UBYTE_,
        NPY_SHORT_,
        NPY_USHORT_,
        NPY_INT_,
        NPY_UINT_,
        NPY_LONG_,
        NPY_ULONG_,
        NPY_LONGLONG_,
        NPY_ULONGLONG_,
        NPY_FLOAT_,
        NPY_DOUBLE_,
        NPY_LONGDOUBLE_,
        NPY_CFLOAT_,
        NPY_CDOUBLE_,
        NPY_CLONGDOUBLE_,
        NPY_OBJECT_ = 17,
        NPY_STRING_,
        NPY_UNICODE_,
        NPY_VOID_,
        // Platform-dependent normalization
        NPY_INT8_ = NPY_BYTE_,
        NPY_UINT8_ = NPY_UBYTE_,
        NPY_INT16_ = NPY_SHORT_,
        NPY_UINT16_ = NPY_USHORT_,
        // `npy_common.h` defines the integer aliases. In order, it checks:
        // NPY_BITSOF_LONG, NPY_BITSOF_LONGLONG, NPY_BITSOF_INT, NPY_BITSOF_SHORT, NPY_BITSOF_CHAR
        // and assigns the alias to the first matching size, so we should check in this order.
        NPY_INT32_
        = platform_lookup<std::int32_t, long, int, short>(NPY_LONG_, NPY_INT_, NPY_SHORT_),
        NPY_UINT32_ = platform_lookup<std::uint32_t, unsigned long, unsigned int, unsigned short>(
            NPY_ULONG_, NPY_UINT_, NPY_USHORT_),
        NPY_INT64_
        = platform_lookup<std::int64_t, long, long long, int>(NPY_LONG_, NPY_LONGLONG_, NPY_INT_),
        NPY_UINT64_
        = platform_lookup<std::uint64_t, unsigned long, unsigned long long, unsigned int>(
            NPY_ULONG_, NPY_ULONGLONG_, NPY_UINT_),
        NPY_FLOAT32_ = platform_lookup<float, double, float, long double>(
            NPY_DOUBLE_, NPY_FLOAT_, NPY_LONGDOUBLE_),
        NPY_FLOAT64_ = platform_lookup<double, double, float, long double>(
            NPY_DOUBLE_, NPY_FLOAT_, NPY_LONGDOUBLE_),
        NPY_COMPLEX64_
        = platform_lookup<std::complex<float>,
                          std::complex<double>,
                          std::complex<float>,
                          std::complex<long double>>(NPY_DOUBLE_, NPY_FLOAT_, NPY_LONGDOUBLE_),
        NPY_COMPLEX128_
        = platform_lookup<std::complex<double>,
                          std::complex<double>,
                          std::complex<float>,
                          std::complex<long double>>(NPY_DOUBLE_, NPY_FLOAT_, NPY_LONGDOUBLE_),
        NPY_CHAR_ = std::is_signed<char>::value ? NPY_BYTE_ : NPY_UBYTE_,
    };

    unsigned int PyArray_RUNTIME_VERSION_;

    struct PyArray_Dims {
        Py_intptr_t *ptr;
        int len;
    };

    static npy_api &get() {
        PYBIND11_CONSTINIT static gil_safe_call_once_and_store<npy_api> storage;
        return storage.call_once_and_store_result(lookup).get_stored();
    }

    bool PyArray_Check_(PyObject *obj) const {
        return PyObject_TypeCheck(obj, PyArray_Type_) != 0;
    }
    bool PyArrayDescr_Check_(PyObject *obj) const {
        return PyObject_TypeCheck(obj, PyArrayDescr_Type_) != 0;
    }

    unsigned int (*PyArray_GetNDArrayCFeatureVersion_)();
    PyObject *(*PyArray_DescrFromType_)(int);
    PyObject *(*PyArray_TypeObjectFromType_)(int);
    PyObject *(*PyArray_NewFromDescr_)(PyTypeObject *,
                                       PyObject *,
                                       int,
                                       Py_intptr_t const *,
                                       Py_intptr_t const *,
                                       void *,
                                       int,
                                       PyObject *);
    // Unused. Not removed because that affects ABI of the class.
    PyObject *(*PyArray_DescrNewFromType_)(int);
    int (*PyArray_CopyInto_)(PyObject *, PyObject *);
    PyObject *(*PyArray_NewCopy_)(PyObject *, int);
    PyTypeObject *PyArray_Type_;
    PyTypeObject *PyVoidArrType_Type_;
    PyTypeObject *PyArrayDescr_Type_;
    PyObject *(*PyArray_DescrFromScalar_)(PyObject *);
    PyObject *(*PyArray_Scalar_)(void *, PyObject *, PyObject *);
    void (*PyArray_ScalarAsCtype_)(PyObject *, void *);
    PyObject *(*PyArray_FromAny_)(PyObject *, PyObject *, int, int, int, PyObject *);
    int (*PyArray_DescrConverter_)(PyObject *, PyObject **);
    bool (*PyArray_EquivTypes_)(PyObject *, PyObject *);
    PyObject *(*PyArray_Squeeze_)(PyObject *);
    // Unused. Not removed because that affects ABI of the class.
    int (*PyArray_SetBaseObject_)(PyObject *, PyObject *);
    PyObject *(*PyArray_Resize_)(PyObject *, PyArray_Dims *, int, int);
    PyObject *(*PyArray_Newshape_)(PyObject *, PyArray_Dims *, int);
    PyObject *(*PyArray_View_)(PyObject *, PyObject *, PyObject *);

private:
    enum functions {
        API_PyArray_GetNDArrayCFeatureVersion = 211,
        API_PyArray_Type = 2,
        API_PyArrayDescr_Type = 3,
        API_PyVoidArrType_Type = 39,
        API_PyArray_DescrFromType = 45,
        API_PyArray_TypeObjectFromType = 46,
        API_PyArray_DescrFromScalar = 57,
        API_PyArray_Scalar = 60,
        API_PyArray_ScalarAsCtype = 62,
        API_PyArray_FromAny = 69,
        API_PyArray_Resize = 80,
        // CopyInto was slot 82 and 50 was effectively an alias. NumPy 2 removed 82.
        API_PyArray_CopyInto = 50,
        API_PyArray_NewCopy = 85,
        API_PyArray_NewFromDescr = 94,
        API_PyArray_DescrNewFromType = 96,
        API_PyArray_Newshape = 135,
        API_PyArray_Squeeze = 136,
        API_PyArray_View = 137,
        API_PyArray_DescrConverter = 174,
        API_PyArray_EquivTypes = 182,
        API_PyArray_SetBaseObject = 282
    };

    static npy_api lookup() {
        module_ m = detail::import_numpy_core_submodule("multiarray");
        auto c = m.attr("_ARRAY_API");
        void **api_ptr = (void **) PyCapsule_GetPointer(c.ptr(), nullptr);
        if (api_ptr == nullptr) {
            raise_from(PyExc_SystemError, "FAILURE obtaining numpy _ARRAY_API pointer.");
            throw error_already_set();
        }
        npy_api api;
#define DECL_NPY_API(Func) api.Func##_ = (decltype(api.Func##_)) api_ptr[API_##Func];
        DECL_NPY_API(PyArray_GetNDArrayCFeatureVersion);
        api.PyArray_RUNTIME_VERSION_ = api.PyArray_GetNDArrayCFeatureVersion_();
        if (api.PyArray_RUNTIME_VERSION_ < 0x7) {
            pybind11_fail("pybind11 numpy support requires numpy >= 1.7.0");
        }
        DECL_NPY_API(PyArray_Type);
        DECL_NPY_API(PyVoidArrType_Type);
        DECL_NPY_API(PyArrayDescr_Type);
        DECL_NPY_API(PyArray_DescrFromType);
        DECL_NPY_API(PyArray_TypeObjectFromType);
        DECL_NPY_API(PyArray_DescrFromScalar);
        DECL_NPY_API(PyArray_Scalar);
        DECL_NPY_API(PyArray_ScalarAsCtype);
        DECL_NPY_API(PyArray_FromAny);
        DECL_NPY_API(PyArray_Resize);
        DECL_NPY_API(PyArray_CopyInto);
        DECL_NPY_API(PyArray_NewCopy);
        DECL_NPY_API(PyArray_NewFromDescr);
        DECL_NPY_API(PyArray_DescrNewFromType);
        DECL_NPY_API(PyArray_Newshape);
        DECL_NPY_API(PyArray_Squeeze);
        DECL_NPY_API(PyArray_View);
        DECL_NPY_API(PyArray_DescrConverter);
        DECL_NPY_API(PyArray_EquivTypes);
        DECL_NPY_API(PyArray_SetBaseObject);

#undef DECL_NPY_API
        return api;
    }
};

template <typename T>
struct is_complex : std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

template <typename T, typename = void>
struct npy_format_descriptor_name;

template <typename T>
struct npy_format_descriptor_name<T, enable_if_t<std::is_integral<T>::value>> {
    static constexpr auto name = const_name<std::is_same<T, bool>::value>(
        const_name("numpy.bool"),
        const_name<std::is_signed<T>::value>("numpy.int", "numpy.uint")
            + const_name<sizeof(T) * 8>());
};

template <typename T>
struct npy_format_descriptor_name<T, enable_if_t<std::is_floating_point<T>::value>> {
    static constexpr auto name = const_name < std::is_same<T, float>::value
                                 || std::is_same<T, const float>::value
                                 || std::is_same<T, double>::value
                                 || std::is_same<T, const double>::value
                                        > (const_name("numpy.float") + const_name<sizeof(T) * 8>(),
                                           const_name("numpy.longdouble"));
};

template <typename T>
struct npy_format_descriptor_name<T, enable_if_t<is_complex<T>::value>> {
    static constexpr auto name = const_name < std::is_same<typename T::value_type, float>::value
                                 || std::is_same<typename T::value_type, const float>::value
                                 || std::is_same<typename T::value_type, double>::value
                                 || std::is_same<typename T::value_type, const double>::value
                                        > (const_name("numpy.complex")
                                               + const_name<sizeof(typename T::value_type) * 16>(),
                                           const_name("numpy.clongdouble"));
};

template <typename T>
struct numpy_scalar_info {};

#define PYBIND11_NUMPY_SCALAR_IMPL(ctype_, typenum_)                                              \
    template <>                                                                                   \
    struct numpy_scalar_info<ctype_> {                                                            \
        static constexpr auto name = npy_format_descriptor_name<ctype_>::name;                    \
        static constexpr int typenum = npy_api::typenum_##_;                                      \
    }

// boolean type
PYBIND11_NUMPY_SCALAR_IMPL(bool, NPY_BOOL);

// character types
PYBIND11_NUMPY_SCALAR_IMPL(char, NPY_CHAR);
PYBIND11_NUMPY_SCALAR_IMPL(signed char, NPY_BYTE);
PYBIND11_NUMPY_SCALAR_IMPL(unsigned char, NPY_UBYTE);

// signed integer types
PYBIND11_NUMPY_SCALAR_IMPL(std::int16_t, NPY_INT16);
PYBIND11_NUMPY_SCALAR_IMPL(std::int32_t, NPY_INT32);
PYBIND11_NUMPY_SCALAR_IMPL(std::int64_t, NPY_INT64);

// unsigned integer types
PYBIND11_NUMPY_SCALAR_IMPL(std::uint16_t, NPY_UINT16);
PYBIND11_NUMPY_SCALAR_IMPL(std::uint32_t, NPY_UINT32);
PYBIND11_NUMPY_SCALAR_IMPL(std::uint64_t, NPY_UINT64);

// floating point types
PYBIND11_NUMPY_SCALAR_IMPL(float, NPY_FLOAT);
PYBIND11_NUMPY_SCALAR_IMPL(double, NPY_DOUBLE);
PYBIND11_NUMPY_SCALAR_IMPL(long double, NPY_LONGDOUBLE);

// complex types
PYBIND11_NUMPY_SCALAR_IMPL(std::complex<float>, NPY_CFLOAT);
PYBIND11_NUMPY_SCALAR_IMPL(std::complex<double>, NPY_CDOUBLE);
PYBIND11_NUMPY_SCALAR_IMPL(std::complex<long double>, NPY_CLONGDOUBLE);

#undef PYBIND11_NUMPY_SCALAR_IMPL

// This table normalizes typenums by mapping NPY_INT_, NPY_LONG, ... to NPY_INT32_, NPY_INT64, ...
// This is needed to correctly handle situations where multiple typenums map to the same type,
// e.g. NPY_LONG_ may be equivalent to NPY_INT_ or NPY_LONGLONG_ despite having a different
// typenum. The normalized typenum should always match the values used in npy_format_descriptor.
// If you change this code, please review `enum constants` above.
static constexpr int normalized_dtype_num[npy_api::NPY_VOID_ + 1] = {
    // NPY_BOOL_ =>
    npy_api::NPY_BOOL_,
    // NPY_BYTE_ =>
    npy_api::NPY_BYTE_,
    // NPY_UBYTE_ =>
    npy_api::NPY_UBYTE_,
    // NPY_SHORT_ =>
    npy_api::NPY_INT16_,
    // NPY_USHORT_ =>
    npy_api::NPY_UINT16_,
    // NPY_INT_ =>
    sizeof(int) == sizeof(std::int16_t)   ? npy_api::NPY_INT16_
    : sizeof(int) == sizeof(std::int32_t) ? npy_api::NPY_INT32_
    : sizeof(int) == sizeof(std::int64_t) ? npy_api::NPY_INT64_
                                          : npy_api::NPY_INT_,
    // NPY_UINT_ =>
    sizeof(unsigned int) == sizeof(std::uint16_t)   ? npy_api::NPY_UINT16_
    : sizeof(unsigned int) == sizeof(std::uint32_t) ? npy_api::NPY_UINT32_
    : sizeof(unsigned int) == sizeof(std::uint64_t) ? npy_api::NPY_UINT64_
                                                    : npy_api::NPY_UINT_,
    // NPY_LONG_ =>
    sizeof(long) == sizeof(std::int16_t)   ? npy_api::NPY_INT16_
    : sizeof(long) == sizeof(std::int32_t) ? npy_api::NPY_INT32_
    : sizeof(long) == sizeof(std::int64_t) ? npy_api::NPY_INT64_
                                           : npy_api::NPY_LONG_,
    // NPY_ULONG_ =>
    sizeof(unsigned long) == sizeof(std::uint16_t)   ? npy_api::NPY_UINT16_
    : sizeof(unsigned long) == sizeof(std::uint32_t) ? npy_api::NPY_UINT32_
    : sizeof(unsigned long) == sizeof(std::uint64_t) ? npy_api::NPY_UINT64_
                                                     : npy_api::NPY_ULONG_,
    // NPY_LONGLONG_ =>
    sizeof(long long) == sizeof(std::int16_t)   ? npy_api::NPY_INT16_
    : sizeof(long long) == sizeof(std::int32_t) ? npy_api::NPY_INT32_
    : sizeof(long long) == sizeof(std::int64_t) ? npy_api::NPY_INT64_
                                                : npy_api::NPY_LONGLONG_,
    // NPY_ULONGLONG_ =>
    sizeof(unsigned long long) == sizeof(std::uint16_t)   ? npy_api::NPY_UINT16_
    : sizeof(unsigned long long) == sizeof(std::uint32_t) ? npy_api::NPY_UINT32_
    : sizeof(unsigned long long) == sizeof(std::uint64_t) ? npy_api::NPY_UINT64_
                                                          : npy_api::NPY_ULONGLONG_,
    // NPY_FLOAT_ =>
    npy_api::NPY_FLOAT_,
    // NPY_DOUBLE_ =>
    npy_api::NPY_DOUBLE_,
    // NPY_LONGDOUBLE_ =>
    npy_api::NPY_LONGDOUBLE_,
    // NPY_CFLOAT_ =>
    npy_api::NPY_CFLOAT_,
    // NPY_CDOUBLE_ =>
    npy_api::NPY_CDOUBLE_,
    // NPY_CLONGDOUBLE_ =>
    npy_api::NPY_CLONGDOUBLE_,
    // NPY_OBJECT_ =>
    npy_api::NPY_OBJECT_,
    // NPY_STRING_ =>
    npy_api::NPY_STRING_,
    // NPY_UNICODE_ =>
    npy_api::NPY_UNICODE_,
    // NPY_VOID_ =>
    npy_api::NPY_VOID_,
};

inline PyArray_Proxy *array_proxy(void *ptr) { return reinterpret_cast<PyArray_Proxy *>(ptr); }

inline const PyArray_Proxy *array_proxy(const void *ptr) {
    return reinterpret_cast<const PyArray_Proxy *>(ptr);
}

inline PyArrayDescr_Proxy *array_descriptor_proxy(PyObject *ptr) {
    return reinterpret_cast<PyArrayDescr_Proxy *>(ptr);
}

inline const PyArrayDescr_Proxy *array_descriptor_proxy(const PyObject *ptr) {
    return reinterpret_cast<const PyArrayDescr_Proxy *>(ptr);
}

inline const PyArrayDescr1_Proxy *array_descriptor1_proxy(const PyObject *ptr) {
    return reinterpret_cast<const PyArrayDescr1_Proxy *>(ptr);
}

inline const PyArrayDescr2_Proxy *array_descriptor2_proxy(const PyObject *ptr) {
    return reinterpret_cast<const PyArrayDescr2_Proxy *>(ptr);
}

inline bool check_flags(const void *ptr, int flag) {
    return (flag == (array_proxy(ptr)->flags & flag));
}

template <typename T>
struct is_std_array : std::false_type {};
template <typename T, size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
struct array_info_scalar {
    using type = T;
    static constexpr bool is_array = false;
    static constexpr bool is_empty = false;
    static constexpr auto extents = const_name("");
    static void append_extents(list & /* shape */) {}
};
// Computes underlying type and a comma-separated list of extents for array
// types (any mix of std::array and built-in arrays). An array of char is
// treated as scalar because it gets special handling.
template <typename T>
struct array_info : array_info_scalar<T> {};
template <typename T, size_t N>
struct array_info<std::array<T, N>> {
    using type = typename array_info<T>::type;
    static constexpr bool is_array = true;
    static constexpr bool is_empty = (N == 0) || array_info<T>::is_empty;
    static constexpr size_t extent = N;

    // appends the extents to shape
    static void append_extents(list &shape) {
        shape.append(N);
        array_info<T>::append_extents(shape);
    }

    static constexpr auto extents = const_name<array_info<T>::is_array>(
        ::pybind11::detail::concat(const_name<N>(), array_info<T>::extents), const_name<N>());
};
// For numpy we have special handling for arrays of characters, so we don't include
// the size in the array extents.
template <size_t N>
struct array_info<char[N]> : array_info_scalar<char[N]> {};
template <size_t N>
struct array_info<std::array<char, N>> : array_info_scalar<std::array<char, N>> {};
template <typename T, size_t N>
struct array_info<T[N]> : array_info<std::array<T, N>> {};
template <typename T>
using remove_all_extents_t = typename array_info<T>::type;

template <typename T>
using is_pod_struct
    = all_of<std::is_standard_layout<T>, // since we're accessing directly in memory
                                         // we need a standard layout type
#if defined(__GLIBCXX__)                                                                          \
    && (__GLIBCXX__ < 20150422 || __GLIBCXX__ == 20150426 || __GLIBCXX__ == 20150623              \
        || __GLIBCXX__ == 20150626 || __GLIBCXX__ == 20160803)
             // libstdc++ < 5 (including versions 4.8.5, 4.9.3 and 4.9.4 which were released after
             // 5) don't implement is_trivially_copyable, so approximate it
             std::is_trivially_destructible<T>,
             satisfies_any_of<T, std::has_trivial_copy_constructor, std::has_trivial_copy_assign>,
#else
             std::is_trivially_copyable<T>,
#endif
             satisfies_none_of<T,
                               std::is_reference,
                               std::is_array,
                               is_std_array,
                               std::is_arithmetic,
                               is_complex,
                               std::is_enum>>;

// Replacement for std::is_pod (deprecated in C++20)
template <typename T>
using is_pod = all_of<std::is_standard_layout<T>, std::is_trivial<T>>;

template <ssize_t Dim = 0, typename Strides>
ssize_t byte_offset_unsafe(const Strides &) {
    return 0;
}
template <ssize_t Dim = 0, typename Strides, typename... Ix>
ssize_t byte_offset_unsafe(const Strides &strides, ssize_t i, Ix... index) {
    return i * strides[Dim] + byte_offset_unsafe<Dim + 1>(strides, index...);
}

/**
 * Proxy class providing unsafe, unchecked const access to array data.  This is constructed through
 * the `unchecked<T, N>()` method of `array` or the `unchecked<N>()` method of `array_t<T>`. `Dims`
 * will be -1 for dimensions determined at runtime.
 */
template <typename T, ssize_t Dims>
class unchecked_reference {
protected:
    static constexpr bool Dynamic = Dims < 0;
    const unsigned char *data_;
    // Storing the shape & strides in local variables (i.e. these arrays) allows the compiler to
    // make large performance gains on big, nested loops, but requires compile-time dimensions
    conditional_t<Dynamic, const ssize_t *, std::array<ssize_t, (size_t) Dims>> shape_, strides_;
    const ssize_t dims_;

    friend class pybind11::array;
    // Constructor for compile-time dimensions:
    template <bool Dyn = Dynamic>
    unchecked_reference(const void *data,
                        const ssize_t *shape,
                        const ssize_t *strides,
                        enable_if_t<!Dyn, ssize_t>)
        : data_{reinterpret_cast<const unsigned char *>(data)}, dims_{Dims} {
        for (size_t i = 0; i < (size_t) dims_; i++) {
            shape_[i] = shape[i];
            strides_[i] = strides[i];
        }
    }
    // Constructor for runtime dimensions:
    template <bool Dyn = Dynamic>
    unchecked_reference(const void *data,
                        const ssize_t *shape,
                        const ssize_t *strides,
                        enable_if_t<Dyn, ssize_t> dims)
        : data_{reinterpret_cast<const unsigned char *>(data)}, shape_{shape}, strides_{strides},
          dims_{dims} {}

public:
    /**
     * Unchecked const reference access to data at the given indices.  For a compile-time known
     * number of dimensions, this requires the correct number of arguments; for run-time
     * dimensionality, this is not checked (and so is up to the caller to use safely).
     */
    template <typename... Ix>
    const T &operator()(Ix... index) const {
        static_assert(ssize_t{sizeof...(Ix)} == Dims || Dynamic,
                      "Invalid number of indices for unchecked array reference");
        return *reinterpret_cast<const T *>(data_
                                            + byte_offset_unsafe(strides_, ssize_t(index)...));
    }
    /**
     * Unchecked const reference access to data; this operator only participates if the reference
     * is to a 1-dimensional array.  When present, this is exactly equivalent to `obj(index)`.
     */
    template <ssize_t D = Dims, typename = enable_if_t<D == 1 || Dynamic>>
    const T &operator[](ssize_t index) const {
        return operator()(index);
    }

    /// Pointer access to the data at the given indices.
    template <typename... Ix>
    const T *data(Ix... ix) const {
        return &operator()(ssize_t(ix)...);
    }

    /// Returns the item size, i.e. sizeof(T)
    constexpr static ssize_t itemsize() { return sizeof(T); }

    /// Returns the shape (i.e. size) of dimension `dim`
    ssize_t shape(ssize_t dim) const { return shape_[(size_t) dim]; }

    /// Returns the number of dimensions of the array
    ssize_t ndim() const { return dims_; }

    /// Returns the total number of elements in the referenced array, i.e. the product of the
    /// shapes
    template <bool Dyn = Dynamic>
    enable_if_t<!Dyn, ssize_t> size() const {
        return std::accumulate(
            shape_.begin(), shape_.end(), (ssize_t) 1, std::multiplies<ssize_t>());
    }
    template <bool Dyn = Dynamic>
    enable_if_t<Dyn, ssize_t> size() const {
        return std::accumulate(shape_, shape_ + ndim(), (ssize_t) 1, std::multiplies<ssize_t>());
    }

    /// Returns the total number of bytes used by the referenced data.  Note that the actual span
    /// in memory may be larger if the referenced array has non-contiguous strides (e.g. for a
    /// slice).
    ssize_t nbytes() const { return size() * itemsize(); }
};

template <typename T, ssize_t Dims>
class unchecked_mutable_reference : public unchecked_reference<T, Dims> {
    friend class pybind11::array;
    using ConstBase = unchecked_reference<T, Dims>;
    using ConstBase::ConstBase;
    using ConstBase::Dynamic;

public:
    // Bring in const-qualified versions from base class
    using ConstBase::operator();
    using ConstBase::operator[];

    /// Mutable, unchecked access to data at the given indices.
    template <typename... Ix>
    T &operator()(Ix... index) {
        static_assert(ssize_t{sizeof...(Ix)} == Dims || Dynamic,
                      "Invalid number of indices for unchecked array reference");
        return const_cast<T &>(ConstBase::operator()(index...));
    }
    /**
     * Mutable, unchecked access data at the given index; this operator only participates if the
     * reference is to a 1-dimensional array (or has runtime dimensions).  When present, this is
     * exactly equivalent to `obj(index)`.
     */
    template <ssize_t D = Dims, typename = enable_if_t<D == 1 || Dynamic>>
    T &operator[](ssize_t index) {
        return operator()(index);
    }

    /// Mutable pointer access to the data at the given indices.
    template <typename... Ix>
    T *mutable_data(Ix... ix) {
        return &operator()(ssize_t(ix)...);
    }
};

template <typename T, ssize_t Dim>
struct type_caster<unchecked_reference<T, Dim>> {
    static_assert(Dim == 0 && Dim > 0 /* always fail */,
                  "unchecked array proxy object is not castable");
};
template <typename T, ssize_t Dim>
struct type_caster<unchecked_mutable_reference<T, Dim>>
    : type_caster<unchecked_reference<T, Dim>> {};

template <typename T>
struct type_caster<numpy_scalar<T>> {
    using value_type = T;
    using type_info = numpy_scalar_info<T>;

    PYBIND11_TYPE_CASTER(numpy_scalar<T>, type_info::name);

    static handle &target_type() {
        static handle tp = npy_api::get().PyArray_TypeObjectFromType_(type_info::typenum);
        return tp;
    }

    static handle &target_dtype() {
        static handle tp = npy_api::get().PyArray_DescrFromType_(type_info::typenum);
        return tp;
    }

    bool load(handle src, bool) {
        if (isinstance(src, target_type())) {
            npy_api::get().PyArray_ScalarAsCtype_(src.ptr(), &value.value);
            return true;
        }
        return false;
    }

    static handle cast(numpy_scalar<T> src, return_value_policy, handle) {
        return npy_api::get().PyArray_Scalar_(&src.value, target_dtype().ptr(), nullptr);
    }
};

PYBIND11_NAMESPACE_END(detail)

template <typename T>
struct numpy_scalar {
    using value_type = T;

    value_type value;

    numpy_scalar() = default;
    explicit numpy_scalar(value_type value) : value(value) {}

    explicit operator value_type() const { return value; }
    numpy_scalar &operator=(value_type value) {
        this->value = value;
        return *this;
    }

    friend bool operator==(const numpy_scalar &a, const numpy_scalar &b) {
        return a.value == b.value;
    }

    friend bool operator!=(const numpy_scalar &a, const numpy_scalar &b) { return !(a == b); }
};

template <typename T>
numpy_scalar<T> make_scalar(T value) {
    return numpy_scalar<T>(value);
}

class dtype : public object {
public:
    PYBIND11_OBJECT_DEFAULT(dtype, object, detail::npy_api::get().PyArrayDescr_Check_)

    explicit dtype(const buffer_info &info) {
        dtype descr(_dtype_from_pep3118()(pybind11::str(info.format)));
        // If info.itemsize == 0, use the value calculated from the format string
        m_ptr = descr.strip_padding(info.itemsize != 0 ? info.itemsize : descr.itemsize())
                    .release()
                    .ptr();
    }

    explicit dtype(const pybind11::str &format) : dtype(from_args(format)) {}

    explicit dtype(const std::string &format) : dtype(pybind11::str(format)) {}

    explicit dtype(const char *format) : dtype(pybind11::str(format)) {}

    dtype(list names, list formats, list offsets, ssize_t itemsize) {
        dict args;
        args["names"] = std::move(names);
        args["formats"] = std::move(formats);
        args["offsets"] = std::move(offsets);
        args["itemsize"] = pybind11::int_(itemsize);
        m_ptr = from_args(args).release().ptr();
    }

    /// Return dtype for the given typenum (one of the NPY_TYPES).
    /// https://numpy.org/devdocs/reference/c-api/array.html#c.PyArray_DescrFromType
    explicit dtype(int typenum)
        : object(detail::npy_api::get().PyArray_DescrFromType_(typenum), stolen_t{}) {
        if (m_ptr == nullptr) {
            throw error_already_set();
        }
    }

    /// This is essentially the same as calling numpy.dtype(args) in Python.
    static dtype from_args(const object &args) {
        PyObject *ptr = nullptr;
        if ((detail::npy_api::get().PyArray_DescrConverter_(args.ptr(), &ptr) == 0) || !ptr) {
            throw error_already_set();
        }
        return reinterpret_steal<dtype>(ptr);
    }

    /// Return dtype associated with a C++ type.
    template <typename T>
    static dtype of() {
        return detail::npy_format_descriptor<typename std::remove_cv<T>::type>::dtype();
    }

    /// Return the type number associated with a C++ type.
    /// This is the constexpr equivalent of `dtype::of<T>().num()`.
    template <typename T>
    static constexpr int num_of() {
        return detail::npy_format_descriptor<typename std::remove_cv<T>::type>::value;
    }

    /// Size of the data type in bytes.
    ssize_t itemsize() const {
        if (detail::npy_api::get().PyArray_RUNTIME_VERSION_ < 0x12) {
            return detail::array_descriptor1_proxy(m_ptr)->elsize;
        }
        return detail::array_descriptor2_proxy(m_ptr)->elsize;
    }

    /// Returns true for structured data types.
    bool has_fields() const {
        if (detail::npy_api::get().PyArray_RUNTIME_VERSION_ < 0x12) {
            return detail::array_descriptor1_proxy(m_ptr)->names != nullptr;
        }
        const auto *proxy = detail::array_descriptor2_proxy(m_ptr);
        if (proxy->type_num < 0 || proxy->type_num >= 2056) {
            return false;
        }
        return proxy->names != nullptr;
    }

    /// Single-character code for dtype's kind.
    /// For example, floating point types are 'f' and integral types are 'i'.
    char kind() const { return detail::array_descriptor_proxy(m_ptr)->kind; }

    /// Single-character for dtype's type.
    /// For example, ``float`` is 'f', ``double`` 'd', ``int`` 'i', and ``long`` 'l'.
    char char_() const {
        // Note: The signature, `dtype::char_` follows the naming of NumPy's
        // public Python API (i.e., ``dtype.char``), rather than its internal
        // C API (``PyArray_Descr::type``).
        return detail::array_descriptor_proxy(m_ptr)->type;
    }

    /// Type number of dtype. Note that different values may be returned for equivalent types,
    /// e.g. even though ``long`` may be equivalent to ``int`` or ``long long``, they still have
    /// different type numbers. Consider using `normalized_num` to avoid this.
    int num() const {
        // Note: The signature, `dtype::num` follows the naming of NumPy's public
        // Python API (i.e., ``dtype.num``), rather than its internal
        // C API (``PyArray_Descr::type_num``).
        return detail::array_descriptor_proxy(m_ptr)->type_num;
    }

    /// Type number of dtype, normalized to match the return value of `num_of` for equivalent
    /// types. This function can be used to write switch statements that correctly handle
    /// equivalent types with different type numbers.
    int normalized_num() const {
        int value = num();
        if (value >= 0 && value <= detail::npy_api::NPY_VOID_) {
            return detail::normalized_dtype_num[value];
        }
        return value;
    }

    /// Single character for byteorder
    char byteorder() const { return detail::array_descriptor_proxy(m_ptr)->byteorder; }

    /// Alignment of the data type
    ssize_t alignment() const {
        if (detail::npy_api::get().PyArray_RUNTIME_VERSION_ < 0x12) {
            return detail::array_descriptor1_proxy(m_ptr)->alignment;
        }
        return detail::array_descriptor2_proxy(m_ptr)->alignment;
    }

    /// Flags for the array descriptor
    std::uint64_t flags() const {
        if (detail::npy_api::get().PyArray_RUNTIME_VERSION_ < 0x12) {
            return (unsigned char) detail::array_descriptor1_proxy(m_ptr)->flags;
        }
        return detail::array_descriptor2_proxy(m_ptr)->flags;
    }

private:
    static object &_dtype_from_pep3118() {
        PYBIND11_CONSTINIT static gil_safe_call_once_and_store<object> storage;
        return storage
            .call_once_and_store_result([]() {
                return detail::import_numpy_core_submodule("_internal")
                    .attr("_dtype_from_pep3118");
            })
            .get_stored();
    }

    dtype strip_padding(ssize_t itemsize) {
        // Recursively strip all void fields with empty names that are generated for
        // padding fields (as of NumPy v1.11).
        if (!has_fields()) {
            return *this;
        }

        struct field_descr {
            pybind11::str name;
            object format;
            pybind11::int_ offset;
            field_descr(pybind11::str &&name, object &&format, pybind11::int_ &&offset)
                : name{std::move(name)}, format{std::move(format)}, offset{std::move(offset)} {};
        };
        auto field_dict = attr("fields").cast<dict>();
        std::vector<field_descr> field_descriptors;
        field_descriptors.reserve(field_dict.size());

        for (auto field : field_dict.attr("items")()) {
            auto spec = field.cast<tuple>();
            auto name = spec[0].cast<pybind11::str>();
            auto spec_fo = spec[1].cast<tuple>();
            auto format = spec_fo[0].cast<dtype>();
            auto offset = spec_fo[1].cast<pybind11::int_>();
            if ((len(name) == 0u) && format.kind() == 'V') {
                continue;
            }
            field_descriptors.emplace_back(
                std::move(name), format.strip_padding(format.itemsize()), std::move(offset));
        }

        std::sort(field_descriptors.begin(),
                  field_descriptors.end(),
                  [](const field_descr &a, const field_descr &b) {
                      return a.offset.cast<int>() < b.offset.cast<int>();
                  });

        list names, formats, offsets;
        for (auto &descr : field_descriptors) {
            names.append(std::move(descr.name));
            formats.append(std::move(descr.format));
            offsets.append(std::move(descr.offset));
        }
        return dtype(std::move(names), std::move(formats), std::move(offsets), itemsize);
    }
};

class array : public buffer {
public:
    PYBIND11_OBJECT_CVT(array, buffer, detail::npy_api::get().PyArray_Check_, raw_array)

    enum {
        c_style = detail::npy_api::NPY_ARRAY_C_CONTIGUOUS_,
        f_style = detail::npy_api::NPY_ARRAY_F_CONTIGUOUS_,
        forcecast = detail::npy_api::NPY_ARRAY_FORCECAST_
    };

    array() : array(0, static_cast<const double *>(nullptr)) {}

    using ShapeContainer = detail::any_container<ssize_t>;
    using StridesContainer = detail::any_container<ssize_t>;

    // Constructs an array taking shape/strides from arbitrary container types
    array(const pybind11::dtype &dt,
          ShapeContainer shape,
          StridesContainer strides,
          const void *ptr = nullptr,
          handle base = handle()) {

        if (strides->empty()) {
            *strides = detail::c_strides(*shape, dt.itemsize());
        }

        auto ndim = shape->size();
        if (ndim != strides->size()) {
            pybind11_fail("NumPy: shape ndim doesn't match strides ndim");
        }
        auto descr = dt;

        int flags = 0;
        if (base && ptr) {
            if (isinstance<array>(base)) {
                /* Copy flags from base (except ownership bit) */
                flags = reinterpret_borrow<array>(base).flags()
                        & ~detail::npy_api::NPY_ARRAY_OWNDATA_;
            } else {
                /* Writable by default, easy to downgrade later on if needed */
                flags = detail::npy_api::NPY_ARRAY_WRITEABLE_;
            }
        }

        auto &api = detail::npy_api::get();
        auto tmp = reinterpret_steal<object>(api.PyArray_NewFromDescr_(
            api.PyArray_Type_,
            descr.release().ptr(),
            (int) ndim,
            // Use reinterpret_cast for PyPy on Windows (remove if fixed, checked on 7.3.1)
            reinterpret_cast<Py_intptr_t *>(shape->data()),
            reinterpret_cast<Py_intptr_t *>(strides->data()),
            const_cast<void *>(ptr),
            flags,
            nullptr));
        if (!tmp) {
            throw error_already_set();
        }
        if (ptr) {
            if (base) {
                api.PyArray_SetBaseObject_(tmp.ptr(), base.inc_ref().ptr());
            } else {
                tmp = reinterpret_steal<object>(
                    api.PyArray_NewCopy_(tmp.ptr(), -1 /* any order */));
            }
        }
        m_ptr = tmp.release().ptr();
    }

    array(const pybind11::dtype &dt,
          ShapeContainer shape,
          const void *ptr = nullptr,
          handle base = handle())
        : array(dt, std::move(shape), {}, ptr, base) {}

    template <typename T,
              typename
              = detail::enable_if_t<std::is_integral<T>::value && !std::is_same<bool, T>::value>>
    array(const pybind11::dtype &dt, T count, const void *ptr = nullptr, handle base = handle())
        : array(dt, {{count}}, ptr, base) {}

    template <typename T>
    array(ShapeContainer shape, StridesContainer strides, const T *ptr, handle base = handle())
        : array(pybind11::dtype::of<T>(),
                std::move(shape),
                std::move(strides),
                reinterpret_cast<const void *>(ptr),
                base) {}

    template <typename T>
    array(ShapeContainer shape, const T *ptr, handle base = handle())
        : array(std::move(shape), {}, ptr, base) {}

    template <typename T>
    explicit array(ssize_t count, const T *ptr, handle base = handle())
        : array({count}, {}, ptr, base) {}

    explicit array(const buffer_info &info, handle base = handle())
        : array(pybind11::dtype(info), info.shape, info.strides, info.ptr, base) {}

    /// Array descriptor (dtype)
    pybind11::dtype dtype() const {
        return reinterpret_borrow<pybind11::dtype>(detail::array_proxy(m_ptr)->descr);
    }

    /// Total number of elements
    ssize_t size() const {
        return std::accumulate(shape(), shape() + ndim(), (ssize_t) 1, std::multiplies<ssize_t>());
    }

    /// Byte size of a single element
    ssize_t itemsize() const { return dtype().itemsize(); }

    /// Total number of bytes
    ssize_t nbytes() const { return size() * itemsize(); }

    /// Number of dimensions
    ssize_t ndim() const { return detail::array_proxy(m_ptr)->nd; }

    /// Base object
    object base() const { return reinterpret_borrow<object>(detail::array_proxy(m_ptr)->base); }

    /// Dimensions of the array
    const ssize_t *shape() const { return detail::array_proxy(m_ptr)->dimensions; }

#ifdef PYBIND11_HAS_SPAN
    /// Dimensions of the array as a span
    std::span<const ssize_t, std::dynamic_extent> shape_span() const {
        return std::span(shape(), static_cast<std::size_t>(ndim()));
    }
#endif

    /// Dimension along a given axis
    ssize_t shape(ssize_t dim) const {
        if (dim >= ndim()) {
            fail_dim_check(dim, "invalid axis");
        }
        return shape()[dim];
    }

    /// Strides of the array
    const ssize_t *strides() const { return detail::array_proxy(m_ptr)->strides; }

#ifdef PYBIND11_HAS_SPAN
    /// Strides of the array as a span
    std::span<const ssize_t, std::dynamic_extent> strides_span() const {
        return std::span(strides(), static_cast<std::size_t>(ndim()));
    }
#endif

    /// Stride along a given axis
    ssize_t strides(ssize_t dim) const {
        if (dim >= ndim()) {
            fail_dim_check(dim, "invalid axis");
        }
        return strides()[dim];
    }

    /// Return the NumPy array flags
    int flags() const { return detail::array_proxy(m_ptr)->flags; }

    /// If set, the array is writeable (otherwise the buffer is read-only)
    bool writeable() const {
        return detail::check_flags(m_ptr, detail::npy_api::NPY_ARRAY_WRITEABLE_);
    }

    /// If set, the array owns the data (will be freed when the array is deleted)
    bool owndata() const {
        return detail::check_flags(m_ptr, detail::npy_api::NPY_ARRAY_OWNDATA_);
    }

    /// Pointer to the contained data. If index is not provided, points to the
    /// beginning of the buffer. May throw if the index would lead to out of bounds access.
    template <typename... Ix>
    const void *data(Ix... index) const {
        return static_cast<const void *>(detail::array_proxy(m_ptr)->data + offset_at(index...));
    }

    /// Mutable pointer to the contained data. If index is not provided, points to the
    /// beginning of the buffer. May throw if the index would lead to out of bounds access.
    /// May throw if the array is not writeable.
    template <typename... Ix>
    void *mutable_data(Ix... index) {
        check_writeable();
        return static_cast<void *>(detail::array_proxy(m_ptr)->data + offset_at(index...));
    }

    /// Byte offset from beginning of the array to a given index (full or partial).
    /// May throw if the index would lead to out of bounds access.
    template <typename... Ix>
    ssize_t offset_at(Ix... index) const {
        if ((ssize_t) sizeof...(index) > ndim()) {
            fail_dim_check(sizeof...(index), "too many indices for an array");
        }
        return byte_offset(ssize_t(index)...);
    }

    ssize_t offset_at() const { return 0; }

    /// Item count from beginning of the array to a given index (full or partial).
    /// May throw if the index would lead to out of bounds access.
    template <typename... Ix>
    ssize_t index_at(Ix... index) const {
        return offset_at(index...) / itemsize();
    }

    /**
     * Returns a proxy object that provides access to the array's data without bounds or
     * dimensionality checking.  Will throw if the array is missing the `writeable` flag.  Use with
     * care: the array must not be destroyed or reshaped for the duration of the returned object,
     * and the caller must take care not to access invalid dimensions or dimension indices.
     */
    template <typename T, ssize_t Dims = -1>
    detail::unchecked_mutable_reference<T, Dims> mutable_unchecked() & {
        if (Dims >= 0 && ndim() != Dims) {
            throw std::domain_error("array has incorrect number of dimensions: "
                                    + std::to_string(ndim()) + "; expected "
                                    + std::to_string(Dims));
        }
        return detail::unchecked_mutable_reference<T, Dims>(
            mutable_data(), shape(), strides(), ndim());
    }

    /**
     * Returns a proxy object that provides const access to the array's data without bounds or
     * dimensionality checking.  Unlike `mutable_unchecked()`, this does not require that the
     * underlying array have the `writable` flag.  Use with care: the array must not be destroyed
     * or reshaped for the duration of the returned object, and the caller must take care not to
     * access invalid dimensions or dimension indices.
     */
    template <typename T, ssize_t Dims = -1>
    detail::unchecked_reference<T, Dims> unchecked() const & {
        if (Dims >= 0 && ndim() != Dims) {
            throw std::domain_error("array has incorrect number of dimensions: "
                                    + std::to_string(ndim()) + "; expected "
                                    + std::to_string(Dims));
        }
        return detail::unchecked_reference<T, Dims>(data(), shape(), strides(), ndim());
    }

    /// Return a new view with all of the dimensions of length 1 removed
    array squeeze() {
        auto &api = detail::npy_api::get();
        return reinterpret_steal<array>(api.PyArray_Squeeze_(m_ptr));
    }

    /// Resize array to given shape
    /// If refcheck is true and more that one reference exist to this array
    /// then resize will succeed only if it makes a reshape, i.e. original size doesn't change
    void resize(ShapeContainer new_shape, bool refcheck = true) {
        detail::npy_api::PyArray_Dims d
            = {// Use reinterpret_cast for PyPy on Windows (remove if fixed, checked on 7.3.1)
               reinterpret_cast<Py_intptr_t *>(new_shape->data()),
               int(new_shape->size())};
        // try to resize, set ordering param to -1 cause it's not used anyway
        auto new_array = reinterpret_steal<object>(
            detail::npy_api::get().PyArray_Resize_(m_ptr, &d, int(refcheck), -1));
        if (!new_array) {
            throw error_already_set();
        }
        if (isinstance<array>(new_array)) {
            *this = std::move(new_array);
        }
    }

    /// Optional `order` parameter omitted, to be added as needed.
    array reshape(ShapeContainer new_shape) {
        detail::npy_api::PyArray_Dims d
            = {reinterpret_cast<Py_intptr_t *>(new_shape->data()), int(new_shape->size())};
        auto new_array
            = reinterpret_steal<array>(detail::npy_api::get().PyArray_Newshape_(m_ptr, &d, 0));
        if (!new_array) {
            throw error_already_set();
        }
        return new_array;
    }

    /// Create a view of an array in a different data type.
    /// This function may fundamentally reinterpret the data in the array.
    /// It is the responsibility of the caller to ensure that this is safe.
    /// Only supports the `dtype` argument, the `type` argument is omitted,
    /// to be added as needed.
    array view(const std::string &dtype) {
        auto &api = detail::npy_api::get();
        auto new_view = reinterpret_steal<array>(api.PyArray_View_(
            m_ptr, dtype::from_args(pybind11::str(dtype)).release().ptr(), nullptr));
        if (!new_view) {
            throw error_already_set();
        }
        return new_view;
    }

    /// Ensure that the argument is a NumPy array
    /// In case of an error, nullptr is returned and the Python error is cleared.
    static array ensure(handle h, int ExtraFlags = 0) {
        auto result = reinterpret_steal<array>(raw_array(h.ptr(), ExtraFlags));
        if (!result) {
            PyErr_Clear();
        }
        return result;
    }

protected:
    template <typename, typename>
    friend struct detail::npy_format_descriptor;

    void fail_dim_check(ssize_t dim, const std::string &msg) const {
        throw index_error(msg + ": " + std::to_string(dim) + " (ndim = " + std::to_string(ndim())
                          + ')');
    }

    template <typename... Ix>
    ssize_t byte_offset(Ix... index) const {
        check_dimensions(index...);
        return detail::byte_offset_unsafe(strides(), ssize_t(index)...);
    }

    void check_writeable() const {
        if (!writeable()) {
            throw std::domain_error("array is not writeable");
        }
    }

    template <typename... Ix>
    void check_dimensions(Ix... index) const {
        check_dimensions_impl(ssize_t(0), shape(), ssize_t(index)...);
    }

    void check_dimensions_impl(ssize_t, const ssize_t *) const {}

    template <typename... Ix>
    void check_dimensions_impl(ssize_t axis, const ssize_t *shape, ssize_t i, Ix... index) const {
        if (i >= *shape) {
            throw index_error(std::string("index ") + std::to_string(i)
                              + " is out of bounds for axis " + std::to_string(axis)
                              + " with size " + std::to_string(*shape));
        }
        check_dimensions_impl(axis + 1, shape + 1, index...);
    }

    /// Create array from any object -- always returns a new reference
    static PyObject *raw_array(PyObject *ptr, int ExtraFlags = 0) {
        if (ptr == nullptr) {
            set_error(PyExc_ValueError, "cannot create a pybind11::array from a nullptr");
            return nullptr;
        }
        return detail::npy_api::get().PyArray_FromAny_(
            ptr, nullptr, 0, 0, detail::npy_api::NPY_ARRAY_ENSUREARRAY_ | ExtraFlags, nullptr);
    }
};

template <typename T, int ExtraFlags = array::forcecast>
class array_t : public array {
private:
    struct private_ctor {};
    // Delegating constructor needed when both moving and accessing in the same constructor
    array_t(private_ctor,
            ShapeContainer &&shape,
            StridesContainer &&strides,
            const T *ptr,
            handle base)
        : array(std::move(shape), std::move(strides), ptr, base) {}

public:
    static_assert(!detail::array_info<T>::is_array, "Array types cannot be used with array_t");

    using value_type = T;

    array_t() : array(0, static_cast<const T *>(nullptr)) {}
    array_t(handle h, borrowed_t) : array(h, borrowed_t{}) {}
    array_t(handle h, stolen_t) : array(h, stolen_t{}) {}

    PYBIND11_DEPRECATED("Use array_t<T>::ensure() instead")
    array_t(handle h, bool is_borrowed) : array(raw_array_t(h.ptr()), stolen_t{}) {
        if (!m_ptr) {
            PyErr_Clear();
        }
        if (!is_borrowed) {
            Py_XDECREF(h.ptr());
        }
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    array_t(const object &o) : array(raw_array_t(o.ptr()), stolen_t{}) {
        if (!m_ptr) {
            throw error_already_set();
        }
    }

    explicit array_t(const buffer_info &info, handle base = handle()) : array(info, base) {}

    array_t(ShapeContainer shape,
            StridesContainer strides,
            const T *ptr = nullptr,
            handle base = handle())
        : array(std::move(shape), std::move(strides), ptr, base) {}

    explicit array_t(ShapeContainer shape, const T *ptr = nullptr, handle base = handle())
        : array_t(private_ctor{},
                  std::move(shape),
                  (ExtraFlags & f_style) != 0 ? detail::f_strides(*shape, itemsize())
                                              : detail::c_strides(*shape, itemsize()),
                  ptr,
                  base) {}

    explicit array_t(ssize_t count, const T *ptr = nullptr, handle base = handle())
        : array({count}, {}, ptr, base) {}

    constexpr ssize_t itemsize() const { return sizeof(T); }

    template <typename... Ix>
    ssize_t index_at(Ix... index) const {
        return offset_at(index...) / itemsize();
    }

    template <typename... Ix>
    const T *data(Ix... index) const {
        return static_cast<const T *>(array::data(index...));
    }

    template <typename... Ix>
    T *mutable_data(Ix... index) {
        return static_cast<T *>(array::mutable_data(index...));
    }

    // Reference to element at a given index
    template <typename... Ix>
    const T &at(Ix... index) const {
        if ((ssize_t) sizeof...(index) != ndim()) {
            fail_dim_check(sizeof...(index), "index dimension mismatch");
        }
        return *(static_cast<const T *>(array::data())
                 + byte_offset(ssize_t(index)...) / itemsize());
    }

    // Mutable reference to element at a given index
    template <typename... Ix>
    T &mutable_at(Ix... index) {
        if ((ssize_t) sizeof...(index) != ndim()) {
            fail_dim_check(sizeof...(index), "index dimension mismatch");
        }
        return *(static_cast<T *>(array::mutable_data())
                 + byte_offset(ssize_t(index)...) / itemsize());
    }

    /**
     * Returns a proxy object that provides access to the array's data without bounds or
     * dimensionality checking.  Will throw if the array is missing the `writeable` flag.  Use with
     * care: the array must not be destroyed or reshaped for the duration of the returned object,
     * and the caller must take care not to access invalid dimensions or dimension indices.
     */
    template <ssize_t Dims = -1>
    detail::unchecked_mutable_reference<T, Dims> mutable_unchecked() & {
        return array::mutable_unchecked<T, Dims>();
    }

    /**
     * Returns a proxy object that provides const access to the array's data without bounds or
     * dimensionality checking.  Unlike `mutable_unchecked()`, this does not require that the
     * underlying array have the `writable` flag.  Use with care: the array must not be destroyed
     * or reshaped for the duration of the returned object, and the caller must take care not to
     * access invalid dimensions or dimension indices.
     */
    template <ssize_t Dims = -1>
    detail::unchecked_reference<T, Dims> unchecked() const & {
        return array::unchecked<T, Dims>();
    }

    /// Ensure that the argument is a NumPy array of the correct dtype (and if not, try to convert
    /// it).  In case of an error, nullptr is returned and the Python error is cleared.
    static array_t ensure(handle h) {
        auto result = reinterpret_steal<array_t>(raw_array_t(h.ptr()));
        if (!result) {
            PyErr_Clear();
        }
        return result;
    }

    static bool check_(handle h) {
        const auto &api = detail::npy_api::get();
        return api.PyArray_Check_(h.ptr())
               && api.PyArray_EquivTypes_(detail::array_proxy(h.ptr())->descr,
                                          dtype::of<T>().ptr())
               && detail::check_flags(h.ptr(), ExtraFlags & (array::c_style | array::f_style));
    }

protected:
    /// Create array from any object -- always returns a new reference
    static PyObject *raw_array_t(PyObject *ptr) {
        if (ptr == nullptr) {
            set_error(PyExc_ValueError, "cannot create a pybind11::array_t from a nullptr");
            return nullptr;
        }
        return detail::npy_api::get().PyArray_FromAny_(ptr,
                                                       dtype::of<T>().release().ptr(),
                                                       0,
                                                       0,
                                                       detail::npy_api::NPY_ARRAY_ENSUREARRAY_
                                                           | ExtraFlags,
                                                       nullptr);
    }
};

template <typename T>
struct format_descriptor<T, detail::enable_if_t<detail::is_pod_struct<T>::value>> {
    static std::string format() {
        return detail::npy_format_descriptor<typename std::remove_cv<T>::type>::format();
    }
};

template <size_t N>
struct format_descriptor<char[N]> {
    static std::string format() { return std::to_string(N) + 's'; }
};
template <size_t N>
struct format_descriptor<std::array<char, N>> {
    static std::string format() { return std::to_string(N) + 's'; }
};

template <typename T>
struct format_descriptor<T, detail::enable_if_t<std::is_enum<T>::value>> {
    static std::string format() {
        return format_descriptor<
            typename std::remove_cv<typename std::underlying_type<T>::type>::type>::format();
    }
};

template <typename T>
struct format_descriptor<T, detail::enable_if_t<detail::array_info<T>::is_array>> {
    static std::string format() {
        using namespace detail;
        static constexpr auto extents = const_name("(") + array_info<T>::extents + const_name(")");
        return extents.text + format_descriptor<remove_all_extents_t<T>>::format();
    }
};

PYBIND11_NAMESPACE_BEGIN(detail)
template <typename T, int ExtraFlags>
struct pyobject_caster<array_t<T, ExtraFlags>> {
    using type = array_t<T, ExtraFlags>;

    bool load(handle src, bool convert) {
        if (!convert && !type::check_(src)) {
            return false;
        }
        value = type::ensure(src);
        return static_cast<bool>(value);
    }

    static handle cast(const handle &src, return_value_policy /* policy */, handle /* parent */) {
        return src.inc_ref();
    }
    PYBIND11_TYPE_CASTER(type, handle_type_name<type>::name);
};

template <typename T>
struct compare_buffer_info<T, detail::enable_if_t<detail::is_pod_struct<T>::value>> {
    static bool compare(const buffer_info &b) {
        return npy_api::get().PyArray_EquivTypes_(dtype::of<T>().ptr(), dtype(b).ptr());
    }
};

template <typename T>
struct npy_format_descriptor<
    T,
    enable_if_t<satisfies_any_of<T, std::is_arithmetic, is_complex>::value>>
    : npy_format_descriptor_name<T> {
private:
    // NB: the order here must match the one in common.h
    constexpr static const int values[15] = {npy_api::NPY_BOOL_,
                                             npy_api::NPY_BYTE_,
                                             npy_api::NPY_UBYTE_,
                                             npy_api::NPY_INT16_,
                                             npy_api::NPY_UINT16_,
                                             npy_api::NPY_INT32_,
                                             npy_api::NPY_UINT32_,
                                             npy_api::NPY_INT64_,
                                             npy_api::NPY_UINT64_,
                                             npy_api::NPY_FLOAT_,
                                             npy_api::NPY_DOUBLE_,
                                             npy_api::NPY_LONGDOUBLE_,
                                             npy_api::NPY_CFLOAT_,
                                             npy_api::NPY_CDOUBLE_,
                                             npy_api::NPY_CLONGDOUBLE_};

public:
    static constexpr int value = values[detail::is_fmt_numeric<T>::index];

    static pybind11::dtype dtype() { return pybind11::dtype(/*typenum*/ value); }
};

template <typename T>
struct npy_format_descriptor<
    T,
    enable_if_t<is_same_ignoring_cvref<T, PyObject *>::value
                || ((std::is_same<T, handle>::value || std::is_same<T, object>::value)
                    && sizeof(T) == sizeof(PyObject *))>> {
    static constexpr auto name = const_name("numpy.object_");

    static constexpr int value = npy_api::NPY_OBJECT_;

    static pybind11::dtype dtype() { return pybind11::dtype(/*typenum*/ value); }
};

#define PYBIND11_DECL_CHAR_FMT                                                                    \
    static constexpr auto name = const_name("S") + const_name<N>();                               \
    static pybind11::dtype dtype() {                                                              \
        return pybind11::dtype(std::string("S") + std::to_string(N));                             \
    }
template <size_t N>
struct npy_format_descriptor<char[N]> {
    PYBIND11_DECL_CHAR_FMT
};
template <size_t N>
struct npy_format_descriptor<std::array<char, N>> {
    PYBIND11_DECL_CHAR_FMT
};
#undef PYBIND11_DECL_CHAR_FMT

template <typename T>
struct npy_format_descriptor<T, enable_if_t<array_info<T>::is_array>> {
private:
    using base_descr = npy_format_descriptor<typename array_info<T>::type>;

public:
    static_assert(!array_info<T>::is_empty, "Zero-sized arrays are not supported");

    static constexpr auto name
        = const_name("(") + array_info<T>::extents + const_name(")") + base_descr::name;
    static pybind11::dtype dtype() {
        list shape;
        array_info<T>::append_extents(shape);
        return pybind11::dtype::from_args(
            pybind11::make_tuple(base_descr::dtype(), std::move(shape)));
    }
};

template <typename T>
struct npy_format_descriptor<T, enable_if_t<std::is_enum<T>::value>> {
private:
    using base_descr = npy_format_descriptor<typename std::underlying_type<T>::type>;

public:
    static constexpr auto name = base_descr::name;
    static pybind11::dtype dtype() { return base_descr::dtype(); }
};

struct field_descriptor {
    const char *name;
    ssize_t offset;
    ssize_t size;
    std::string format;
    dtype descr;
};

PYBIND11_NOINLINE void register_structured_dtype(any_container<field_descriptor> fields,
                                                 const std::type_info &tinfo,
                                                 ssize_t itemsize,
                                                 bool (*direct_converter)(PyObject *, void *&)) {

    auto &numpy_internals = get_numpy_internals();
    if (numpy_internals.get_type_info(tinfo, false)) {
        pybind11_fail("NumPy: dtype is already registered");
    }

    // Use ordered fields because order matters as of NumPy 1.14:
    // https://docs.scipy.org/doc/numpy/release.html#multiple-field-indexing-assignment-of-structured-arrays
    std::vector<field_descriptor> ordered_fields(std::move(fields));
    std::sort(
        ordered_fields.begin(),
        ordered_fields.end(),
        [](const field_descriptor &a, const field_descriptor &b) { return a.offset < b.offset; });

    list names, formats, offsets;
    for (auto &field : ordered_fields) {
        if (!field.descr) {
            pybind11_fail(std::string("NumPy: unsupported field dtype: `") + field.name + "` @ "
                          + tinfo.name());
        }
        names.append(pybind11::str(field.name));
        formats.append(field.descr);
        offsets.append(pybind11::int_(field.offset));
    }
    auto *dtype_ptr
        = pybind11::dtype(std::move(names), std::move(formats), std::move(offsets), itemsize)
              .release()
              .ptr();

    // There is an existing bug in NumPy (as of v1.11): trailing bytes are
    // not encoded explicitly into the format string. This will supposedly
    // get fixed in v1.12; for further details, see these:
    // - https://github.com/numpy/numpy/issues/7797
    // - https://github.com/numpy/numpy/pull/7798
    // Because of this, we won't use numpy's logic to generate buffer format
    // strings and will just do it ourselves.
    ssize_t offset = 0;
    std::ostringstream oss;
    // mark the structure as unaligned with '^', because numpy and C++ don't
    // always agree about alignment (particularly for complex), and we're
    // explicitly listing all our padding. This depends on none of the fields
    // overriding the endianness. Putting the ^ in front of individual fields
    // isn't guaranteed to work due to https://github.com/numpy/numpy/issues/9049
    oss << "^T{";
    for (auto &field : ordered_fields) {
        if (field.offset > offset) {
            oss << (field.offset - offset) << 'x';
        }
        oss << field.format << ':' << field.name << ':';
        offset = field.offset + field.size;
    }
    if (itemsize > offset) {
        oss << (itemsize - offset) << 'x';
    }
    oss << '}';
    auto format_str = oss.str();

    // Smoke test: verify that NumPy properly parses our buffer format string
    auto &api = npy_api::get();
    auto arr = array(buffer_info(nullptr, itemsize, format_str, 1));
    if (!api.PyArray_EquivTypes_(dtype_ptr, arr.dtype().ptr())) {
        pybind11_fail("NumPy: invalid buffer descriptor!");
    }

    auto tindex = std::type_index(tinfo);
    numpy_internals.registered_dtypes[tindex] = {dtype_ptr, std::move(format_str)};
    with_internals([tindex, &direct_converter](internals &internals) {
        internals.direct_conversions[tindex].push_back(direct_converter);
    });
}

template <typename T, typename SFINAE>
struct npy_format_descriptor {
    static_assert(is_pod_struct<T>::value,
                  "Attempt to use a non-POD or unimplemented POD type as a numpy dtype");

    static constexpr auto name = make_caster<T>::name;

    static pybind11::dtype dtype() { return reinterpret_borrow<pybind11::dtype>(dtype_ptr()); }

    static std::string format() {
        static auto format_str = get_numpy_internals().get_type_info<T>(true)->format_str;
        return format_str;
    }

    static void register_dtype(any_container<field_descriptor> fields) {
        register_structured_dtype(std::move(fields),
                                  typeid(typename std::remove_cv<T>::type),
                                  sizeof(T),
                                  &direct_converter);
    }

private:
    static PyObject *dtype_ptr() {
        static PyObject *ptr = get_numpy_internals().get_type_info<T>(true)->dtype_ptr;
        return ptr;
    }

    static bool direct_converter(PyObject *obj, void *&value) {
        auto &api = npy_api::get();
        if (!PyObject_TypeCheck(obj, api.PyVoidArrType_Type_)) {
            return false;
        }
        if (auto descr = reinterpret_steal<object>(api.PyArray_DescrFromScalar_(obj))) {
            if (api.PyArray_EquivTypes_(dtype_ptr(), descr.ptr())) {
                value = ((PyVoidScalarObject_Proxy *) obj)->obval;
                return true;
            }
        }
        return false;
    }
};

#ifdef __CLION_IDE__ // replace heavy macro with dummy code for the IDE (doesn't affect code)
#    define PYBIND11_NUMPY_DTYPE(Type, ...) ((void) 0)
#    define PYBIND11_NUMPY_DTYPE_EX(Type, ...) ((void) 0)
#else

#    define PYBIND11_FIELD_DESCRIPTOR_EX(T, Field, Name)                                          \
        ::pybind11::detail::field_descriptor {                                                    \
            Name, offsetof(T, Field), sizeof(decltype(std::declval<T>().Field)),                  \
                ::pybind11::format_descriptor<decltype(std::declval<T>().Field)>::format(),       \
                ::pybind11::detail::npy_format_descriptor<                                        \
                    decltype(std::declval<T>().Field)>::dtype()                                   \
        }

// Extract name, offset and format descriptor for a struct field
#    define PYBIND11_FIELD_DESCRIPTOR(T, Field) PYBIND11_FIELD_DESCRIPTOR_EX(T, Field, #Field)

// The main idea of this macro is borrowed from https://github.com/swansontec/map-macro
// (C) William Swanson, Paul Fultz
#    define PYBIND11_EVAL0(...) __VA_ARGS__
#    define PYBIND11_EVAL1(...) PYBIND11_EVAL0(PYBIND11_EVAL0(PYBIND11_EVAL0(__VA_ARGS__)))
#    define PYBIND11_EVAL2(...) PYBIND11_EVAL1(PYBIND11_EVAL1(PYBIND11_EVAL1(__VA_ARGS__)))
#    define PYBIND11_EVAL3(...) PYBIND11_EVAL2(PYBIND11_EVAL2(PYBIND11_EVAL2(__VA_ARGS__)))
#    define PYBIND11_EVAL4(...) PYBIND11_EVAL3(PYBIND11_EVAL3(PYBIND11_EVAL3(__VA_ARGS__)))
#    define PYBIND11_EVAL(...) PYBIND11_EVAL4(PYBIND11_EVAL4(PYBIND11_EVAL4(__VA_ARGS__)))
#    define PYBIND11_MAP_END(...)
#    define PYBIND11_MAP_OUT
#    define PYBIND11_MAP_COMMA ,
#    define PYBIND11_MAP_GET_END() 0, PYBIND11_MAP_END
#    define PYBIND11_MAP_NEXT0(test, next, ...) next PYBIND11_MAP_OUT
#    define PYBIND11_MAP_NEXT1(test, next) PYBIND11_MAP_NEXT0(test, next, 0)
#    define PYBIND11_MAP_NEXT(test, next) PYBIND11_MAP_NEXT1(PYBIND11_MAP_GET_END test, next)
#    if defined(_MSC_VER)                                                                         \
        && !defined(__clang__) // MSVC is not as eager to expand macros, hence this workaround
#        define PYBIND11_MAP_LIST_NEXT1(test, next)                                               \
            PYBIND11_EVAL0(PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0))
#    else
#        define PYBIND11_MAP_LIST_NEXT1(test, next)                                               \
            PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0)
#    endif
#    define PYBIND11_MAP_LIST_NEXT(test, next)                                                    \
        PYBIND11_MAP_LIST_NEXT1(PYBIND11_MAP_GET_END test, next)
#    define PYBIND11_MAP_LIST0(f, t, x, peek, ...)                                                \
        f(t, x) PYBIND11_MAP_LIST_NEXT(peek, PYBIND11_MAP_LIST1)(f, t, peek, __VA_ARGS__)
#    define PYBIND11_MAP_LIST1(f, t, x, peek, ...)                                                \
        f(t, x) PYBIND11_MAP_LIST_NEXT(peek, PYBIND11_MAP_LIST0)(f, t, peek, __VA_ARGS__)
// PYBIND11_MAP_LIST(f, t, a1, a2, ...) expands to f(t, a1), f(t, a2), ...
#    define PYBIND11_MAP_LIST(f, t, ...)                                                          \
        PYBIND11_EVAL(PYBIND11_MAP_LIST1(f, t, __VA_ARGS__, (), 0))

#    define PYBIND11_NUMPY_DTYPE(Type, ...)                                                       \
        ::pybind11::detail::npy_format_descriptor<Type>::register_dtype(                          \
            ::std::vector<::pybind11::detail::field_descriptor>{                                  \
                PYBIND11_MAP_LIST(PYBIND11_FIELD_DESCRIPTOR, Type, __VA_ARGS__)})

#    if defined(_MSC_VER) && !defined(__clang__)
#        define PYBIND11_MAP2_LIST_NEXT1(test, next)                                              \
            PYBIND11_EVAL0(PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0))
#    else
#        define PYBIND11_MAP2_LIST_NEXT1(test, next)                                              \
            PYBIND11_MAP_NEXT0(test, PYBIND11_MAP_COMMA next, 0)
#    endif
#    define PYBIND11_MAP2_LIST_NEXT(test, next)                                                   \
        PYBIND11_MAP2_LIST_NEXT1(PYBIND11_MAP_GET_END test, next)
#    define PYBIND11_MAP2_LIST0(f, t, x1, x2, peek, ...)                                          \
        f(t, x1, x2) PYBIND11_MAP2_LIST_NEXT(peek, PYBIND11_MAP2_LIST1)(f, t, peek, __VA_ARGS__)
#    define PYBIND11_MAP2_LIST1(f, t, x1, x2, peek, ...)                                          \
        f(t, x1, x2) PYBIND11_MAP2_LIST_NEXT(peek, PYBIND11_MAP2_LIST0)(f, t, peek, __VA_ARGS__)
// PYBIND11_MAP2_LIST(f, t, a1, a2, ...) expands to f(t, a1, a2), f(t, a3, a4), ...
#    define PYBIND11_MAP2_LIST(f, t, ...)                                                         \
        PYBIND11_EVAL(PYBIND11_MAP2_LIST1(f, t, __VA_ARGS__, (), 0))

#    define PYBIND11_NUMPY_DTYPE_EX(Type, ...)                                                    \
        ::pybind11::detail::npy_format_descriptor<Type>::register_dtype(                          \
            ::std::vector<::pybind11::detail::field_descriptor>{                                  \
                PYBIND11_MAP2_LIST(PYBIND11_FIELD_DESCRIPTOR_EX, Type, __VA_ARGS__)})

#endif // __CLION_IDE__

class common_iterator {
public:
    using container_type = std::vector<ssize_t>;
    using value_type = container_type::value_type;
    using size_type = container_type::size_type;

    common_iterator() = default;

    common_iterator(void *ptr, const container_type &strides, const container_type &shape)
        : p_ptr(reinterpret_cast<char *>(ptr)), m_strides(strides.size()) {
        m_strides.back() = static_cast<value_type>(strides.back());
        for (size_type i = m_strides.size() - 1; i != 0; --i) {
            size_type j = i - 1;
            auto s = static_cast<value_type>(shape[i]);
            m_strides[j] = strides[j] + m_strides[i] - strides[i] * s;
        }
    }

    void increment(size_type dim) { p_ptr += m_strides[dim]; }

    void *data() const { return p_ptr; }

private:
    char *p_ptr{nullptr};
    container_type m_strides;
};

template <size_t N>
class multi_array_iterator {
public:
    using container_type = std::vector<ssize_t>;

    multi_array_iterator(const std::array<buffer_info, N> &buffers, const container_type &shape)
        : m_shape(shape.size()), m_index(shape.size(), 0), m_common_iterator() {

        // Manual copy to avoid conversion warning if using std::copy
        for (size_t i = 0; i < shape.size(); ++i) {
            m_shape[i] = shape[i];
        }

        container_type strides(shape.size());
        for (size_t i = 0; i < N; ++i) {
            init_common_iterator(buffers[i], shape, m_common_iterator[i], strides);
        }
    }

    multi_array_iterator &operator++() {
        for (size_t j = m_index.size(); j != 0; --j) {
            size_t i = j - 1;
            if (++m_index[i] != m_shape[i]) {
                increment_common_iterator(i);
                break;
            }
            m_index[i] = 0;
        }
        return *this;
    }

    template <size_t K, class T = void>
    T *data() const {
        return reinterpret_cast<T *>(m_common_iterator[K].data());
    }

private:
    using common_iter = common_iterator;

    void init_common_iterator(const buffer_info &buffer,
                              const container_type &shape,
                              common_iter &iterator,
                              container_type &strides) {
        auto buffer_shape_iter = buffer.shape.rbegin();
        auto buffer_strides_iter = buffer.strides.rbegin();
        auto shape_iter = shape.rbegin();
        auto strides_iter = strides.rbegin();

        while (buffer_shape_iter != buffer.shape.rend()) {
            if (*shape_iter == *buffer_shape_iter) {
                *strides_iter = *buffer_strides_iter;
            } else {
                *strides_iter = 0;
            }

            ++buffer_shape_iter;
            ++buffer_strides_iter;
            ++shape_iter;
            ++strides_iter;
        }

        std::fill(strides_iter, strides.rend(), 0);
        iterator = common_iter(buffer.ptr, strides, shape);
    }

    void increment_common_iterator(size_t dim) {
        for (auto &iter : m_common_iterator) {
            iter.increment(dim);
        }
    }

    container_type m_shape;
    container_type m_index;
    std::array<common_iter, N> m_common_iterator;
};

enum class broadcast_trivial { non_trivial, c_trivial, f_trivial };

// Populates the shape and number of dimensions for the set of buffers.  Returns a
// broadcast_trivial enum value indicating whether the broadcast is "trivial"--that is, has each
// buffer being either a singleton or a full-size, C-contiguous (`c_trivial`) or Fortran-contiguous
// (`f_trivial`) storage buffer; returns `non_trivial` otherwise.
template <size_t N>
broadcast_trivial
broadcast(const std::array<buffer_info, N> &buffers, ssize_t &ndim, std::vector<ssize_t> &shape) {
    ndim = std::accumulate(
        buffers.begin(), buffers.end(), ssize_t(0), [](ssize_t res, const buffer_info &buf) {
            return std::max(res, buf.ndim);
        });

    shape.clear();
    shape.resize((size_t) ndim, 1);

    // Figure out the output size, and make sure all input arrays conform (i.e. are either size 1
    // or the full size).
    for (size_t i = 0; i < N; ++i) {
        auto res_iter = shape.rbegin();
        auto end = buffers[i].shape.rend();
        for (auto shape_iter = buffers[i].shape.rbegin(); shape_iter != end;
             ++shape_iter, ++res_iter) {
            const auto &dim_size_in = *shape_iter;
            auto &dim_size_out = *res_iter;

            // Each input dimension can either be 1 or `n`, but `n` values must match across
            // buffers
            if (dim_size_out == 1) {
                dim_size_out = dim_size_in;
            } else if (dim_size_in != 1 && dim_size_in != dim_size_out) {
                pybind11_fail("pybind11::vectorize: incompatible size/dimension of inputs!");
            }
        }
    }

    bool trivial_broadcast_c = true;
    bool trivial_broadcast_f = true;
    for (size_t i = 0; i < N && (trivial_broadcast_c || trivial_broadcast_f); ++i) {
        if (buffers[i].size == 1) {
            continue;
        }

        // Require the same number of dimensions:
        if (buffers[i].ndim != ndim) {
            return broadcast_trivial::non_trivial;
        }

        // Require all dimensions be full-size:
        if (!std::equal(buffers[i].shape.cbegin(), buffers[i].shape.cend(), shape.cbegin())) {
            return broadcast_trivial::non_trivial;
        }

        // Check for C contiguity (but only if previous inputs were also C contiguous)
        if (trivial_broadcast_c) {
            ssize_t expect_stride = buffers[i].itemsize;
            auto end = buffers[i].shape.crend();
            for (auto shape_iter = buffers[i].shape.crbegin(),
                      stride_iter = buffers[i].strides.crbegin();
                 trivial_broadcast_c && shape_iter != end;
                 ++shape_iter, ++stride_iter) {
                if (expect_stride == *stride_iter) {
                    expect_stride *= *shape_iter;
                } else {
                    trivial_broadcast_c = false;
                }
            }
        }

        // Check for Fortran contiguity (if previous inputs were also F contiguous)
        if (trivial_broadcast_f) {
            ssize_t expect_stride = buffers[i].itemsize;
            auto end = buffers[i].shape.cend();
            for (auto shape_iter = buffers[i].shape.cbegin(),
                      stride_iter = buffers[i].strides.cbegin();
                 trivial_broadcast_f && shape_iter != end;
                 ++shape_iter, ++stride_iter) {
                if (expect_stride == *stride_iter) {
                    expect_stride *= *shape_iter;
                } else {
                    trivial_broadcast_f = false;
                }
            }
        }
    }

    return trivial_broadcast_c   ? broadcast_trivial::c_trivial
           : trivial_broadcast_f ? broadcast_trivial::f_trivial
                                 : broadcast_trivial::non_trivial;
}

template <typename T>
struct vectorize_arg {
    static_assert(!std::is_rvalue_reference<T>::value,
                  "Functions with rvalue reference arguments cannot be vectorized");
    // The wrapped function gets called with this type:
    using call_type = remove_reference_t<T>;
    // Is this a vectorized argument?
    static constexpr bool vectorize
        = satisfies_any_of<call_type, std::is_arithmetic, is_complex, is_pod>::value
          && satisfies_none_of<call_type,
                               std::is_pointer,
                               std::is_array,
                               is_std_array,
                               std::is_enum>::value
          && (!std::is_reference<T>::value
              || (std::is_lvalue_reference<T>::value && std::is_const<call_type>::value));
    // Accept this type: an array for vectorized types, otherwise the type as-is:
    using type = conditional_t<vectorize, array_t<remove_cv_t<call_type>, array::forcecast>, T>;
};

// py::vectorize when a return type is present
template <typename Func, typename Return, typename... Args>
struct vectorize_returned_array {
    using Type = array_t<Return>;

    static Type create(broadcast_trivial trivial, const std::vector<ssize_t> &shape) {
        if (trivial == broadcast_trivial::f_trivial) {
            return array_t<Return, array::f_style>(shape);
        }
        return array_t<Return>(shape);
    }

    static Return *mutable_data(Type &array) { return array.mutable_data(); }

    static Return call(Func &f, Args &...args) { return f(args...); }

    static void call(Return *out, size_t i, Func &f, Args &...args) { out[i] = f(args...); }
};

// py::vectorize when a return type is not present
template <typename Func, typename... Args>
struct vectorize_returned_array<Func, void, Args...> {
    using Type = none;

    static Type create(broadcast_trivial, const std::vector<ssize_t> &) { return none(); }

    static void *mutable_data(Type &) { return nullptr; }

    static detail::void_type call(Func &f, Args &...args) {
        f(args...);
        return {};
    }

    static void call(void *, size_t, Func &f, Args &...args) { f(args...); }
};

template <typename Func, typename Return, typename... Args>
struct vectorize_helper {

// NVCC for some reason breaks if NVectorized is private
#ifdef __CUDACC__
public:
#else
private:
#endif

    static constexpr size_t N = sizeof...(Args);
    static constexpr size_t NVectorized = constexpr_sum(vectorize_arg<Args>::vectorize...);
    static_assert(
        NVectorized >= 1,
        "pybind11::vectorize(...) requires a function with at least one vectorizable argument");

public:
    template <typename T,
              // SFINAE to prevent shadowing the copy constructor.
              typename = detail::enable_if_t<
                  !std::is_same<vectorize_helper, typename std::decay<T>::type>::value>>
    explicit vectorize_helper(T &&f) : f(std::forward<T>(f)) {}

    object operator()(typename vectorize_arg<Args>::type... args) {
        return run(args...,
                   make_index_sequence<N>(),
                   select_indices<vectorize_arg<Args>::vectorize...>(),
                   make_index_sequence<NVectorized>());
    }

private:
    remove_reference_t<Func> f;

    // Internal compiler error in MSVC 19.16.27025.1 (Visual Studio 2017 15.9.4), when compiling
    // with "/permissive-" flag when arg_call_types is manually inlined.
    using arg_call_types = std::tuple<typename vectorize_arg<Args>::call_type...>;
    template <size_t Index>
    using param_n_t = typename std::tuple_element<Index, arg_call_types>::type;

    using returned_array = vectorize_returned_array<Func, Return, Args...>;

    // Runs a vectorized function given arguments tuple and three index sequences:
    //     - Index is the full set of 0 ... (N-1) argument indices;
    //     - VIndex is the subset of argument indices with vectorized parameters, letting us access
    //       vectorized arguments (anything not in this sequence is passed through)
    //     - BIndex is a incremental sequence (beginning at 0) of the same size as VIndex, so that
    //       we can store vectorized buffer_infos in an array (argument VIndex has its buffer at
    //       index BIndex in the array).
    template <size_t... Index, size_t... VIndex, size_t... BIndex>
    object run(typename vectorize_arg<Args>::type &...args,
               index_sequence<Index...> i_seq,
               index_sequence<VIndex...> vi_seq,
               index_sequence<BIndex...> bi_seq) {

        // Pointers to values the function was called with; the vectorized ones set here will start
        // out as array_t<T> pointers, but they will be changed them to T pointers before we make
        // call the wrapped function.  Non-vectorized pointers are left as-is.
        std::array<void *, N> params{{reinterpret_cast<void *>(&args)...}};

        // The array of `buffer_info`s of vectorized arguments:
        std::array<buffer_info, NVectorized> buffers{
            {reinterpret_cast<array *>(params[VIndex])->request()...}};

        /* Determine dimensions parameters of output array */
        ssize_t nd = 0;
        std::vector<ssize_t> shape(0);
        auto trivial = broadcast(buffers, nd, shape);
        auto ndim = (size_t) nd;

        size_t size
            = std::accumulate(shape.begin(), shape.end(), (size_t) 1, std::multiplies<size_t>());

        // If all arguments are 0-dimension arrays (i.e. single values) return a plain value (i.e.
        // not wrapped in an array).
        if (size == 1 && ndim == 0) {
            PYBIND11_EXPAND_SIDE_EFFECTS(params[VIndex] = buffers[BIndex].ptr);
            return cast(
                returned_array::call(f, *reinterpret_cast<param_n_t<Index> *>(params[Index])...));
        }

        auto result = returned_array::create(trivial, shape);

        PYBIND11_WARNING_PUSH
#ifdef PYBIND11_DETECTED_CLANG_WITH_MISLEADING_CALL_STD_MOVE_EXPLICITLY_WARNING
        PYBIND11_WARNING_DISABLE_CLANG("-Wreturn-std-move")
#endif

        if (size == 0) {
            return result;
        }

        /* Call the function */
        auto *mutable_data = returned_array::mutable_data(result);
        if (trivial == broadcast_trivial::non_trivial) {
            apply_broadcast(buffers, params, mutable_data, size, shape, i_seq, vi_seq, bi_seq);
        } else {
            apply_trivial(buffers, params, mutable_data, size, i_seq, vi_seq, bi_seq);
        }

        return result;
        PYBIND11_WARNING_POP
    }

    template <size_t... Index, size_t... VIndex, size_t... BIndex>
    void apply_trivial(std::array<buffer_info, NVectorized> &buffers,
                       std::array<void *, N> &params,
                       Return *out,
                       size_t size,
                       index_sequence<Index...>,
                       index_sequence<VIndex...>,
                       index_sequence<BIndex...>) {

        // Initialize an array of mutable byte references and sizes with references set to the
        // appropriate pointer in `params`; as we iterate, we'll increment each pointer by its size
        // (except for singletons, which get an increment of 0).
        std::array<std::pair<unsigned char *&, const size_t>, NVectorized> vecparams{
            {std::pair<unsigned char *&, const size_t>(
                reinterpret_cast<unsigned char *&>(params[VIndex] = buffers[BIndex].ptr),
                buffers[BIndex].size == 1 ? 0 : sizeof(param_n_t<VIndex>))...}};

        for (size_t i = 0; i < size; ++i) {
            returned_array::call(
                out, i, f, *reinterpret_cast<param_n_t<Index> *>(params[Index])...);
            for (auto &x : vecparams) {
                x.first += x.second;
            }
        }
    }

    template <size_t... Index, size_t... VIndex, size_t... BIndex>
    void apply_broadcast(std::array<buffer_info, NVectorized> &buffers,
                         std::array<void *, N> &params,
                         Return *out,
                         size_t size,
                         const std::vector<ssize_t> &output_shape,
                         index_sequence<Index...>,
                         index_sequence<VIndex...>,
                         index_sequence<BIndex...>) {

        multi_array_iterator<NVectorized> input_iter(buffers, output_shape);

        for (size_t i = 0; i < size; ++i, ++input_iter) {
            PYBIND11_EXPAND_SIDE_EFFECTS((params[VIndex] = input_iter.template data<BIndex>()));
            returned_array::call(
                out, i, f, *reinterpret_cast<param_n_t<Index> *>(std::get<Index>(params))...);
        }
    }
};

template <typename Func, typename Return, typename... Args>
vectorize_helper<Func, Return, Args...> vectorize_extractor(const Func &f, Return (*)(Args...)) {
    return detail::vectorize_helper<Func, Return, Args...>(f);
}

template <typename T, int Flags>
struct handle_type_name<array_t<T, Flags>> {
    static constexpr auto name
        = io_name("typing.Annotated[numpy.typing.ArrayLike, ", "numpy.typing.NDArray[")
          + npy_format_descriptor<T>::name + const_name("]");
};

PYBIND11_NAMESPACE_END(detail)

// Vanilla pointer vectorizer:
template <typename Return, typename... Args>
detail::vectorize_helper<Return (*)(Args...), Return, Args...> vectorize(Return (*f)(Args...)) {
    return detail::vectorize_helper<Return (*)(Args...), Return, Args...>(f);
}

// lambda vectorizer:
template <typename Func, detail::enable_if_t<detail::is_lambda<Func>::value, int> = 0>
auto vectorize(Func &&f)
    -> decltype(detail::vectorize_extractor(std::forward<Func>(f),
                                            (detail::function_signature_t<Func> *) nullptr)) {
    return detail::vectorize_extractor(std::forward<Func>(f),
                                       (detail::function_signature_t<Func> *) nullptr);
}

// Vectorize a class method (non-const):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...)>())),
              Return,
              Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...)) {
    return Helper(std::mem_fn(f));
}

// Vectorize a class method (const):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) const>())),
              Return,
              const Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) const) {
    return Helper(std::mem_fn(f));
}

// Intentionally no &&/const&& overloads: vectorized method calls operate on the bound Python
// instance and should not consume/move-from self.
// Vectorize a class method (non-const, lvalue ref-qualified):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) &>())),
              Return,
              Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) &) {
    return Helper(std::mem_fn(f));
}

// Vectorize a class method (const, lvalue ref-qualified):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) const &>())),
              Return,
              const Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) const &) {
    return Helper(std::mem_fn(f));
}

#ifdef __cpp_noexcept_function_type
// Vectorize a class method (non-const, noexcept):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) noexcept>())),
              Return,
              Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) noexcept) {
    return Helper(std::mem_fn(f));
}

// Vectorize a class method (const, noexcept):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) const noexcept>())),
              Return,
              const Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) const noexcept) {
    return Helper(std::mem_fn(f));
}

// Vectorize a class method (non-const, lvalue ref-qualified, noexcept):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) & noexcept>())),
              Return,
              Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) & noexcept) {
    return Helper(std::mem_fn(f));
}

// Vectorize a class method (const, lvalue ref-qualified, noexcept):
template <typename Return,
          typename Class,
          typename... Args,
          typename Helper = detail::vectorize_helper<
              decltype(std::mem_fn(std::declval<Return (Class::*)(Args...) const & noexcept>())),
              Return,
              const Class *,
              Args...>>
Helper vectorize(Return (Class::*f)(Args...) const & noexcept) {
    return Helper(std::mem_fn(f));
}
#endif

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
