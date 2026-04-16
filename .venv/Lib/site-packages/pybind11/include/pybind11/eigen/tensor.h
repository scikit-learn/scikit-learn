/*
    pybind11/eigen/tensor.h: Transparent conversion for Eigen tensors

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <pybind11/numpy.h>

#include "common.h"

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
static_assert(__GNUC__ > 5, "Eigen Tensor support in pybind11 requires GCC > 5.0");
#endif

// Disable warnings for Eigen
PYBIND11_WARNING_PUSH
PYBIND11_WARNING_DISABLE_MSVC(4554)
PYBIND11_WARNING_DISABLE_MSVC(4127)
#if defined(__MINGW32__)
PYBIND11_WARNING_DISABLE_GCC("-Wmaybe-uninitialized")
#endif

#include <unsupported/Eigen/CXX11/Tensor>

PYBIND11_WARNING_POP

static_assert(EIGEN_VERSION_AT_LEAST(3, 3, 0),
              "Eigen Tensor support in pybind11 requires Eigen >= 3.3.0");

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_WARNING_DISABLE_MSVC(4127)

PYBIND11_NAMESPACE_BEGIN(detail)

inline bool is_tensor_aligned(const void *data) {
    return (reinterpret_cast<std::size_t>(data) % EIGEN_DEFAULT_ALIGN_BYTES) == 0;
}

template <typename T>
constexpr int compute_array_flag_from_tensor() {
    static_assert((static_cast<int>(T::Layout) == static_cast<int>(Eigen::RowMajor))
                      || (static_cast<int>(T::Layout) == static_cast<int>(Eigen::ColMajor)),
                  "Layout must be row or column major");
    return (static_cast<int>(T::Layout) == static_cast<int>(Eigen::RowMajor)) ? array::c_style
                                                                              : array::f_style;
}

template <typename T>
struct eigen_tensor_helper {};

template <typename Scalar_, int NumIndices_, int Options_, typename IndexType>
struct eigen_tensor_helper<Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>> {
    using Type = Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>;
    using ValidType = void;

    static Eigen::DSizes<typename Type::Index, Type::NumIndices> get_shape(const Type &f) {
        return f.dimensions();
    }

    static constexpr bool
    is_correct_shape(const Eigen::DSizes<typename Type::Index, Type::NumIndices> & /*shape*/) {
        return true;
    }

    template <typename T>
    struct helper {};

    template <size_t... Is>
    struct helper<index_sequence<Is...>> {
        static constexpr auto value = ::pybind11::detail::concat(const_name(((void) Is, "?"))...);
    };

    static constexpr auto dimensions_descriptor
        = helper<decltype(make_index_sequence<Type::NumIndices>())>::value;

    template <typename... Args>
    static Type *alloc(Args &&...args) {
        return new Type(std::forward<Args>(args)...);
    }

    static void free(Type *tensor) { delete tensor; }
};

template <typename Scalar_, typename std::ptrdiff_t... Indices, int Options_, typename IndexType>
struct eigen_tensor_helper<
    Eigen::TensorFixedSize<Scalar_, Eigen::Sizes<Indices...>, Options_, IndexType>> {
    using Type = Eigen::TensorFixedSize<Scalar_, Eigen::Sizes<Indices...>, Options_, IndexType>;
    using ValidType = void;

    static constexpr Eigen::DSizes<typename Type::Index, Type::NumIndices>
    get_shape(const Type & /*f*/) {
        return get_shape();
    }

    static constexpr Eigen::DSizes<typename Type::Index, Type::NumIndices> get_shape() {
        return Eigen::DSizes<typename Type::Index, Type::NumIndices>(Indices...);
    }

    static bool
    is_correct_shape(const Eigen::DSizes<typename Type::Index, Type::NumIndices> &shape) {
        return get_shape() == shape;
    }

    static constexpr auto dimensions_descriptor
        = ::pybind11::detail::concat(const_name<Indices>()...);

    template <typename... Args>
    static Type *alloc(Args &&...args) {
        Eigen::aligned_allocator<Type> allocator;
        return ::new (allocator.allocate(1)) Type(std::forward<Args>(args)...);
    }

    static void free(Type *tensor) {
        Eigen::aligned_allocator<Type> allocator;
        tensor->~Type();
        allocator.deallocate(tensor, 1);
    }
};

template <typename Type, bool ShowDetails, bool NeedsWriteable = false>
struct get_tensor_descriptor {
    static constexpr auto details
        = const_name<NeedsWriteable>(", \"flags.writeable\"", "")
          + const_name<static_cast<int>(Type::Layout) == static_cast<int>(Eigen::RowMajor)>(
              ", \"flags.c_contiguous\"", ", \"flags.f_contiguous\"");
    static constexpr auto value
        = const_name("typing.Annotated[")
          + io_name("numpy.typing.ArrayLike, ", "numpy.typing.NDArray[")
          + npy_format_descriptor<typename Type::Scalar>::name + io_name("", "]")
          + const_name(", \"[") + eigen_tensor_helper<remove_cv_t<Type>>::dimensions_descriptor
          + const_name("]\"") + const_name<ShowDetails>(details, const_name("")) + const_name("]");
};

// When EIGEN_AVOID_STL_ARRAY is defined, Eigen::DSizes<T, 0> does not have the begin() member
// function. Falling back to a simple loop works around this issue.
//
// We need to disable the type-limits warning for the inner loop when size = 0.

PYBIND11_WARNING_PUSH
PYBIND11_WARNING_DISABLE_GCC("-Wtype-limits")

template <typename T, int size>
std::vector<T> convert_dsizes_to_vector(const Eigen::DSizes<T, size> &arr) {
    std::vector<T> result(size);

    for (size_t i = 0; i < size; i++) {
        result[i] = arr[i];
    }

    return result;
}

template <typename T, int size>
Eigen::DSizes<T, size> get_shape_for_array(const array &arr) {
    Eigen::DSizes<T, size> result;
    const T *shape = arr.shape();
    for (size_t i = 0; i < size; i++) {
        result[i] = shape[i];
    }

    return result;
}

PYBIND11_WARNING_POP

template <typename Type>
struct type_caster<Type, typename eigen_tensor_helper<Type>::ValidType> {
    static_assert(!std::is_pointer<typename Type::Scalar>::value,
                  PYBIND11_EIGEN_MESSAGE_POINTER_TYPES_ARE_NOT_SUPPORTED);
    using Helper = eigen_tensor_helper<Type>;
    static constexpr auto temp_name = get_tensor_descriptor<Type, false>::value;
    PYBIND11_TYPE_CASTER(Type, temp_name);

    bool load(handle src, bool convert) {
        if (!convert) {
            if (!isinstance<array>(src)) {
                return false;
            }
            array temp = array::ensure(src);
            if (!temp) {
                return false;
            }

            if (!temp.dtype().is(dtype::of<typename Type::Scalar>())) {
                return false;
            }
        }

        array_t<typename Type::Scalar, compute_array_flag_from_tensor<Type>()> arr(
            reinterpret_borrow<object>(src));

        if (arr.ndim() != Type::NumIndices) {
            return false;
        }
        auto shape = get_shape_for_array<typename Type::Index, Type::NumIndices>(arr);

        if (!Helper::is_correct_shape(shape)) {
            return false;
        }

#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
        auto data_pointer = arr.data();
#else
        // Handle Eigen bug
        auto data_pointer = const_cast<typename Type::Scalar *>(arr.data());
#endif

        if (is_tensor_aligned(arr.data())) {
            value = Eigen::TensorMap<const Type, Eigen::Aligned>(data_pointer, shape);
        } else {
            value = Eigen::TensorMap<const Type>(data_pointer, shape);
        }

        return true;
    }

    static handle cast(Type &&src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::reference
            || policy == return_value_policy::reference_internal) {
            pybind11_fail("Cannot use a reference return value policy for an rvalue");
        }
        return cast_impl(&src, return_value_policy::move, parent);
    }

    static handle cast(const Type &&src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::reference
            || policy == return_value_policy::reference_internal) {
            pybind11_fail("Cannot use a reference return value policy for an rvalue");
        }
        return cast_impl(&src, return_value_policy::move, parent);
    }

    static handle cast(Type &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic
            || policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::copy;
        }
        return cast_impl(&src, policy, parent);
    }

    static handle cast(const Type &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic
            || policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::copy;
        }
        return cast(&src, policy, parent);
    }

    static handle cast(Type *src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic) {
            policy = return_value_policy::take_ownership;
        } else if (policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::reference;
        }
        return cast_impl(src, policy, parent);
    }

    static handle cast(const Type *src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic) {
            policy = return_value_policy::take_ownership;
        } else if (policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::reference;
        }
        return cast_impl(src, policy, parent);
    }

    template <typename C>
    static handle cast_impl(C *src, return_value_policy policy, handle parent) {
        object parent_object;
        bool writeable = false;
        switch (policy) {
            case return_value_policy::move:
                if (std::is_const<C>::value) {
                    pybind11_fail("Cannot move from a constant reference");
                }

                src = Helper::alloc(std::move(*src));

                parent_object
                    = capsule(src, [](void *ptr) { Helper::free(reinterpret_cast<Type *>(ptr)); });
                writeable = true;
                break;

            case return_value_policy::take_ownership:
                if (std::is_const<C>::value) {
                    // This cast is ugly, and might be UB in some cases, but we don't have an
                    // alternative here as we must free that memory
                    Helper::free(const_cast<Type *>(src));
                    pybind11_fail("Cannot take ownership of a const reference");
                }

                parent_object
                    = capsule(src, [](void *ptr) { Helper::free(reinterpret_cast<Type *>(ptr)); });
                writeable = true;
                break;

            case return_value_policy::copy:
                writeable = true;
                break;

            case return_value_policy::reference:
                parent_object = none();
                writeable = !std::is_const<C>::value;
                break;

            case return_value_policy::reference_internal:
                // Default should do the right thing
                if (!parent) {
                    pybind11_fail("Cannot use reference internal when there is no parent");
                }
                parent_object = reinterpret_borrow<object>(parent);
                writeable = !std::is_const<C>::value;
                break;

            default:
                pybind11_fail("pybind11 bug in eigen.h, please file a bug report");
        }

        auto result = array_t<typename Type::Scalar, compute_array_flag_from_tensor<Type>()>(
            convert_dsizes_to_vector(Helper::get_shape(*src)), src->data(), parent_object);

        if (!writeable) {
            array_proxy(result.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;
        }

        return result.release();
    }
};

template <typename StoragePointerType,
          bool needs_writeable,
          enable_if_t<!needs_writeable, bool> = true>
StoragePointerType get_array_data_for_type(array &arr) {
#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
    return reinterpret_cast<StoragePointerType>(arr.data());
#else
    // Handle Eigen bug
    return reinterpret_cast<StoragePointerType>(const_cast<void *>(arr.data()));
#endif
}

template <typename StoragePointerType,
          bool needs_writeable,
          enable_if_t<needs_writeable, bool> = true>
StoragePointerType get_array_data_for_type(array &arr) {
    return reinterpret_cast<StoragePointerType>(arr.mutable_data());
}

template <typename T, typename = void>
struct get_storage_pointer_type;

template <typename MapType>
struct get_storage_pointer_type<MapType, void_t<typename MapType::StoragePointerType>> {
    using SPT = typename MapType::StoragePointerType;
};

template <typename MapType>
struct get_storage_pointer_type<MapType, void_t<typename MapType::PointerArgType>> {
    using SPT = typename MapType::PointerArgType;
};

template <typename Type, int Options>
struct type_caster<Eigen::TensorMap<Type, Options>,
                   typename eigen_tensor_helper<remove_cv_t<Type>>::ValidType> {
    static_assert(!std::is_pointer<typename Type::Scalar>::value,
                  PYBIND11_EIGEN_MESSAGE_POINTER_TYPES_ARE_NOT_SUPPORTED);
    using MapType = Eigen::TensorMap<Type, Options>;
    using Helper = eigen_tensor_helper<remove_cv_t<Type>>;

    bool load(handle src, bool /*convert*/) {
        // Note that we have a lot more checks here as we want to make sure to avoid copies
        if (!isinstance<array>(src)) {
            return false;
        }
        auto arr = reinterpret_borrow<array>(src);
        if ((arr.flags() & compute_array_flag_from_tensor<Type>()) == 0) {
            return false;
        }

        if (!arr.dtype().is(dtype::of<typename Type::Scalar>())) {
            return false;
        }

        if (arr.ndim() != Type::NumIndices) {
            return false;
        }

        constexpr bool is_aligned = (Options & Eigen::Aligned) != 0;

        if (is_aligned && !is_tensor_aligned(arr.data())) {
            return false;
        }

        auto shape = get_shape_for_array<typename Type::Index, Type::NumIndices>(arr);

        if (!Helper::is_correct_shape(shape)) {
            return false;
        }

        if (needs_writeable && !arr.writeable()) {
            return false;
        }

        auto result = get_array_data_for_type<typename get_storage_pointer_type<MapType>::SPT,
                                              needs_writeable>(arr);

        value.reset(new MapType(std::move(result), std::move(shape)));

        return true;
    }

    static handle cast(MapType &&src, return_value_policy policy, handle parent) {
        return cast_impl(&src, policy, parent);
    }

    static handle cast(const MapType &&src, return_value_policy policy, handle parent) {
        return cast_impl(&src, policy, parent);
    }

    static handle cast(MapType &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic
            || policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::copy;
        }
        return cast_impl(&src, policy, parent);
    }

    static handle cast(const MapType &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic
            || policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::copy;
        }
        return cast(&src, policy, parent);
    }

    static handle cast(MapType *src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic) {
            policy = return_value_policy::take_ownership;
        } else if (policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::reference;
        }
        return cast_impl(src, policy, parent);
    }

    static handle cast(const MapType *src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic) {
            policy = return_value_policy::take_ownership;
        } else if (policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::reference;
        }
        return cast_impl(src, policy, parent);
    }

    template <typename C>
    static handle cast_impl(C *src, return_value_policy policy, handle parent) {
        object parent_object;
        constexpr bool writeable = !std::is_const<C>::value;
        switch (policy) {
            case return_value_policy::reference:
                parent_object = none();
                break;

            case return_value_policy::reference_internal:
                // Default should do the right thing
                if (!parent) {
                    pybind11_fail("Cannot use reference internal when there is no parent");
                }
                parent_object = reinterpret_borrow<object>(parent);
                break;

            default:
                // move, take_ownership don't make any sense for a ref/map:
                pybind11_fail("Invalid return_value_policy for Eigen Map type, must be either "
                              "reference or reference_internal");
        }

        auto result = array_t<typename Type::Scalar, compute_array_flag_from_tensor<Type>()>(
            convert_dsizes_to_vector(Helper::get_shape(*src)),
            src->data(),
            std::move(parent_object));

        if (!writeable) {
            array_proxy(result.ptr())->flags &= ~detail::npy_api::NPY_ARRAY_WRITEABLE_;
        }

        return result.release();
    }

#if EIGEN_VERSION_AT_LEAST(3, 4, 0)

    static constexpr bool needs_writeable = !std::is_const<typename std::remove_pointer<
        typename get_storage_pointer_type<MapType>::SPT>::type>::value;
#else
    // Handle Eigen bug
    static constexpr bool needs_writeable = !std::is_const<Type>::value;
#endif

protected:
    // TODO: Move to std::optional once std::optional has more support
    std::unique_ptr<MapType> value;

public:
    // return_descr forces the use of NDArray instead of ArrayLike since refs can only reference
    // arrays
    static constexpr auto name
        = return_descr(get_tensor_descriptor<Type, true, needs_writeable>::value);
    explicit operator MapType *() { return value.get(); }
    explicit operator MapType &() { return *value; }
    explicit operator MapType &&() && { return std::move(*value); }

    template <typename T_>
    using cast_op_type = ::pybind11::detail::movable_cast_op_type<T_>;
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
