// Copyright (c) 2021 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "common.h"

#include <type_traits>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <typename To, typename From, typename SFINAE = void>
struct dynamic_raw_ptr_cast_is_possible : std::false_type {};

template <typename To, typename From>
struct dynamic_raw_ptr_cast_is_possible<
    To,
    From,
    detail::enable_if_t<!std::is_same<To, void>::value && std::is_polymorphic<From>::value>>
    : std::true_type {};

template <typename To,
          typename From,
          detail::enable_if_t<!dynamic_raw_ptr_cast_is_possible<To, From>::value, int> = 0>
To *dynamic_raw_ptr_cast_if_possible(From * /*ptr*/) {
    return nullptr;
}

template <typename To,
          typename From,
          detail::enable_if_t<dynamic_raw_ptr_cast_is_possible<To, From>::value, int> = 0>
To *dynamic_raw_ptr_cast_if_possible(From *ptr) {
    return dynamic_cast<To *>(ptr);
}

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
