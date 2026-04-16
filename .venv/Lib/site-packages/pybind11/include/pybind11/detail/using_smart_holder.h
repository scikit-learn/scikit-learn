// Copyright (c) 2024 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "common.h"
#include "struct_smart_holder.h"

#include <type_traits>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

using pybind11::memory::smart_holder;

PYBIND11_NAMESPACE_BEGIN(detail)

template <typename H>
using is_smart_holder = std::is_same<H, smart_holder>;

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
