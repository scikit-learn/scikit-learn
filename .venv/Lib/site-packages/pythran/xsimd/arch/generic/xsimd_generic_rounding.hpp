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

#ifndef XSIMD_GENERIC_ROUNDING_HPP
#define XSIMD_GENERIC_ROUNDING_HPP

#include "./xsimd_generic_details.hpp"

namespace xsimd
{

    namespace kernel
    {

        using namespace types;

        // ceil
        template <class A, class T>
        XSIMD_INLINE batch<T, A> ceil(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            batch<T, A> truncated_self = trunc(self);
            return select(truncated_self < self, truncated_self + 1, truncated_self);
        }

        // floor
        template <class A, class T>
        XSIMD_INLINE batch<T, A> floor(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            batch<T, A> truncated_self = trunc(self);
            return select(truncated_self > self, truncated_self - 1, truncated_self);
        }

        // round
        template <class A, class T>
        XSIMD_INLINE batch<T, A> round(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            auto v = abs(self);
            auto c = ceil(v);
            auto cp = select(c - 0.5 > v, c - 1, c);
            return select(v > constants::maxflint<batch<T, A>>(), self, copysign(cp, self));
        }

        // trunc
        template <class A, class T, class = typename std::enable_if<std::is_integral<T>::value, void>::type>
        XSIMD_INLINE batch<T, A> trunc(batch<T, A> const& self, requires_arch<generic>) noexcept
        {
            return self;
        }
        template <class A>
        XSIMD_INLINE batch<float, A> trunc(batch<float, A> const& self, requires_arch<generic>) noexcept
        {
            return select(abs(self) < constants::maxflint<batch<float, A>>(), to_float(to_int(self)), self);
        }
        template <class A>
        XSIMD_INLINE batch<double, A> trunc(batch<double, A> const& self, requires_arch<generic>) noexcept
        {
            return select(abs(self) < constants::maxflint<batch<double, A>>(), to_float(to_int(self)), self);
        }

    }

}

#endif
