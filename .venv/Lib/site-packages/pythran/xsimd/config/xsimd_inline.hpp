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

#ifndef XSIMD_INLINE_HPP
#define XSIMD_INLINE_HPP

#if defined(__GNUC__)
#define XSIMD_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define XSIMD_INLINE inline __forceinline
#else
#define XSIMD_INLINE inline
#endif

#endif
