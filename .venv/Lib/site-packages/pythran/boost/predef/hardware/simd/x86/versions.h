/*
Copyright Charly Chevalier 2015
Copyright Joel Falcou 2015
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_PREDEF_HARDWARE_SIMD_X86_VERSIONS_H
#define BOOST_PREDEF_HARDWARE_SIMD_X86_VERSIONS_H

#include <boost/predef/version_number.h>

/*`
 Those defines represent x86 SIMD extensions versions.

 [note You *MUST* compare them with the predef `BOOST_HW_SIMD_X86`.]
 */

// ---------------------------------

/*`
 [heading `BOOST_HW_SIMD_X86_MMX_VERSION`]

 The [@https://en.wikipedia.org/wiki/MMX_(instruction_set) MMX] x86 extension
 version number.

 Version number is: *0.99.0*.
 */
#define BOOST_HW_SIMD_X86_MMX_VERSION BOOST_VERSION_NUMBER(0, 99, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSE_VERSION`]

 The [@https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions SSE] x86 extension
 version number.

 Version number is: *1.0.0*.
 */
#define BOOST_HW_SIMD_X86_SSE_VERSION BOOST_VERSION_NUMBER(1, 0, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSE2_VERSION`]

 The [@https://en.wikipedia.org/wiki/SSE2 SSE2] x86 extension version number.

 Version number is: *2.0.0*.
 */
#define BOOST_HW_SIMD_X86_SSE2_VERSION BOOST_VERSION_NUMBER(2, 0, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSE3_VERSION`]

 The [@https://en.wikipedia.org/wiki/SSE3 SSE3] x86 extension version number.

 Version number is: *3.0.0*.
 */
#define BOOST_HW_SIMD_X86_SSE3_VERSION BOOST_VERSION_NUMBER(3, 0, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSSE3_VERSION`]

 The [@https://en.wikipedia.org/wiki/SSSE3 SSSE3] x86 extension version number.

 Version number is: *3.1.0*.
 */
#define BOOST_HW_SIMD_X86_SSSE3_VERSION BOOST_VERSION_NUMBER(3, 1, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSE4_1_VERSION`]

 The [@https://en.wikipedia.org/wiki/SSE4#SSE4.1 SSE4_1] x86 extension version
 number.

 Version number is: *4.1.0*.
 */
#define BOOST_HW_SIMD_X86_SSE4_1_VERSION BOOST_VERSION_NUMBER(4, 1, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_SSE4_2_VERSION`]

 The [@https://en.wikipedia.org/wiki/SSE4##SSE4.2 SSE4_2] x86 extension version
 number.

 Version number is: *4.2.0*.
 */
#define BOOST_HW_SIMD_X86_SSE4_2_VERSION BOOST_VERSION_NUMBER(4, 2, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_AVX_VERSION`]

 The [@https://en.wikipedia.org/wiki/Advanced_Vector_Extensions AVX] x86
 extension version number.

 Version number is: *5.0.0*.
 */
#define BOOST_HW_SIMD_X86_AVX_VERSION BOOST_VERSION_NUMBER(5, 0, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_FMA3_VERSION`]

 The [@https://en.wikipedia.org/wiki/FMA_instruction_set FMA3] x86 extension
 version number.

 Version number is: *5.2.0*.
 */
#define BOOST_HW_SIMD_X86_FMA3_VERSION BOOST_VERSION_NUMBER(5, 2, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_AVX2_VERSION`]

 The [@https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2 AVX2]
 x86 extension version number.

 Version number is: *5.3.0*.
 */
#define BOOST_HW_SIMD_X86_AVX2_VERSION BOOST_VERSION_NUMBER(5, 3, 0)

/*`
 [heading `BOOST_HW_SIMD_X86_MIC_VERSION`]

 The [@https://en.wikipedia.org/wiki/Xeon_Phi MIC] (Xeon Phi) x86 extension
 version number.

 Version number is: *9.0.0*.
 */
#define BOOST_HW_SIMD_X86_MIC_VERSION BOOST_VERSION_NUMBER(9, 0, 0)

#endif
