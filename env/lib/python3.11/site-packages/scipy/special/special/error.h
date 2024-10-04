#pragma once

// should be included from config.h, but that won't work until we've cleanly separated out the C and C++ parts of the
// code
#ifdef __CUDACC__
#define SPECFUN_HOST_DEVICE __host__ __device__
#else
#define SPECFUN_HOST_DEVICE
#endif

typedef enum {
    SF_ERROR_OK = 0,    /* no error */
    SF_ERROR_SINGULAR,  /* singularity encountered */
    SF_ERROR_UNDERFLOW, /* floating point underflow */
    SF_ERROR_OVERFLOW,  /* floating point overflow */
    SF_ERROR_SLOW,      /* too many iterations required */
    SF_ERROR_LOSS,      /* loss of precision */
    SF_ERROR_NO_RESULT, /* no result obtained */
    SF_ERROR_DOMAIN,    /* out of domain */
    SF_ERROR_ARG,       /* invalid input parameter */
    SF_ERROR_OTHER,     /* unclassified error */
    SF_ERROR__LAST
} sf_error_t;

#ifdef __cplusplus

#include <complex>

namespace special {

#ifndef SP_SPECFUN_ERROR
SPECFUN_HOST_DEVICE inline void set_error(const char *func_name, sf_error_t code, const char *fmt, ...) {
    // nothing
}
#else
void set_error(const char *func_name, sf_error_t code, const char *fmt, ...);
#endif

template <typename T>
void set_error_and_nan(const char *name, sf_error_t code, T &value) {
    if (code != SF_ERROR_OK) {
        set_error(name, code, nullptr);

        if (code == SF_ERROR_DOMAIN || code == SF_ERROR_OVERFLOW || code == SF_ERROR_NO_RESULT) {
            value = std::numeric_limits<T>::quiet_NaN();
        }
    }
}

template <typename T>
void set_error_and_nan(const char *name, sf_error_t code, std::complex<T> &value) {
    if (code != SF_ERROR_OK) {
        set_error(name, code, nullptr);

        if (code == SF_ERROR_DOMAIN || code == SF_ERROR_OVERFLOW || code == SF_ERROR_NO_RESULT) {
            value.real(std::numeric_limits<T>::quiet_NaN());
            value.imag(std::numeric_limits<T>::quiet_NaN());
        }
    }
}

} // namespace special

#endif
