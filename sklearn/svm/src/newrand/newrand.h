/*
   Creation, 2020:
   - New random number generator using a mersenne twister + tweaked lemire
     postprocessor. This fixed a convergence issue on windows targets for
     libsvm and liblinear.
     Sylvain Marie, Schneider Electric
     See <https://github.com/scikit-learn/scikit-learn/pull/13511#issuecomment-481729756>
 */
#ifndef _NEWRAND_H
#define _NEWRAND_H

#ifdef __cplusplus
#include <random>  // needed for cython to generate a .cpp file from newrand.h
extern "C" {
#endif

// Scikit-Learn-specific random number generator replacing `rand()` originally
// used in LibSVM / LibLinear, to ensure the same behaviour on windows-linux,
// with increased speed
// - (1) Init a `mt_rand` object
std::mt19937 mt_rand(std::mt19937::default_seed);

// - (2) public `set_seed()` function that should be used instead of `srand()` to set a new seed.
void set_seed(unsigned custom_seed) {
    mt_rand.seed(custom_seed);
}

// - (3) New internal `bounded_rand_int` function, used instead of rand() everywhere.
inline uint32_t bounded_rand_int(uint32_t range) {
    // "LibSVM / LibLinear Original way" - make a 31bit positive
    // random number and use modulo to make it fit in the range
    // return abs( (int)mt_rand()) % range;

    // "Better way": tweaked Lemire post-processor
    // from http://www.pcg-random.org/posts/bounded-rands.html
    uint32_t x = mt_rand();
    uint64_t m = uint64_t(x) * uint64_t(range);
    uint32_t l = uint32_t(m);
    if (l < range) {
        uint32_t t = -range;
        if (t >= range) {
            t -= range;
            if (t >= range)
                t %= range;
        }
        while (l < t) {
            x = mt_rand();
            m = uint64_t(x) * uint64_t(range);
            l = uint32_t(m);
        }
    }
    return m >> 32;
}

#ifdef __cplusplus
}
#endif

#endif /* _NEWRAND_H */
