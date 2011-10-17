/* Fast low level math functions */

#include <math.h>


/****************************************************************************
 * Fast logarithm (2-based) approximation.
 *
 * This function is a brutal approximation and will easy have errors up
 * to 10 %
 *
 * Notes
 * ======
 *
 * The strategy relies on two steps. First the exponent of the floating
 * point representation of the number is retrieved. This is cheap because
 * there is no log involved, just accessing the IEEE float
 * representation. This exponent gives us the log2 rounded below.
 *
 * To get additional accuracy, a linear interpolation between the value
 * rounding below and the value rounding up is used.
 *
 * Reference
 * ==========
 *
 * The core idea is taken from Schraudolph 1998,
 * "A Fast, Compact Approximation of the Exponential Function" 
 * http://nic.linotune.com/pubs/Schraudolph99.pdf
 *
 ****************************************************************************/
static inline float fast_log2(float f) 
{
    /* Grab the raw representation the floating point */
    int exponent;
    float mantissa = frexpf(f, &exponent);
    /* The exponent is one too much (exponent of 1 is 1) and the mantissa
     * is scaled between -.5 and .5 while we want it scaled from 0 to 2 */
    return exponent - 1 + 2*(mantissa - .5);
}


/* Fast logarithm approximation. */
static inline float fast_log(float val)
{
    /* 0.69314718f is log(2) */
    return (fast_log2(val) * 0.69314718f);
}

