/* Fast low level math functions */

/* Fast logarithm (2-based) approximation.
 *
 * This function is a brutal approximation and will easy have errors up
 * to 10 %
 * */
inline float fast_log2(float f) 
{
    /* Grab the mantisse of the floating point */
    int i = (*(int *)&f);
    /* Approximate log on the mantisse */
    return (((i&0x7f800000)>>23)-0x7f)+(i&0x007fffff)/(float)0x800000;
}


/* Fast logarithm approximation. */
inline float fast_log(float val)
{
    return (fast_log2(val) * 0.69314718f);
}

