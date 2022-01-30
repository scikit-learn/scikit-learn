"""
Common location for shared fused types
"""

from numpy cimport (
    float32_t,
    float64_t,
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)

# All numeric types except complex
ctypedef fused numeric_t:
    int8_t
    int16_t
    int32_t
    int64_t

    uint8_t
    uint16_t
    uint32_t
    uint64_t

    float32_t
    float64_t

# All numeric types + object, doesn't include complex
ctypedef fused numeric_object_t:
    numeric_t
    object

# i64 + u64 + all float types
ctypedef fused iu_64_floating_t:
    float64_t
    float32_t
    int64_t
    uint64_t

# i64 + u64 + all float types + object
ctypedef fused iu_64_floating_obj_t:
    iu_64_floating_t
    object
