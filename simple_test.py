#!/usr/bin/env python3
"""
Simple test to understand the _most_frequent function issue with incomparable types
"""

import numpy as np
from collections import Counter

def original_most_frequent(array, extra_value, n_repeat):
    """Original implementation that fails with incomparable types"""
    # Compute the most frequent value in array only
    if array.size > 0:
        if array.dtype == object:
            # scipy.stats.mode is slow with object dtype array.
            # Python Counter is more efficient
            counter = Counter(array)
            most_frequent_count = counter.most_common(1)[0][1]
            # tie breaking similarly to scipy.stats.mode
            # THIS IS THE PROBLEMATIC LINE:
            most_frequent_value = min(
                value
                for value, count in counter.items()
                if count == most_frequent_count
            )
        else:
            # For non-object arrays, use scipy stats mode
            from scipy.stats import mode
            mode_result = mode(array, keepdims=True)
            most_frequent_value = mode_result.mode[0]
            most_frequent_count = mode_result.count[0]
    else:
        most_frequent_value = 0
        most_frequent_count = 0

    # Compare to array + [extra_value] * n_repeat
    if most_frequent_count == 0 and n_repeat == 0:
        return np.nan
    elif most_frequent_count < n_repeat:
        return extra_value
    elif most_frequent_count > n_repeat:
        return most_frequent_value
    elif most_frequent_count == n_repeat:
        # tie breaking similarly to scipy.stats.mode
        # THIS IS ANOTHER PROBLEMATIC LINE:
        return min(most_frequent_value, extra_value)

def test_incomparable_types():
    """Test cases that demonstrate the issue"""
    
    print("Testing incomparable types in _most_frequent function")
    
    # Test case 1: None vs string (tie situation)
    try:
        array = np.array([None, "string"], dtype=object)
        result = original_most_frequent(array, "another_string", 1)
        print(f"Test 1 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 1 FAILED: {type(e).__name__}: {e}")
    
    # Test case 2: Different object types with ties
    try:
        array = np.array([1, "string"], dtype=object)  # int vs string
        result = original_most_frequent(array, None, 1)
        print(f"Test 2 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 2 FAILED: {type(e).__name__}: {e}")
    
    # Test case 3: Complex objects
    try:
        array = np.array([[], "string"], dtype=object)  # list vs string
        result = original_most_frequent(array, {}, 1)  # dict as extra value
        print(f"Test 3 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 3 FAILED: {type(e).__name__}: {e}")

if __name__ == "__main__":
    test_incomparable_types()
