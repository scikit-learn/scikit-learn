#!/usr/bin/env python3
"""
Fixed version of _most_frequent function that handles incomparable types
"""

import numpy as np
from collections import Counter

def safe_min_with_fallback(values):
    """
    Safely find minimum value, handling incomparable types.
    
    For incomparable types, fall back to deterministic selection based on:
    1. String representation
    2. Object id as last resort
    """
    try:
        return min(values)
    except TypeError:
        # If direct comparison fails, use string representation for comparison
        try:
            return min(values, key=lambda x: str(x))
        except (TypeError, ValueError):
            # Last resort: use object id for deterministic but arbitrary selection
            return min(values, key=lambda x: id(x))

def fixed_most_frequent(array, extra_value, n_repeat):
    """
    Fixed version of _most_frequent that handles incomparable types
    """
    # Compute the most frequent value in array only
    if array.size > 0:
        if array.dtype == object:
            # scipy.stats.mode is slow with object dtype array.
            # Python Counter is more efficient
            counter = Counter(array)
            most_frequent_count = counter.most_common(1)[0][1]
            
            # Get all values with the most frequent count
            tied_values = [
                value for value, count in counter.items()
                if count == most_frequent_count
            ]
            
            # Use safe comparison for tie breaking
            most_frequent_value = safe_min_with_fallback(tied_values)
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
        # tie breaking - use safe comparison
        return safe_min_with_fallback([most_frequent_value, extra_value])

def test_fixed_version():
    """Test the fixed version with incomparable types"""
    
    print("Testing FIXED _most_frequent function with incomparable types")
    
    # Test case 1: None vs string (tie situation)
    try:
        array = np.array([None, "string"], dtype=object)
        result = fixed_most_frequent(array, "another_string", 1)
        print(f"Test 1 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 1 FAILED: {type(e).__name__}: {e}")
    
    # Test case 2: Different object types with ties
    try:
        array = np.array([1, "string"], dtype=object)  # int vs string
        result = fixed_most_frequent(array, None, 1)
        print(f"Test 2 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 2 FAILED: {type(e).__name__}: {e}")
    
    # Test case 3: Normal comparable case (should still work)
    try:
        array = np.array([1, 2, 1], dtype=object)  # integers
        result = fixed_most_frequent(array, 3, 1)
        print(f"Test 3 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 3 FAILED: {type(e).__name__}: {e}")
        
    # Test case 4: String ties (should work normally)
    try:
        array = np.array(["a", "b"], dtype=object)
        result = fixed_most_frequent(array, "c", 1)
        print(f"Test 4 SUCCESS: {result}")
    except Exception as e:
        print(f"Test 4 FAILED: {type(e).__name__}: {e}")

if __name__ == "__main__":
    test_fixed_version()
