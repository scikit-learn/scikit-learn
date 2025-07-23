#!/usr/bin/env python3
"""
Simple demonstration of the issue #31717 and its fix.
"""

import numpy as np
from collections import Counter

def original_most_frequent(array, extra_value, n_repeat):
    """Original implementation that would fail with incomparable types in ties"""
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

            # This would FAIL with incomparable types
            most_frequent_value = min(tied_values)  # This line causes TypeError
        else:
            # Use scipy mode for non-object dtypes
            most_frequent_value = array[0]  # Simplified
            most_frequent_count = 1
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
        # tie breaking - this would also fail with incomparable types
        return min(most_frequent_value, extra_value)  # This line causes TypeError


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
    """Fixed implementation that handles incomparable types in ties"""
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
            # Use scipy mode for non-object dtypes
            most_frequent_value = array[0]  # Simplified
            most_frequent_count = 1
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


def test_incomparable_types():
    """Test case to demonstrate the issue and fix"""
    print("Testing incomparable types that cause ties...")
    
    # Create array with incomparable types
    array = np.array([None, "string"], dtype=object)
    extra_value = 42  # Another incomparable type
    n_repeat = 1  # This creates a tie scenario
    
    print(f"Array: {array}")
    print(f"Extra value: {extra_value}")
    print(f"N repeat: {n_repeat}")
    print()
    
    # Test original implementation (would fail)
    print("Testing original implementation:")
    try:
        result = original_most_frequent(array, extra_value, n_repeat)
        print(f"SUCCESS: {result}")
    except Exception as e:
        print(f"FAILURE: {type(e).__name__}: {e}")
    
    print()
    
    # Test fixed implementation 
    print("Testing fixed implementation:")
    try:
        result = fixed_most_frequent(array, extra_value, n_repeat)
        print(f"SUCCESS: {result}")
    except Exception as e:
        print(f"FAILURE: {type(e).__name__}: {e}")


def test_safe_min_function():
    """Test the safe_min_with_fallback function directly"""
    print("\n" + "="*50)
    print("Testing safe_min_with_fallback function:")
    
    # Test case 1: Normal comparable values
    values1 = [3, 1, 4, 2]
    result1 = safe_min_with_fallback(values1)
    print(f"Comparable values {values1}: {result1}")
    
    # Test case 2: Incomparable types 
    values2 = [None, "string", 42]
    result2 = safe_min_with_fallback(values2)
    print(f"Incomparable values {values2}: {result2}")
    
    # Test case 3: Mix of types that can be string-compared
    values3 = [1, "2", None]
    result3 = safe_min_with_fallback(values3)
    print(f"Mixed types {values3}: {result3}")


if __name__ == "__main__":
    test_incomparable_types()
    test_safe_min_function()
