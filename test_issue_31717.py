#!/usr/bin/env python3
"""
Test to reproduce issue #31717: SimpleImputer fails in "most_frequent" if incomparable types only if ties

The issue occurs when:
1. Using SimpleImputer with strategy="most_frequent" 
2. Data contains incomparable types (e.g., different object types)
3. There are ties in the frequency counts
4. The `min()` function fails on incomparable values during tie breaking
"""

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.impute._base import _most_frequent

def test_issue_reproduction():
    """Test case to reproduce the issue with incomparable types causing ties"""
    
    print("Testing issue #31717: SimpleImputer fails with incomparable types in ties")
    
    # Test case 1: Incomparable object types with ties
    # This should fail with the current implementation when there are ties
    try:
        # Create data with incomparable types that will have ties
        # Using a mix of different object types that can't be compared with min()
        X = np.array([
            [None, "a", 1],  # None is incomparable with other types
            [None, "b", 2], 
            ["string", "a", 3],  # string vs None - incomparable 
            ["string", "b", 4],
        ], dtype=object)
        
        print("Input data:")
        print(X)
        print("Data types:", [type(x) for row in X for x in row])
        
        imputer = SimpleImputer(missing_values=None, strategy="most_frequent")
        result = imputer.fit_transform(X)
        print("SUCCESS: Imputation completed without error")
        print("Result:")
        print(result)
        
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")
        
    print("\n" + "="*60 + "\n")
    
    # Test case 2: Direct test of _most_frequent function
    print("Testing _most_frequent function directly:")
    
    try:
        # Test with incomparable types that create ties
        array = np.array([None, "string"], dtype=object)
        result = _most_frequent(array, "another_string", 1)  # This creates a tie
        print(f"SUCCESS: _most_frequent returned: {result}")
        
    except Exception as e:
        print(f"ERROR in _most_frequent: {type(e).__name__}: {e}")
        
    print("\n" + "="*60 + "\n")
    
    # Test case 3: Test with numeric vs string mix (also incomparable)
    print("Testing with numeric vs string mix:")
    
    try:
        X = np.array([
            [np.nan, 1, "a"],
            [np.nan, 2, "b"], 
            [5, 1, "a"],  # 5 vs 1 (comparable), but in ties it might mix with strings
            [10, 2, "b"],
        ], dtype=object)
        
        print("Input data:")
        print(X)
        
        imputer = SimpleImputer(missing_values=np.nan, strategy="most_frequent")
        result = imputer.fit_transform(X)
        print("SUCCESS: Imputation completed")
        print("Result:")
        print(result)
        
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")

if __name__ == "__main__":
    test_issue_reproduction()
