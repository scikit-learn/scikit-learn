# Import necessary packages
import pandas as pd
import numpy as np
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import FixOutliers
import pytest

# Test FixOutliers class
def test_FixOutliers():

    # Load Boston housing dataset
    data = load_boston()
    X = pd.DataFrame(data['data'], columns=data['feature_names'])
    y = data['target']

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    # Test 1: Test approach parameter
    with pytest.raises(ValueError):
        fix_outliers = FixOutliers(approach='invalid_approach')

    # Test 2: Test treatment parameter
    with pytest.raises(ValueError):
        fix_outliers = FixOutliers()
        fix_outliers.fit_transform(X_train, 'CRIM', treatment='invalid_treatment')

    # Test 3: Test imputation parameter
    with pytest.raises(ValueError):
        fix_outliers = FixOutliers()
        fix_outliers.fit_transform(X_train, 'CRIM', treatment='impute', imputation='invalid_imputation')

    # Test 4: Test consider_outliers parameter
    with pytest.raises(ValueError):
        fix_outliers = FixOutliers()
        fix_outliers.fit_transform(X_train, 'CRIM', treatment='impute', consider_outliers='invalid_bool')

    # Test 5: Test removal of outliers using z-score approach
    fix_outliers = FixOutliers(approach='z_score')
    X_train_cleaned = fix_outliers.fit_transform(X_train, 'CRIM', treatment='remove', threshold=2)
    assert (X_train_cleaned['CRIM'] < 2).all()

    # Test 6: Test imputation of outliers using z-score approach
    fix_outliers = FixOutliers(approach='z_score')
    X_train_imputed = fix_outliers.fit_transform(X_train, 'CRIM', treatment='impute', imputation='mean', consider_outliers=True, threshold=2)
    assert np.isclose(X_train_imputed['CRIM'].mean(), X_train['CRIM'].mean(), rtol=0.01)

    # Test 7: Test removal of outliers using interquartile range approach
    fix_outliers = FixOutliers(approach='interquartile_range')
    X_train_cleaned = fix_outliers.fit_transform(X_train, 'CRIM', treatment='remove')
    Q1 = X_train['CRIM'].quantile(0.25)
    Q3 = X_train['CRIM'].quantile(0.75)
    IQR = Q3 - Q1
    assert (X_train_cleaned['CRIM'] >= Q1 - 1.5 * IQR).all()
    assert (X_train_cleaned['CRIM'] <= Q3 + 1.5 * IQR).all()

    # Test 8: Test imputation of outliers using interquartile range approach
    fix_outliers = FixOutliers(approach='interquartile_range')
    X_train_imputed = fix_outliers.fit_transform(X_train, 'CRIM', treatment='impute', imputation='median', consider_outliers=True)
    assert np.isclose(X_train_imputed['CRIM'].median(), X_train['CRIM'].median(), rtol=0.01)
