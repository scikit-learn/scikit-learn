"""
Quick benchmark to compare vanilla scikit-learn RandomForestClassifier vs sklearnex-accelerated RandomForest.
Run: python benchmarks/benchmark_histogram_vs_intel.py

Requirements: scikit-learn, optionally sklearnex (pip install scikit-learn-intelex or pip install sklearnex)
"""
import time
import numpy as np
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

TRY_INTEL = True

try:
    if TRY_INTEL:
        # sklearnex patches sklearn to use Intel optimized kernels where available.
        import sklearnex
        sklearnex.patch_sklearn()
        using_intel = True
    else:
        using_intel = False
except Exception:
    using_intel = False


def run_rf(X, y, n_estimators=100, n_jobs=-1, random_state=0):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    clf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=n_jobs, random_state=random_state)
    t0 = time.perf_counter()
    clf.fit(X_train, y_train)
    t1 = time.perf_counter()
    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    return (t1 - t0), acc


def main():
    print('Note: this benchmark measures wall-clock training time for RandomForestClassifier.fit')
    sizes = [50_000, 100_000, 200_000]
    n_features = 50
    n_informative = 30

    for n_samples in sizes:
        print('\n=== n_samples = {:,} ==='.format(n_samples))
        X, y = make_classification(n_samples=n_samples, n_features=n_features,
                                   n_informative=n_informative, n_redundant=0,
                                   n_classes=2, random_state=0)
        t, acc = run_rf(X, y, n_estimators=100)
        print('training_time(s) = {:.3f}, test_acc = {:.4f} (using_intel={})'.format(t, acc, using_intel))

    if using_intel:
        print('\nIntel-optimized sklearn used (sklearnex.patch_sklearn()).')
    else:
        print('\nIntel-optimized sklearn NOT used (sklearnex not installed or patch failed).')


if __name__ == '__main__':
    main()