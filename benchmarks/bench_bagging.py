"""
Benchmark ``BaggingClassifier`` and ``BaggingRegressor`` with parallelism.
"""

from time import time

import numpy as np

from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.model_selection import train_test_split

n_samples = 5_000_000
dim = 10
n_classes = 10
X = np.random.randn(n_samples, dim)
y = np.random.randint(0, n_classes, (n_samples,))
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.99, random_state=42
)

# Warm up process executor pool if any:
clf = BaggingClassifier(n_jobs=-1)
clf.fit(X_train, y_train).predict(X_test)

print("BaggingClassifer")
clf = BaggingClassifier(n_jobs=-1)
start = time()
clf.fit(X_train, y_train)
print("Fit", time() - start)
start = time()
clf.predict(X_test)
print("Predict", time() - start)

clf = BaggingRegressor(n_jobs=-1)
print("\nBaggingRegressor")
start = time()
clf.fit(X_train, y_train)
print("Fit", time() - start)
start = time()
clf.predict(X_test)
print("Predict", time() - start)
