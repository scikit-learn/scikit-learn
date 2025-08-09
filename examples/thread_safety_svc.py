"""
Example to demonstrate thread-safety issue with SVC(probability=True).
Run multiple threads to see possible inconsistent behavior or errors.
"""

import threading
from sklearn.svm import SVC
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split

def train_svc():
    X, y = load_iris(return_X_y=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    clf = SVC(probability=True)
    clf.fit(X_train, y_train)
    preds = clf.predict_proba(X_test)
    print(preds[:5])

threads = []
for i in range(4):  # Running 4 threads in parallel
    t = threading.Thread(target=train_svc)
    threads.append(t)
    t.start()

for t in threads:
    t.join()
