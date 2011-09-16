"""
To run this, you'll need to have installed.

  * scikit-learn 

Does two benchmarks

First, we fix a training set, increase the number of
samples to classify and plot number of classified samples as a
function of time.

In the second benchmark, we increase the number of dimensions of the
training set, classify a sample and plot the time taken as a function
of the number of dimensions.
"""
import numpy as np
import pylab as pl
import gc
from datetime import datetime
from time import time

from sklearn import datasets
from sklearn import metrics

from sklearn.utils import shuffle

# to store the results
scikit_classifier_results = []
scikit_regressor_results = []

mu_second = 0.0 + 10**6 # number of microseconds in a second


def bench_scikit_tree_classifier(X, Y):
    """
    bench with scikit-learn decision tree classifier
    """
    import sklearn
    from sklearn.tree import DecisionTreeClassifier

    #gc.collect()

    # start time
    #tstart = datetime.now()
    clf = DecisionTreeClassifier()
    clf.fit(X, Y).predict(X)
    #delta = (datetime.now() - tstart)
    # stop time

    #scikit_classifier_results.append(delta.seconds + delta.microseconds/mu_second)

def bench_scikit_tree_regressor(X, Y):
    """
    bench with scikit-learn decision tree regressor
    """
    import sklearn
    from sklearn.tree import DecisionTreeRegressor

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = DecisionTreeRegressor()
    clf.fit(X, Y).predict(X)
    delta = (datetime.now() - tstart)
    # stop time

    scikit_regressor_results.append(delta.seconds + delta.microseconds/mu_second)

def profile(n_samples=1000, dim=50, K=2):
    np.random.seed(13)
    X = np.random.randn(n_samples, dim)
    Y = np.random.randint(0, K, (n_samples,))
    bench_scikit_tree_classifier(X, Y)


def repeat(f):
    def wrapper(*args, **kargs):
        scores = []
        for i in range(10):
            scores.append(f(*args, random_state=i, **kargs))
        return np.mean(scores)
    return wrapper


def bench_madelon():
    X_train = np.loadtxt("/home/pprett/corpora/madelon/madelon_train.data")
    y_train = np.loadtxt("/home/pprett/corpora/madelon/madelon_train.labels")
    X_test = np.loadtxt("/home/pprett/corpora/madelon/madelon_valid.data")
    y_test = np.loadtxt("/home/pprett/corpora/madelon/madelon_valid.labels")
    from sklearn.tree import DecisionTreeClassifier
    clf = DecisionTreeClassifier(max_depth=100, min_split=5)
    t0 = datetime.now()
    clf.fit(X_train, y_train)
    delta = (datetime.now() - t0)
    score = np.mean(clf.predict(X_test) == y_test)
#    clf.export_to_graphviz()
    return score

def bench_arcene():
    X_train = np.loadtxt("/home/pprett/corpora/arcene/arcene_train.data")
    y_train = np.loadtxt("/home/pprett/corpora/arcene/arcene_train.labels")
    X_test = np.loadtxt("/home/pprett/corpora/arcene/arcene_valid.data")
    y_test = np.loadtxt("/home/pprett/corpora/arcene/arcene_valid.labels")
    from sklearn.tree import DecisionTreeClassifier
    clf = DecisionTreeClassifier(max_depth=100, min_split=5)
    t0 = datetime.now()
    clf.fit(X_train, y_train)
    delta = (datetime.now() - t0)
    score = np.mean(clf.predict(X_test) == y_test)
    return score

@repeat
def bench_boston(random_state=None):
    from sklearn import datasets
    from sklearn.utils import shuffle
    boston = datasets.load_boston()
    X, y = shuffle(boston.data, boston.target, random_state=random_state)
    offset = int(X.shape[0] * 0.9)
    X_train = X[:offset]
    y_train = y[:offset]
    X_test = X[offset:]
    y_test = y[offset:]
    from sklearn.tree import DecisionTreeRegressor
    clf = DecisionTreeRegressor()
    t0 = datetime.now()
    clf.fit(X_train, y_train)
    delta = (datetime.now() - t0)
    mse = np.mean((clf.predict(X_test) - y_test)**2.0)
    #clf.export_to_graphviz()

    return mse

@repeat
def bench_friedman1(random_state=None):
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    from sklearn.tree import DecisionTreeRegressor
    clf = DecisionTreeRegressor()
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test)**2.0)
    return mse

@repeat
def bench_friedman2(random_state=None):
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    from sklearn.tree import DecisionTreeRegressor
    clf = DecisionTreeRegressor()
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test)**2.0)    
    return mse

@repeat
def bench_friedman3(random_state=None):
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    from sklearn.tree import DecisionTreeRegressor
    clf = DecisionTreeRegressor()
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test)**2.0)
    return mse

if __name__ == '__main__':
    #bench_madelon()
    #bench_arcene()
    print "Boston", bench_boston()
    print "Friedman#1", bench_friedman1()
    print "Friedman#2", bench_friedman2()
    print "Friedman#3", bench_friedman3()
