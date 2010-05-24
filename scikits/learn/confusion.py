import numpy as np

def confusion_matrix(y, y_):
    """
    compute confusion matrix
    to evaluate the accuracy of a classification result

    By definition a confusion matrix cm is such that

    cm[i,j] is equal to the number of observations known to be in group i
    but predicted to be in group j

    Parameters
    ==========

    y : array
        true targets

    y_ : array
        estimated targets

    """

    # removing possible NaNs in targets (they are ignored)
    clean_y = y[np.isfinite(y)].ravel()
    clean_y_ = y_[np.isfinite(y_)].ravel()

    labels = np.r_[np.unique(clean_y).ravel(),np.unique(clean_y_).ravel()]
    labels = np.unique(labels)
    n_labels = labels.size

    cm = np.empty((n_labels,n_labels))
    for i, label_i in enumerate(labels):
        for j, label_j in enumerate(labels):
            cm[i,j] = np.sum(np.logical_and(y==label_i, y_==label_j))

    return cm

if __name__ == '__main__':
    import pylab as pl
    from scikits.learn import svm, datasets
    import random
    random.seed(0)

    # import some data to play with
    iris = datasets.load_iris()
    X = iris.data
    y = iris.target
    n_samples, n_features = X.shape
    p = range(n_samples)
    random.shuffle(p)
    X, y = X[p], y[p]
    half = int(n_samples/2)

    classifier = svm.SVC(kernel='linear')
    y_ = classifier.fit(X[:half],y[:half]).predict(X[half:])

    cm = confusion_matrix(y[half:], y_)

    print cm

    pl.matshow(cm)
    pl.title('Confusion matrix')
    pl.colorbar()
    pl.show()