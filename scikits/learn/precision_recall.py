import numpy as np

def precision_recall(y, probas_):
    """compute Precision-Recall"""
    y = y.ravel()
    probas_ = probas_.ravel()
    thresholds = np.sort(np.unique(probas_))
    precision = []
    recall = []
    for t in thresholds:
        true_pos = np.sum(y[probas_>=t]==1)
        false_pos = np.sum(y[probas_>=t]==0)
        false_neg = np.sum(y[probas_<t]==1)
        precision.append(true_pos / float(true_pos + false_pos))
        recall.append(true_pos / float(true_pos + false_neg))

    precision.append(1.0)
    recall.append(0.0)
    return precision, recall

if __name__ == '__main__':
    import pylab as pl
    from scikits.learn import svm, datasets
    import random
    random.seed(0)
    np.random.seed(0)

    # import some data to play with
    iris = datasets.load_iris()
    X = iris.data
    y = iris.target
    X, y = X[y!=2], y[y!=2]
    n_samples, n_features = X.shape
    p = range(n_samples)
    random.shuffle(p)
    X, y = X[p], y[p]
    half = int(n_samples/2)

    # Add noisy features
    X = np.c_[X,np.random.randn(n_samples, 200*n_features)]

    # Run classifier
    classifier = svm.SVC(kernel='linear', probability=True)
    probas_ = classifier.fit(X[:half],y[:half]).predict_proba(X[half:])

    # Compute ROC curve and area the curve
    precision, recall = precision_recall(y[half:], probas_[:,0])

    pl.figure(-1)
    pl.clf()
    pl.plot(recall, precision, label='Precision-Recall curve')
    pl.xlabel('Recall')
    pl.ylabel('Precision')
    pl.ylim([0.0,1.05])
    pl.xlim([0.0,1.0])
    pl.title('Precision-Recall example')
    pl.legend(loc="lower left")
    pl.show()