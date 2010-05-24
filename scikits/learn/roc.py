import numpy as np

def roc(y, probas_):
    """compute Receiver operating characteristic (ROC)"""
    y = y.ravel()
    probas_ = probas_.ravel()
    n_samples = y.size
    thresholds = np.sort(np.unique(probas_))[::-1]
    tpr = [] # True positive rate
    fpr = [] # False positive rate
    n_pos = float(np.sum(y==1)) # nb of true positive
    n_neg = float(np.sum(y==0)) # nb of true negative
    for t in thresholds:
        tpr.append(np.sum(y[probas_>=t]==1) / n_pos)
        fpr.append(np.sum(y[probas_>=t]==0) / n_neg)

    return fpr, tpr, thresholds

def auc(x, y):
    """Compute Area Under the Curve (AUC)
    using the trapezoidal rule
    """
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    h = np.diff(fpr)
    area = np.sum(h * (y[1:]+y[:-1])) / 2.0
    return area

if __name__ == '__main__':
    import pylab as pl
    from scikits.learn import svm, datasets
    import random
    random.seed(0)

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
    fpr, tpr, thresholds = roc(y[half:], probas_[:,0])
    roc_auc = auc(fpr, tpr)
    print "Area under the ROC curve : %f" % roc_auc

    pl.figure(-1)
    pl.clf()
    pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    pl.hold('on')
    pl.plot([0, 1], [0, 1], 'k--')
    pl.hold('off')
    pl.xlim([0.0,1.0])
    pl.ylim([0.0,1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('Receiver operating characteristic example')
    pl.legend(loc="lower right")
    pl.show()