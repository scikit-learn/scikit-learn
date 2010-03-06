from scikits.learn import svm
import numpy as np
import pylab as pl
from datetime import datetime

def bench(n=200): # (number of points)/2

    np.random.seed(0)
    data1 = np.random.randn(n, 2) + -2.0
    data2 = np.random.randn(n, 2) + 2.0
    data = np.concatenate((data1, data2))
    labels = [-1]*n + [1]*n

    h = .1 # step size
    x = np.arange(-6, 6, h)
    y = np.arange(-6, 6, h)
    xx, yy = np.meshgrid(x, y)
    Z = np.c_[xx.ravel(), yy.ravel()]

    tstart = datetime.now()
    clf = svm.SVM(kernel_type='linear')
    clf.fit(data, labels)
    SV = clf.support_
    Z = clf.predict(Z)
    delta = (datetime.now() - tstart)
    print 'training points: ', len(labels)
    print 'number of support vectors ', SV.shape[0]
    print 'number of samples to classify: ', Z.shape[0]
    print 'and scikit-learn has done this in ...'
    print '\t', delta 

    Z = Z.reshape(np.shape(xx))

    ax = pl.subplot(111)
    pl.pcolormesh(xx, yy, Z)

    pl.scatter(data1[:,0], data1[:,1], c='blue')
    pl.scatter(data2[:,0], data2[:,1], c='red')
    pl.scatter(SV[:,0], SV[:,1], c='black')
    pl.axis('tight')
    pl.show()

if __name__ == '__main__':
    bench()
