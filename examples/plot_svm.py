import numpy as np
import pylab as pl
from scikits.learn.svm import SVM

# import some data to play with
from scikits.learn.datasets.iris import load
SP, SW, PL, PW, LABELS = load()
X = np.c_[SP, SW]
Y = LABELS

def demo_plot_svm(h=.05, kernel_type='linear'):
    """
    Sample usage of Support Vector Machines to classify a sample.
    It will plot the decision surface and the support vectors.

    Parameters
    ----------
    h : float
        step size in the plot. Smaller values of h will give plots
        with higher resolution, but will take longer to compute.

    kernel_type: string
        specifies the kernel to use in support vectors. See the
        docstring of scikits.learn.svm.SVM for a complete list
        of available kernels.
    """
    # we create an instance of SVM and fit out data
    clf = SVM(kernel_type='linear')
    clf.fit(X, Y)

    # Plot the decision boundary. For that, we will asign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    x_min, x_max = SP.min()-1, SP.max()+1
    y_min, y_max = SW.min()-1, SW.max()+1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.pcolormesh(xx, yy, Z)

    # Plot also the training points
    pl.scatter(SP, SW, c=Y)
    # and the support vectors
    pl.scatter(clf.support_[:,0], clf.support_[:, 1], marker='+')
    pl.axis('tight')
    pl.show()


if __name__ == '__main__':
    demo_plot_svm()
