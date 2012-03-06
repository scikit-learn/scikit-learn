=========================
Putting it all together
=========================

..  Imports
    >>> import pylab as pl
    >>> import numpy as np

Pipelining
============

We have seen that some estimators can transform data, and some estimators
can predict variables. We can create combined estimators:

.. image:: ../../auto_examples/images/plot_digits_pipe_1.png
   :target: ../../auto_examples/plot_digits_pipe.html
   :scale: 65
   :align: right

::

    >>> from sklearn import svm, decomposition, datasets

    >>> svc = svm.SVC(kernel='linear')
    >>> pca = decomposition.PCA()
    >>> from sklearn.pipeline import Pipeline
    >>> pipe = Pipeline(steps=[('pca', pca), ('svc', svc)])

    >>> digits = datasets.load_digits()
    >>> X_digits = digits.data
    >>> y_digits = digits.target
    >>> pca.fit(X_digits, y_digits)
    PCA(copy=True, n_components=None, whiten=False)
    >>> pl.plot(pca.explained_variance_) # doctest: +ELLIPSIS
    [<matplotlib.lines.Line2D object at ...>]

Parameters of pipelines can be set using '__' separated parameter names::

    >>> pipe._set_params(pca__n_components=30)
    Pipeline(steps=[('pca', PCA(copy=True, n_components=30, whiten=False)), ('svc', SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3, gamma=0.0,
      kernel='linear', probability=False, scale_C=True, shrinking=True,
      tol=0.001))])

    >>> pca.n_components
    30

    >>> from sklearn.grid_search import GridSearchCV
    >>> n_components = [10, 15, 20, 30, 40, 50, 64]
    >>> Cs = np.logspace(-4, 4, 16)
    >>> estimator = GridSearchCV(pipe,
    ...                          dict(pca__n_components=n_components,
    ...                               svc__C=Cs),
    ...                          n_jobs=-1)
    >>> estimator.fit(X_digits, y_digits) # doctest: +ELLIPSIS
    GridSearchCV(cv=None,...


Face recognition with eigenfaces
=================================

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", aka LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/

.. literalinclude:: ../../auto_examples/applications/plot_face_recognition.py

.. |prediction| image:: ../../auto_examples/applications/images/plot_face_recognition_1.png
   :target: ../../auto_examples/applications/plot_face_recognition.html
   :scale: 50

.. |eigenfaces| image:: ../../auto_examples/applications/images/plot_face_recognition_2.png
   :target: ../../auto_examples/applications/plot_face_recognition.html
   :scale: 50

.. list-table::
   :class: centered

   *

     - |prediction|

     - |eigenfaces|

   * 

     - **Prediction**

     - **Eigenfaces**

Expected results for the top 5 most represented people in the dataset::

                     precision    recall  f1-score   support

  Gerhard_Schroeder       0.91      0.75      0.82        28
    Donald_Rumsfeld       0.84      0.82      0.83        33
         Tony_Blair       0.65      0.82      0.73        34
       Colin_Powell       0.78      0.88      0.83        58
      George_W_Bush       0.93      0.86      0.90       129

        avg / total       0.86      0.84      0.85       282


Open problem: stock market structure
=====================================

Can we predict the variation in stock prices for Google?

.. literalinclude:: ../../auto_examples/applications/plot_stock_market.py
    :lines: 1-167



