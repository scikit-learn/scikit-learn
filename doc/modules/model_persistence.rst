.. _model_persistence:

=================
Model persistence
=================

After training a scikit-learn model, it is desirable to have a way to persist
the model for future use without having to retrain. The following section gives
you an example of how to persist a model with pickle. We'll also review a few
security and maintainability issues when working with pickle serialization.


Persistence example
-------------------

It is possible to save a model in the scikit by using Python's built-in
persistence model, namely `pickle <http://docs.python.org/2/library/pickle.html>`_::

  >>> from sklearn import svm
  >>> from sklearn import datasets
  >>> clf = svm.SVC()
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> clf.fit(X, y)  # doctest: +NORMALIZE_WHITESPACE
  SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
      decision_function_shape=None, degree=3, gamma='auto', kernel='rbf',
      max_iter=-1, probability=False, random_state=None, shrinking=True,
      tol=0.001, verbose=False)

  >>> import pickle
  >>> s = pickle.dumps(clf)
  >>> clf2 = pickle.loads(s)
  >>> clf2.predict(X[0:1])
  array([0])
  >>> y[0]
  0

In the specific case of the scikit, it may be more interesting to use
joblib's replacement of pickle (``joblib.dump`` & ``joblib.load``),
which is more efficient on objects that carry large numpy arrays internally as
is often the case for fitted scikit-learn estimators, but can only pickle to the
disk and not to a string::

  >>> from sklearn.externals import joblib
  >>> joblib.dump(clf, 'filename.pkl') # doctest: +SKIP

Later you can load back the pickled model (possibly in another Python process)
with::

  >>> clf = joblib.load('filename.pkl') # doctest:+SKIP

.. note::

   joblib.dump returns a list of filenames. Each individual numpy array
   contained in the ``clf`` object is serialized as a separate file on the
   filesystem. All files are required in the same folder when reloading the
   model with joblib.load.


Security & maintainability limitations
--------------------------------------

pickle (and joblib by extension), has some issues regarding maintainability
and security. Because of this,

* Never unpickle untrusted data
* Models saved in one version of scikit-learn might not load in another
  version.

In order to rebuild a similar model with future versions of scikit-learn,
additional metadata should be saved along the pickled model:

* The training data, e.g. a reference to a immutable snapshot
* The python source code used to generate the model
* The versions of scikit-learn and its dependencies
* The cross validation score obtained on the training data

This should make it possible to check that the cross-validation score is in the
same range as before.

If you want to know more about these issues and explore other possible
serialization methods, please refer to this
`talk by Alex Gaynor <http://pyvideo.org/video/2566/pickles-are-for-delis-not-software>`_.
