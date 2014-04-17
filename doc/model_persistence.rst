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
persistence model, namely `pickle <http://docs.python.org/library/pickle.html>`_::

  >>> from sklearn import svm
  >>> from sklearn import datasets
  >>> clf = svm.SVC()
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> clf.fit(X, y)  # doctest: +NORMALIZE_WHITESPACE
  SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3, gamma=0.0,
    kernel='rbf', max_iter=-1, probability=False, random_state=None,
    shrinking=True, tol=0.001, verbose=False)

  >>> import pickle
  >>> s = pickle.dumps(clf)
  >>> clf2 = pickle.loads(s)
  >>> clf2.predict(X[0])
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
   contained in the `clf` object is serialized as a separate file on the
   filesystem. All files are required in the same folder when reloading the
   model with joblib.load.


Security & maintainability limitations
--------------------------------------

You must be aware that pickle has some issues regarding maintainability and
security. From the **maintainability** point of view, you should take care the
issues that may arise if you upgrade your sklearn library while still loading a
model that was trained with a previous version, the model may have a code
structure that could not be compatible with newer versions and thus, don't work.
The same issue could also happen if you upgrade numpy or scipy versions.

A good practice is to save the scikit-learn, numpy and scipy versions to know
exactly what versions have been used to generate the model. You can do that, for
example, by executing a ``pip freeze`` command and saving the output to a text
file which should be stored together with your pickles.
Also, save a snapshot of your data to make it possible to retrain the model
if incompatibility issues arise when upgrading the libraries.

Regarding **security** issues, you may know that pickle is implemented with a
stack machine that executes instructions. As a difference with other
serialization methods like JSON, BSON, YAML, etc, which are all data oriented,
pickle is instruction oriented. Pickle serializes objects by persisting a set of
instructions that will be then executed at deserialization time in order to
reconstruct your objects. In fact, as part of the deserialization process,
pickle could call any arbitrary function, which opens up security
vulnerabilities against any malicious data or exploits.

Here is the warning from the official pickle documentation:

.. warning::

    The pickle module is not intended to be secure against erroneous or
    maliciously constructed data.  Never unpickle data received from an untrusted
    or unauthenticated source.
    
If you want to know more about these issues and explore other possible
serialization methods, please refer to this
`talk by Alex Gaynor <http://pyvideo.org/video/2566/pickles-are-for-delis-not-software>`_.  
  
  
  