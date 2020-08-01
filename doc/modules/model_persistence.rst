.. _model_persistence:

=================
Model persistence
=================

After training a scikit-learn model, it is desirable to have a way to persist
the model for future use without having to retrain. The following sections give
you some hints on how to persist a scikit-learn model.

.. note::

   Remember that, once exported in a persistent format, the model could only be
   used for predictions and it cannot be refitted.

Binary serialization
--------------------

It is possible to save a model in scikit-learn by using Python's built-in
persistence model, namely `pickle
<https://docs.python.org/3/library/pickle.html>`_::

  >>> from sklearn import svm
  >>> from sklearn import datasets
  >>> clf = svm.SVC()
  >>> X, y= datasets.load_iris(return_X_y=True)
  >>> clf.fit(X, y)
  SVC()

  >>> import pickle
  >>> s = pickle.dumps(clf)
  >>> clf2 = pickle.loads(s)
  >>> clf2.predict(X[0:1])
  array([0])
  >>> y[0]
  0

In the specific case of scikit-learn, it may be better to use joblib's
replacement of pickle (``dump`` & ``load``), which is more efficient on
objects that carry large numpy arrays internally as is often the case for
fitted scikit-learn estimators, but can only pickle to the disk and not to a
string::

  >>> from joblib import dump, load
  >>> dump(clf, 'filename.joblib') # doctest: +SKIP

Later you can load back the pickled model (possibly in another Python process)
with::

  >>> clf = load('filename.joblib') # doctest:+SKIP

.. note::

   ``dump`` and ``load`` functions also accept file-like object
   instead of filenames. More information on data persistence with Joblib is
   available `here
   <https://joblib.readthedocs.io/en/latest/persistence.html>`_.

.. _persistence_limitations:

Security & maintainability limitations
......................................

pickle (and joblib by extension), has some issues regarding maintainability
and security. Because of this,

* Never unpickle untrusted data as it could lead to malicious code being 
  executed upon loading.
* While models saved using one version of scikit-learn might load in 
  other versions, this is entirely unsupported and inadvisable. It should 
  also be kept in mind that operations performed on such data could give
  different and unexpected results.

In order to rebuild a similar model with future versions of scikit-learn,
additional metadata should be saved along the pickled model:

* The training data, e.g. a reference to an immutable snapshot
* The python source code used to generate the model
* The versions of scikit-learn and its dependencies
* The cross validation score obtained on the training data

This should make it possible to check that the cross-validation score is in the
same range as before.

Since a model internal representation may be different on two different
architectures, dumping a model on one architecture and loading it on
another architecture is not supported.

If you want to know more about these issues and explore other possible
serialization methods, please refer to this
`talk by Alex Gaynor
<https://pyvideo.org/video/2566/pickles-are-for-delis-not-software>`_.

Interoperable formats
---------------------

For production and quality control needs, exporting the model in `Predictive
Model Markup Language (PMML)
<http://dmg.org/pmml/v4-4-1/GeneralStructure.html>`_ or `Open Neural Network
Exchange <https://onnx.ai/>`_ format
would be a better approach than using `pickle`.

PMML is an extension of the `XML
<https://fr.wikipedia.org/wiki/Extensible_Markup_Language>`_ document standard
defined to represent data models together with the data used to generate them.
Being human and machine readable,
PMML is a good option for model validation on different platforms and
long term archiving. On the other hand, as XML in general, its verbosity does
not help in production when performance is critical.

ONNX has been developed to improve the usability of the interoperable
representation of data models. It aims to facilitate the conversion of the data
models between different machine learning frameworks, and to improve their
portability on different computing architectures. More details are available
from the `ONNX tutorial <https://onnx.ai/get-started.html>`_. 
