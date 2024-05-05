.. Places parent toc into the sidebar

:parenttoc: True

.. _model_persistence:

=================
Model persistence
=================

After training a scikit-learn model, it is desirable to have a way to persist
the model for future use without having to retrain. This can be accomplished
using `pickle <https://docs.python.org/3/library/pickle.html>`_, `joblib
<https://joblib.readthedocs.io/en/stable/>`_, `skops
<https://skops.readthedocs.io/en/stable/>`_, `ONNX <https://onnx.ai/>`_,
or `PMML <https://dmg.org/pmml/v4-4-1/GeneralStructure.html>`_. In most cases
`pickle` can be used to persist a trained scikit-learn model. Once all
transitive scikit-learn dependencies have been pinned, the trained model can
then be loaded and executed under conditions similar to those in which it was
originally pinned. The following sections will give you some hints on how to
persist a scikit-learn model and will provide details on what each alternative
can offer.

Workflow Overview
-----------------

In this section we present a general workflow on how to persist a
scikit-learn model. We will demonstrate this with a simple example using
Python's built-in persistence module, namely `pickle
<https://docs.python.org/3/library/pickle.html>`_.

Storing the model in an artifact
................................

Once the model training process in completed, the trained model can be stored
as an artifact with the help of `pickle`. The model can be saved using the
process of serialization, where the Python object hierarchy is converted into
a byte stream. We can persist a trained model in the following manner::

  >>> from sklearn import svm
  >>> from sklearn import datasets
  >>> import pickle
  >>> clf = svm.SVC()
  >>> X, y = datasets.load_iris(return_X_y=True)
  >>> clf.fit(X, y)
  SVC()
  >>> s = pickle.dumps(clf)

Replicating the training environment in production
..................................................

The versions of the dependencies used may differ from training to production.
This may result in unexpected behaviour and errors while using the trained
model. To prevent such situations it is recommended to use the same
dependencies and versions in both the training and production environment.
These transitive dependencies can be pinned with the help of `pip`, `conda`,
`poetry`, `conda-lock`, `pixi`, etc.

.. note::

    To execute a pickled scikit-learn model in a reproducible environment it is
    advisable to pin all transitive scikit-learn dependencies. This prevents
    any incompatibility issues that may arise while trying to load the pickled
    model. You can read more about persisting models with `pickle` over
    :ref:`here <persisting_models_with_pickle>`.

Loading the model artifact
..........................

The saved scikit-learn model can be loaded using `pickle` for future use
without having to re-train the entire model from scratch. The saved model
artifact can be unpickled by converting the byte stream into an object
hierarchy. This can be done with the help of `pickle` as follows::

  >>> clf2 = pickle.loads(s) # doctest:+SKIP
  >>> clf2.predict(X[0:1]) # doctest:+SKIP
  array([0])
  >>> y[0] # doctest:+SKIP
  0

Serving the model artifact
..........................

The last step after training a scikit-learn model is serving the model.
Once the trained model is successfully loaded it can be served to manage
different prediction requests. This can involve deploying the model as a
web service using containerization, or other model deployment strategies,
according to the specifications. In the next sections, we will explore
different approaches to persist a trained scikit-learn model.

.. _persisting_models_with_pickle:

Persisting models with pickle
-----------------------------

As demonstrated in the previous section, `pickle` uses serialization and
deserialization to persist scikit-learn models. Instead of using `dumps` and
`loads`, `dump` and `load` can also be used in the following way::

  >>> from sklearn.tree import DecisionTreeClassifier
  >>> from sklearn import datasets
  >>> clf = DecisionTreeClassifier()
  >>> X, y = datasets.load_iris(return_X_y=True)
  >>> clf.fit(X, y)
  DecisionTreeClassifier()
  >>> from pickle import dump, load
  >>> with open('filename.pkl', 'wb') as f: dump(clf, f) # doctest:+SKIP
  >>> with open('filename.pkl', 'rb') as f: clf2 = load(f) # doctest:+SKIP
  >>> clf2.predict(X[0:1]) # doctest:+SKIP
  array([0])
  >>> y[0]
  0

For applications that involve writing and loading the serialized object to or
from a file, `dump` and `load` can be used instead of `dumps` and `loads`. When
file operations are not required the pickled representation of the object can
be returned as a bytes object with the help of the `dumps` function. The
reconstituted object hierarchy of the pickled data can then be returned using
the `loads` function.

Persisting models with joblib
-----------------------------

In the specific case of scikit-learn, it may be better to use joblib's
replacement of pickle (``dump`` & ``load``), which is more efficient on
objects that carry large numpy arrays internally as is often the case for
fitted scikit-learn estimators, but can only pickle to the disk and not to a
string::

  >>> from joblib import dump, load
  >>> dump(clf, 'filename.joblib') # doctest:+SKIP

Later you can load back the pickled model (possibly in another Python process)
with::

  >>> clf = load('filename.joblib') # doctest:+SKIP

.. note::

   ``dump`` and ``load`` functions also accept file-like object
   instead of filenames. More information on data persistence with Joblib is
   available `here
   <https://joblib.readthedocs.io/en/latest/persistence.html>`_.

|details-start|
**InconsistentVersionWarning**
|details-split|

When an estimator is unpickled with a scikit-learn version that is inconsistent
with the version the estimator was pickled with, a
:class:`~sklearn.exceptions.InconsistentVersionWarning` is raised. This warning
can be caught to obtain the original version the estimator was pickled with::

  from sklearn.exceptions import InconsistentVersionWarning
  warnings.simplefilter("error", InconsistentVersionWarning)

  try:
      est = pickle.loads("model_from_prevision_version.pickle")
  except InconsistentVersionWarning as w:
      print(w.original_sklearn_version)

|details-end|

.. _persistence_limitations:

Security & maintainability limitations for pickle and joblib
------------------------------------------------------------

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

Aside for a few exceptions, pickled models should be portable across
architectures assuming the same versions of dependencies and Python are used.
If you encounter an estimator that is not portable please open an issue on
GitHub. Pickled models are often deployed in production using containers, like
Docker, in order to freeze the environment and dependencies.

If you want to know more about these issues and explore other possible
serialization methods, please refer to this
`talk by Alex Gaynor
<https://pyvideo.org/video/2566/pickles-are-for-delis-not-software>`_.

Persisting models with a more secure format using skops
-------------------------------------------------------

`skops <https://skops.readthedocs.io/en/stable/>`__ provides a more secure
format via the :mod:`skops.io` module. It avoids using :mod:`pickle` and only
loads files which have types and references to functions which are trusted
either by default or by the user.

|details-start|
**Using skops**
|details-split|

The API is very similar to ``pickle``, and
you can persist your models as explain in the `docs
<https://skops.readthedocs.io/en/stable/persistence.html>`__ using
:func:`skops.io.dump` and :func:`skops.io.dumps`::

    import skops.io as sio
    obj = sio.dumps(clf)

And you can load them back using :func:`skops.io.load` and
:func:`skops.io.loads`. However, you need to specify the types which are
trusted by you. You can get existing unknown types in a dumped object / file
using :func:`skops.io.get_untrusted_types`, and after checking its contents,
pass it to the load function::

    unknown_types = sio.get_untrusted_types(data=obj)
    clf = sio.loads(obj, trusted=unknown_types)

If you trust the source of the file / object, you can pass ``trusted=True``::

    clf = sio.loads(obj, trusted=True)

Please report issues and feature requests related to this format on the `skops
issue tracker <https://github.com/skops-dev/skops/issues>`__.

|details-end|

Persisting models with interoperable formats
--------------------------------------------

For reproducibility and quality control needs, when different architectures
and environments should be taken into account, exporting the model in
`Open Neural Network
Exchange <https://onnx.ai/>`_ format or `Predictive Model Markup Language
(PMML) <https://dmg.org/pmml/v4-4-1/GeneralStructure.html>`_ format
might be a better approach than using `pickle` alone.
These are helpful where you may want to use your model for prediction in a
different environment from where the model was trained.

ONNX is a binary serialization of the model. It has been developed to improve
the usability of the interoperable representation of data models.
It aims to facilitate the conversion of the data
models between different machine learning frameworks, and to improve their
portability on different computing architectures. More details are available
from the `ONNX tutorial <https://onnx.ai/get-started.html>`_.
To convert scikit-learn model to ONNX a specific tool `sklearn-onnx
<http://onnx.ai/sklearn-onnx/>`_ has been developed.

PMML is an implementation of the `XML
<https://en.wikipedia.org/wiki/XML>`_ document standard
defined to represent data models together with the data used to generate them.
Being human and machine readable,
PMML is a good option for model validation on different platforms and
long term archiving. On the other hand, as XML in general, its verbosity does
not help in production when performance is critical.
To convert scikit-learn model to PMML you can use for example `sklearn2pmml
<https://github.com/jpmml/sklearn2pmml>`_ distributed under the Affero GPLv3
license.

Summarizing the keypoints
-------------------------

Based on the different approaches for model persistence, the keypoints for each
approach can be summarized as follows:

* `pickle`: It is native to Python and any Python object can be serialized and
  deserialized using `pickle`, including custom Python classes and objects.
  While `pickle` can be used to easily save and load scikit-learn models,
  unpickling of untrusted data might lead to security issues.
* `joblib`: Efficient storage and memory mapping techniques make it faster
  when working with large machine learning models or large numpy arrays. However,
  it may trigger the execution of malicious code while loading untrusted data.
* `skops`: Trained scikit-learn models can be easily shared and put into
  production using `skops`. It is more secure compared to alternate approaches
  as it allows users to load data from trusted sources. It however, does not
  allow for persistence of arbitrary Python code.
* `ONNX`: It provides a uniform format for persisting any machine learning
  or deep learning model (other than scikit-learn) and is useful
  for model inference. It can however, result in compatibility issues with
  different frameworks.
* `PMML`: Platform independent format that can be used to persist models
  and reduce the risk of vendor lock-ins. The complexity and verbosity of
  this format might make it harder to use for larger models.