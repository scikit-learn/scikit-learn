.. Places parent toc into the sidebar

:parenttoc: True

.. _model_persistence:

=================
Model persistence
=================

After training a scikit-learn model, it is desirable to have a way to persist
the model for future use without having to retrain. Based on your use-case,
there are a few different ways to persist a scikit-learn model, and here we
help you decide which suits you best. In order to make a decision, you need to
answer the following questions:

1. Do you need the python object after persistence, or do you only need to
   persist in order to serve the model and get predictions out of it?

If the answer to the above question is that you only need to serve the model
and no further investigation on the python object itself is required, then
:ref:`ONNX <onnx_persistence>` might be the best fit for you. Note that not all
models are supported by ONNX.

2. In case ONNX is not suitable for your use-case, do you absolutely trust the
   source of the model, or are there any security concerns regarding where the
   persisted model comes from?

If you have security concerns, then you should consider using :ref:`skops.io
<skops_persistence>` which gives you back the python object, but unlike
`pickle` based persistence solutions, loading the persisted model doesn't
automatically allow arbitrary code execution. Note that this requires manual
investigation of the persisted file, which :mod:`skops.io` allows you to do.

The other solutions assume you absolutely trust the source of the file to be
loaded, as they are all susceptible to arbitrary code execution upon loading
the persisted file since they all use the pickle protocol under the hood.

3. Do you care about the performance of loading the model, and sharing it
   between processes where a memory mapped object on disk is beneficial?

If yes, then you can consider using :ref:`joblib <pickle_persistence>`. If this
is not a major concern for you, then you can use the built-in :mod:`pickle`
module.

4. Did you try :mod:`pickle` or :mod:`joblib` and found that the model cannot
   be persisted? This can happen for instance when you have user defined
   functions in your model.

If yes, then you can use `cloudpickle`_ which can serialize certain objects
which cannot be serialized by :mod:`pickle` or :mod:`joblib`.


Workflow Overview
-----------------

In a typical workflow, the first step is to train the model using scikit-learn
and scikit-learn compatible libraries. Note that support for scikit-learn and
third party estimators varies across the different persistence methods.

Train and Persist the Model
...........................

Creating an appropriate model depends on your use-case. As an example, here we
train a :class:`sklearn.ensemble.HistGradientBoostingClassifier` on the iris
dataset::

  >>> from sklearn import ensemble
  >>> from sklearn import datasets
  >>> clf = ensemble.HistGradientBoostingClassifier()
  >>> X, y = datasets.load_iris(return_X_y=True)
  >>> clf.fit(X, y)
  HistGradientBoostingClassifier()

Once the model is trained, you can persist it using your desired method, and
then you can load the model in a separate environment and get predictions from
it given input data. Here there are two major paths depending on how you
persist and plan to serve the model:

- :ref:`ONNX <onnx_persistence>`: You need an `ONNX` runtime and an environment
  with appropriate dependencies installed to load the model and use the runtime
  to get predictions. This environment can be minimal and does not necessarily
  even require `python` to be installed.

- :mod:`skops.io`, :mod:`pickle`, :mod:`joblib`, `cloudpickle`_: You need a
  python environment with the appropriate dependencies installed to load the
  model and get predictions from it. This environment should have the same
  packages and the same versions as the environment where the model was
  trained. Note that none of these methods support loading a model trained with
  a different version of scikit-learn, and possibly different versions of other
  dependencies such as `numpy` and `scipy`.


.. _onnx_persistence:

ONNX
----

`ONNX`, or `Open Neural Network Exchange <https://onnx.ai/>`__ format is best
suitable in use-cases where one needs to persist the model and then use the
persisted artifact to get predictions without the need to load the python
object itself. It is also useful in cases where the serving environment needs
to be lean and minimal, since the `ONNX` runtime does not require `python`.

`ONNX` is a binary serialization of the model. It has been developed to improve
the usability of the interoperable representation of data models. It aims to
facilitate the conversion of the data models between different machine learning
frameworks, and to improve their portability on different computing
architectures. More details are available from the `ONNX tutorial
<https://onnx.ai/get-started.html>`__. To convert scikit-learn model to `ONNX`
`sklearn-onnx <http://onnx.ai/sklearn-onnx/>`__ has been developed. However,
not all scikit-learn models are supported, and it is limited to the core
scikit-learn and does not support most third party estimators. One can write a
custom converter for third party or custom estimators, but the documentation to
do that is sparse and it might be challenging to do so.

.. _skops_persistence:

`skops.io`
----------

:mod:`skops.io` avoids using :mod:`pickle` and only loads files which have types
and references to functions which are trusted either by default or by the user.
Therefore it provides a more secure format than :mod:`pickle`, :mod:`joblib`,
and `cloudpickle`_.


|details-start|
**Using skops**
|details-split|

The API is very similar to :mod:`pickle`, and you can persist your models as
explained in the `documentation
<https://skops.readthedocs.io/en/stable/persistence.html>`__ using
:func:`skops.io.dump` and :func:`skops.io.dumps`::

    import skops.io as sio
    obj = sio.dump(clf, "filename.skops")

And you can load them back using :func:`skops.io.load` and
:func:`skops.io.loads`. However, you need to specify the types which are
trusted by you. You can get existing unknown types in a dumped object / file
using :func:`skops.io.get_untrusted_types`, and after checking its contents,
pass it to the load function::

    unknown_types = sio.get_untrusted_types(file="filename.skops")
    # investigate the contents of unknown_types, and only load if you trust
    # everything you see.
    clf = sio.load("filename.skops", trusted=unknown_types)

Please report issues and feature requests related to this format on the `skops
issue tracker <https://github.com/skops-dev/skops/issues>`__.

|details-end|

.. _pickle_persistence:

`pickle`, `joblib`, and `cloudpickle`
-------------------------------------

These three modules / packages, use the `pickle` protocol under the hood, but
come with slight variations:

- :mod:`pickle` is the native Python serialization module. It can serialize and
  deserialize any Python object, including custom Python classes and objects.
- :mod:`joblib` is more efficient than `pickle` when working with large machine
  learning models or large numpy arrays.
- `cloudpickle`_ can serialize certain objects which cannot be serialized
  by :mod:`pickle` or :mod:`joblib`, such as user defined functions and lambda
  functions.

|details-start|
**Using ``pickle``, ``joblib``, or ``cloudpickle``**
|details-split|

Depending on your use-case, you can choose one of these three methods to
persist and load your scikit-learn model, and they all follow the same API::

    # Here you can replace pickle with joblib or cloudpickle
    from pickle import dump
    with open('filename.pkl', 'wb') as f: dump(clf, f)

And later when needed, you can load the same object from the persisted file::

    # Here you can replace pickle with joblib or cloudpickle
    from pickle import load
    with open('filename.pkl', 'rb') as f: clf = load(f)

|details-end|

.. _persistence_limitations:

Security & Maintainability Limitations
--------------------------------------

:mod:`pickle` (and :mod:`joblib` and :mod:`clouldpickle` by extension), has
many documented security vulnerabilities and should only be used if the
artifact is coming from a trusted and verified source.

Also note that arbitrary computations can be represented using the `ONNX`
format, and therefore a sandbox used to serve models using `ONNX` also needs to
safeguard against computational and memory exploits.

Also note that there are no supported ways to load a model trained with a
different version of scikit-learn. While using :mod:`skops.io`, :mod:`joblib`,
:mod:`pickle`, or `cloudpickle`_, models saved using one version of
scikit-learn might load in other versions, however, this is entirely
unsupported and inadvisable. It should also be kept in mind that operations
performed on such data could give different and unexpected results.

In order to rebuild a similar model with future versions of scikit-learn,
additional metadata should be saved along the pickled model:

* The training data, e.g. a reference to an immutable snapshot
* The python source code used to generate the model
* The versions of scikit-learn and its dependencies
* The cross validation score obtained on the training data

This should make it possible to check that the cross-validation score is in the
same range as before.

Aside for a few exceptions, persisted models should be portable across
architectures assuming the same versions of dependencies and Python are used.
If you encounter an estimator that is not portable please open an issue on
GitHub. Persisted models are often deployed in production using containers
like Docker, in order to freeze the environment and dependencies.

If you want to know more about these issues please refer to these talks:

- `Adrin Jalali: Let's exploit pickle, and skops to the rescue! | PyData
  Amsterdam 2023 <https://www.youtube.com/watch?v=9w_H5OSTO9A>`__.
- `Alex Gaynor: Pickles are for Delis, not Software - PyCon 2014
  <https://pyvideo.org/video/2566/pickles-are-for-delis-not-software>`__.


.. _serving_environment:

Replicating the training environment in production
..................................................

If the versions of the dependencies used may differ from training to
production, it may result in unexpected behaviour and errors while using the
trained model. To prevent such situations it is recommended to use the same
dependencies and versions in both the training and production environment.
These transitive dependencies can be pinned with the help of `pip`, `mamba`,
`conda`, `poetry`, `conda-lock`, `pixi`, etc.

|details-start|
**InconsistentVersionWarning**
|details-split|

When an estimator is loaded with a scikit-learn version that is inconsistent
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


Serving the model artifact
..........................

The last step after training a scikit-learn model is serving the model.
Once the trained model is successfully loaded it can be served to manage
different prediction requests. This can involve deploying the model as a
web service using containerization, or other model deployment strategies,
according to the specifications. In the next sections, we will explore
different approaches to persist a trained scikit-learn model.


Summarizing the key points
--------------------------

Based on the different approaches for model persistence, the key points for
each approach can be summarized as follows:

* `ONNX`: It provides a uniform format for persisting any machine learning
  or deep learning model (other than scikit-learn) and is useful
  for model inference. It can however, result in compatibility issues with
  different frameworks.
* :mod:`skops.io`: Trained scikit-learn models can be easily shared and put
  into production using :mod:`skops.io`. It is more secure compared `pickle`
  based alternate approaches as it does not load arbitrary code unless
  explicitly asked by the user.
* :mod:`joblib`: Efficient storage and memory mapping techniques make it faster
  when working with large machine learning models or large numpy arrays.
  However, it may trigger the execution of malicious code while loading
  untrusted data.
* :mod:`pickle`: It is native to Python and any Python object can be serialized
  and deserialized using :mod:`pickle`, including custom Python classes and
  objects. While :mod:`pickle` can be used to easily save and load scikit-learn
  models, it may trigger the execution of malicious code while loading
  untrusted data.
* `cloudpickle`_: It is slower than :mod:`pickle` and :mod:`joblib`, and
  is more insecure than :mod:`pickle` and :mod:`joblib` since it can serialize
  arbitrary code. However, in certain cases it might be a last resort to
  persist certain models.

.. _cloudpickle: https://github.com/cloudpipe/cloudpickle
