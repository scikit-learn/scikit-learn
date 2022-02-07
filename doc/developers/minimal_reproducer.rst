.. _minimal_reproducer:

===========================================
Craft a minimal reproducer for scikit-learn
===========================================


Whether if you are submitting a bug report, designing a pytest for a feature you
are contributing to, or simply posting a question in a forum, being able to
craft minimal, reproducible examples (or minimal, workable examples) is the key
to communicating with the community.

There are very good guidelines on the internet such as `this StackOverflow
document <https://stackoverflow.com/help/mcve>`_ or `this blogpost by Matthew
Rocklin <https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_
on how to craft Minimal Complete Verifiable Examples (MCVE). Our goal is not to
be repetitive with those references but rather to provide a step-by-step guide
on how to narrow down a bug until you have reached the shortest possible code to
reproduce it.

The first step before submitting a bug report to scikit-learn is to read the
`Issue template
<https://github.com/scikit-learn/scikit-learn/blob/main/.github/ISSUE_TEMPLATE/bug_report.yml>`_.
It is already quite informative about the information you will be asked to provide.


.. _good_practices:

Good practices
==============

In this section we will focus on the **Steps/Code to Reproduce** section of
the Issue template. We will start with a snippet of code that already provides a
failing example but that has room for readability improvement. We then craft a
MCVE from it.

**Example**

.. code-block:: python

    # I am currently working in a ML project and when I tried to fit a
    # GradientBoostingRegressor instance to my_data.csv I get a UserWarning:
    # "X has feature names, but DecisionTreeRegressor was fitted without feature
    # names". You can get a copy of my dataset from www.this_link.com and verify
    # my features do have names. The problem seems to arise during fit when I pass
    # an integer to the n_iter_no_change parameter.

    import pandas as pd
    df = pd.read_csv('my_data.csv')
    X = df["feature_name"] # my features do have names
    y = df["target"]

    # We set random_state=42 for the train_test_split
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42)

    from sklearn.ensemble import GradientBoostingRegressor

    # An instance with default n_iter_no_change raises no error nor warnings
    gbdt = GradientBoostingRegressor(random_state=0)
    gbdt.fit(X_train, y_train)

    # the bug appears when I change the value for n_iter_no_change
    gbdt = GradientBoostingRegressor(random_state=0, n_iter_no_change=5)
    gbdt.fit(X_train, y_train)

Provide a failing code example with minimal annotations
-------------------------------------------------------

Please take into account that things may get lost in translation for non-native
english speakers. It is better that all the necessary information can be read
from the code itself. Besides, by this point you already provided a concise
description in the **Describe the bug** section of the Issue template.

**Better example**

.. code-block:: python

    import pandas as pd
    df = pd.read_csv('my_data.csv')
    X = df["feature_name"]
    y = df["target"]

    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42)

    from sklearn.ensemble import GradientBoostingRegressor

    gbdt = GradientBoostingRegressor(random_state=0)
    gbdt.fit(X_train, y_train) # compiles

    gbdt = GradientBoostingRegressor(random_state=0, n_iter_no_change=5)
    gbdt.fit(X_train, y_train) # raises warning


Boil down your script to something as small as possible
-------------------------------------------------------

You have to ask yourself which lines of code are relavant and which are not for
rasing the bug. Erasing unnecessary lines of code will help you and other
contributors narrow down the bug. In this example the warning has nothing to
do with the `train_test_split`, nor the `random_state`.

**Even better example**

.. code-block:: python

    import pandas as pd
    df = pd.read_csv('my_data.csv')
    X = df["feature_name"]
    y = df["target"]

    from sklearn.ensemble import GradientBoostingRegressor

    gbdt = GradientBoostingRegressor()
    gbdt.fit(X, y) # compiles

    gbdt = GradientBoostingRegressor(n_iter_no_change=5)
    gbdt.fit(X, y) # raises warning


**DO NOT** report your data unless it is extremely necessary
------------------------------------------------------------

The idea is to make the code as self-contained as possible. For doing so, you
can use a :ref:`synth_data`. It can be generated using numpy, pandas or the
:mod:`sklearn.datasets` module. Most of the times the bug is not related to a
particular structure of your data. Even if it is, try to find an available
dataset that has similar characteristics to yours and that reproduces the
problem. In this particular case, we are interested in data that has labeled
feature names.

**Even better example**

.. code-block:: python

    import pandas as pd
    df = pd.DataFrame(
        {
            "feature_name": [-12.32, 1.43, 30.00],
            "target": [72, 55, 32],
        }
    )
    X = df["feature_name"]
    y = df["target"]

    from sklearn.ensemble import GradientBoostingRegressor

    gbdt = GradientBoostingRegressor()
    gbdt.fit(X, y) # compiles

    gbdt = GradientBoostingRegressor(n_iter_no_change=5)
    gbdt.fit(X, y) # raises warning

The above steps can be implemented in a different order than the progression we
show in this example. Take into account that having a pipeline that makes sense
as a model is not important when creating a MCVE.

.. _synth_data:


Synthetic dataset
=================

Before chosing a particular synthetic dataset, first you have to identify the
type of problem you are solving: Is it classification, regression, clustering,
etc?

Once that you narrowed down the type of problem, you need to provide a synthetic
dataset accordingly. Most of the times you only need a minimalistic dataset.
Here is a non-exhaustive list of tools that may help you.

Numpy
-----

Numpy tools such as `random.randn
<https://numpy.org/doc/stable/reference/random/generated/numpy.random.randn.html>`_
and `random.randint
<https://numpy.org/doc/stable/reference/random/generated/numpy.random.randint.html>`_
can be used to create dummy numeric data

- regression

Regressions take continuous numeric data as features and target

.. code-block:: python

    import numpy as np

    rng = np.random.RandomState(0)
    n_samples, n_features = 5, 5
    X = rng.randn(n_samples, n_features)
    y = rng.randn(n_samples)

A similar snippet can be used as synthetic data when testing scaling tools such
as :class:`sklearn.preprocessing.StandardScaler`.

- classification

If the bug is not raised during encoding a categorical variable, you can feed
numeric data to a classifier. Just remember to ensure that the target is indeed
an integer.

.. code-block:: python

    import numpy as np

    rng = np.random.RandomState(0)
    n_samples, n_features = 5, 5
    X = rng.randn(n_samples, n_features)
    y = rng.randint(0, 2, n_samples)

If you need to test encoding, you may prefer to start from non-numeric data. In
such case you may use `numpy.random.choice
<https://numpy.org/doc/stable/reference/random/generated/numpy.random.choice.html>`_.

.. code-block:: python

    import numpy as np

    n_samples, n_features = 50, 5
    X = rng.randn(n_samples, n_features)
    y = np.random.choice(
        ["male", "female", "other"], size=n_samples, p=[0.49, 0.49, 0.02]
    )

Pandas
------

Some scikit-learn objets expect pandas dataframes as input. In this case you can
transform numpy arrays into pandas objects using `pandas.DataFrame
<https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`_, or
`pandas.Series
<https://pandas.pydata.org/docs/reference/api/pandas.Series.html>`_.

.. code-block:: python

    import numpy as np
    import pandas as pd

    rng = np.random.RandomState(0)
    n_samples, n_features = 5, 5
    X = rng.randn(n_samples, n_features)
    y = rng.randn(n_samples)
    X = pd.DataFrame(X)
    y = pd.Series(y)

In addition, scikit-learn includes various :ref:`sample_generators` that can be
used to build artificial datasets of controlled size and complexity.

`make_regression`
-----------------

As hinted by the name, :class:`sklearn.datasets.make_regression` produces
regression targets with noise as an optionally-sparse random linear combination
of random features.

.. code-block:: python

    from sklearn.datasets import make_regression

    X, y = make_regression(n_samples=1000, n_features=20)

`make_classification`
---------------------

:class:`sklearn.datasets.make_classification` creates multiclass datasets with multiple Gaussian
clusters per class. Noise can be introduced by means of correlated, redundant or
uninformative features.

.. code-block:: python

    from sklearn.datasets import make_classification

    X, y = make_classification(
        n_features=2, n_redundant=0, n_informative=2, n_clusters_per_class=1
    )

`make_blobs`
------------

Similarly to `make_classification`, :class:`sklearn.datasets.make_blobs` creates multiclass
datasets using normally-distributed clusters of points. It provides greater
control regarding the centers and standard deviations of each cluster, and
therefore it is useful to demonstrate clustering.

.. code-block:: python

    from sklearn.datasets import make_blobs

    X, y = make_blobs(n_samples=10, centers=3, n_features=2)

Dataset loading utilities
-------------------------

You can use the :ref:`datasets` to load and fetch several popular reference
datasets. This option is useful when the bug relates to the particular structure
of the data, e.g. dealing with missing values or image recognition.

.. code-block:: python

    from sklearn.datasets import load_breast_cancer

    cancer = load_breast_cancer()
    X = cancer.data
    y = cancer.target

Formatting
==========

As already mentioned, the key to communication is the readability of the code
and good formatting can really improve it

Use markdown
--------------------------------------------------------------------------------------------------------------------------------------------------------------------

To format code or text into its own distinct block, use triple backticks.
`Markdown
<https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax>`_
supports an optional language identifier to enable syntax highlighting in your
fenced code block. For example::

    ```python
    from sklearn.datasets import make_blobs

    n_samples = 100
    n_components = 3
    X, y = make_blobs(n_samples=n_samples, centers=n_components)
    ```

will render a python formatted snippet as follows

.. code-block:: python

    from sklearn.datasets import make_blobs

    n_samples = 100
    n_components = 3
    X, y = make_blobs(n_samples=n_samples, centers=n_components)

It is not necessary to create several blocks of code when submitting a bug
report. Remember other reviewers are going to copy-paste your code.

In the section `**Actual results**` of your `Issue template
<https://github.com/scikit-learn/scikit-learn/blob/main/.github/ISSUE_TEMPLATE/bug_report.yml>`_.
you are asked to provide the error message including the full traceback of the exception. In this
case, use the `python-traceback` qualifier. For example::

    ```python-traceback
    ---------------------------------------------------------------------------
    TypeError                                 Traceback (most recent call last)
    <ipython-input-1-a674e682c281> in <module>
        4 vectorizer = CountVectorizer(input=docs, analyzer='word')
        5 lda_features = vectorizer.fit_transform(docs)
    ----> 6 lda_model = LatentDirichletAllocation(
        7     n_topics=10,
        8     learning_method='online',

    TypeError: __init__() got an unexpected keyword argument 'n_topics'
    ```

yields the following when rendered:

.. code-block:: python

    ---------------------------------------------------------------------------
    TypeError                                 Traceback (most recent call last)
    <ipython-input-1-a674e682c281> in <module>
        4 vectorizer = CountVectorizer(input=docs, analyzer='word')
        5 lda_features = vectorizer.fit_transform(docs)
    ----> 6 lda_model = LatentDirichletAllocation(
        7     n_topics=10,
        8     learning_method='online',

    TypeError: __init__() got an unexpected keyword argument 'n_topics'


Try to follow the `PEP 8 convention <https://www.python.org/dev/peps/pep-0008/>`_
---------------------------------------------------------------------------------

The convention in a nutshell:
    - Try to limit all lines to a maximum of 79 characters
    - use blank lines to separate groups of related functions
    - blank lines may be omitted between a bunch of related lines of code

**Example**

The MCVE we created in the :ref:`good_practices` section is easier to read than
the equally working MCVE here below

.. code-block:: python

    import pandas as pd
    df = pd.DataFrame({"feature_name": [-12.32, 1.43, 30.00], "target": [72, 55, 32]})
    X = df["feature_name"]
    y = df["target"]
    from sklearn.ensemble import GradientBoostingRegressor
    gbdt = GradientBoostingRegressor()
    gbdt.fit(X, y) # compiles
    gbdt = GradientBoostingRegressor(n_iter_no_change=5)
    gbdt.fit(X, y) # raises warning
