.. _df_output_transform:

===========================================================
Pandas/Polars Output for Transformers with `set_output` API
===========================================================

.. currentmodule:: sklearn

This part of the user guide explains how scikit-learn supports tabular data.


Propagation of Feature Names
============================

By default, scikit-learn :term:`transformers` (estimators with a :meth:`transform`
method) return numpy arrays (sometimes also sparse arrays). Because numpy arrays do
not provide names for the indices of axes/dimensions, prior to version 1.0
the :class:`pipeline.Pipeline` did not know how to propagate feature names:

- The single step estimators did not know how to handle incoming feature names.
- The pipeline did not know how to pass feature names from step to step.

In practice, a lot of use cases start with tabular data like a `pandas dataframe
<https://pandas.pydata.org/docs/user_guide/dsintro.html#dataframe>`_ or a
`polars dataframe <https://docs.pola.rs/api/python/stable/reference/dataframe/index.html>`__
which have column/feature names.

A first step to support this important use case was made by the addition of the
:class:`compose.ColumnTransformer` in :ref:`version 0.20 <changes_0_20>`.
It acts as a gateway to apply different estimators on the different features. Most
notably it understands incoming feature names.

It was then properly solved by `SLEP007: Feature names, their generation and the API
<https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep007/proposal.html>`__
and fully implemented in :ref:`version 1.1 <release_notes_1_1>`, see
the :ref:`release highlights 1.0 <feature_names_in_release_highlights_1_0_0>` and
:ref:`release highlights 1.1 <get_feature_names_out_release_highlights_1_1_0>`.
When an estimator is passed a dataframe during :term:`fit`, the estimator will
set a `feature_names_in_` attribute containing the feature names. It understands pandas
dataframes as well as dataframes with the `Python dataframe interchange protocol
<https://data-apis.org/dataframe-protocol/latest/index.html>`__ `__dataframe__`.
Furthermore, fitted estimators have the method :meth:`get_feature_names_out`. The
`get_feature_names_out` of a transformer returns⸺you guessed it⸺the feature names of
what `transform` returns.


Introducing the `set_output` API
================================

A further major step to support dataframes in a "dataframe in, dataframe out" fashion was
`SLEP018 <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep018/proposal.html>`__,
implemented for pandas dataframes in :ref:`version 1.2 <release_notes_1_2>` and for
polars dataframes in :ref:`version 1.4 <release_notes_1_4>`. It introduced the
`set_output` API to configure transformers to output pandas or polars DataFrames.
The output of transformers can be configured per estimator by calling
the :meth:`set_output` method or globally, by setting `set_config(transform_output="pandas")`.
Set it to `"polars"` instead of `"pandas"` if you want the same thing to happen but with
polars DataFrames.

The usage is basically as follows::

    >>> import numpy as np
    >>> import pandas as pd
    >>> from sklearn.compose import ColumnTransformer
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.preprocessing import OneHotEncoder
    >>> from sklearn.linear_model import LinearRegression

    >>> X = pd.DataFrame(
    ...     {"animals": ["cat", "cat", "dog", "dog"], "numeric": np.linspace(-1, 1, 4)}
    ... )
    >>> y = np.array([-1.5, 0, 0.1, 1.0])
    >>> ct = ColumnTransformer(
    ...     [("categorical", OneHotEncoder(sparse_output=False), ["animals"])],
    ...     remainder="passthrough",
    ... )
    >>> model = make_pipeline(ct, LinearRegression()).fit(X, y)
    >>> model.feature_names_in_
    array(['animals', 'numeric'], dtype=object)
    >>> model[0].get_feature_names_out()
    array(['categorical__animals_cat', 'categorical__animals_dog',
       'remainder__numeric'], dtype=object)
    >>> model[0].transform(X)
    array([[ 1.        ,  0.        , -1.        ],
           [ 1.        ,  0.        , -0.33333333],
           [ 0.        ,  1.        ,  0.33333333],
           [ 0.        ,  1.        ,  1.        ]])

Now the same, but with pandas set as output::

    >>> from sklearn import set_config
    >>> set_config(transform_output="pandas")
    >>> model[0].transform(X)
           c...

.. raw:: html

    <table border="1" class="dataframe">
        <thead>
            <tr style="text-align: right;">
                <th></th>
                <th>categorical__animals_cat</th>
                <th>categorical__animals_dog</th>
                <th>remainder__numeric</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <th>0</th>
                <td>1.0</td>
                <td>0.0</td>
                <td>-1.000000</td>
            </tr>
            <tr>
                <th>1</th>
                <td>1.0</td>
                <td>0.0</td>
                <td>-0.333333</td>
            </tr>
            <tr>
                <th>2</th>
                <td>0.0</td>
                <td>1.0</td>
                <td>0.333333</td>
            </tr>
            <tr>
                <th>3</th>
                <td>0.0</td>
                <td>1.0</td>
                <td>1.000000</td>
            </tr>
        </tbody>
    </table>

To return to the default, simply run::

    >>> set_config(transform_output="default")

A more detailed example can be found in
:ref:`sphx_glr_auto_examples_miscellaneous_plot_set_output.py`.
