
.. _adult_census:

Adult census
============

The samples in this dataset was extracted by Barry Becker from the 1994 Census
database. It consisted of 48842 rows.
A set of reasonably clean records was extracted using the following
conditions: ((AGE>16) && (AGI>100) && (AFNLWGT>1)&& (HRSWK>0)).
Prediction task is to determine whether a person makes over 50K a year.

Each sample has 15 features, described on the
`dataset's homepage <https://archive.ics.uci.edu/ml/datasets/Adult>`_.
Some of the features are integers,
while others are categoricals.

:func:`sklearn.datasets.fetch_adult_census` will load the adult census dataset;
it returns a dictionary-like object with the feature matrix in the ``data``
member, list of feature names in `feature_names`, dataset description in `DESCR`
and the target values in ``target``.
The dataset will be downloaded from the web if necessary.
