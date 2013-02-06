.. include:: includes/big_toc_css.rst

.. _estimator_summary:

=========================
Estimators at a glance
=========================

Classification
--------------

.. list-table::
   :header-rows: 1

   * - Estimator
     - Dense
     - Sparse
     - Scalability
     - Multi-class
     - Multi-task

   * - :class:`.GaussianNB`
     - Yes
     - No
     -
     - Built-in
     - No

   * - :class:`.LinearSVC`
     - Yes
     - Yes
     -
     - One-Vs-Rest
     - No

   * - :class:`.LogisticRegression`
     - Yes
     - Yes
     -
     - One-Vs-Rest
     - No

   * - :class:`.MultinomialNB`
     - Yes
     - Yes
     -
     - Built-in
     - No

   * - :class:`.NuSVC`
     - Yes
     - Yes
     -
     - One-Vs-One
     - No

   * - :class:`.Perceptron`
     - Yes
     - Yes
     - Large n_samples
     - One-Vs-Rest
     - No

   * - :class:`.RidgeClassifier`
     - Yes
     - Yes
     -
     - One-Vs-Rest
     - No

   * - :class:`.RidgeClassifierCV`
     - Yes
     - Yes
     -
     - One-Vs-Rest
     - No

   * - :class:`.SGDClassifier`
     - Yes
     - Yes
     - Large n_samples
     - One-Vs-Rest
     - No

   * - :class:`.SVC`
     - Yes
     - Yes
     -
     - One-Vs-One
     - No

Regression
----------

.. list-table::
   :header-rows: 1

   * - Estimator
     - Dense
     - Sparse
     - Scalability
     - Multi-task

   * - :class:`.ARDRegression`
     - Yes
     - No
     -
     - No

   * - :class:`.BayesianRidge`
     - Yes
     - No
     -
     - No

   * - :class:`.ElasticNet`
     - Yes
     - No
     -
     - No

   * - :class:`.ElasticNetCV`
     - Yes
     - No
     -
     - No

   * - :class:`.Lars`
     - Yes
     - No
     -
     - No

   * - :class:`.LarsCV`
     - Yes
     - No
     -
     - No

   * - :class:`.Lasso`
     - Yes
     - No
     -
     - No

   * - :class:`.LassoCV`
     - Yes
     - No
     -
     - No

   * - :class:`.LassoLars`
     - Yes
     - No
     -
     - No

   * - :class:`.LassoLarsCV`
     - Yes
     - No
     -
     - No

   * - :class:`.LassoLarsIC`
     - Yes
     - No
     -
     - No

   * - :class:`.LinearRegression`
     - Yes
     - Yes
     -
     - Yes

   * - :class:`.NuSVR`
     - Yes
     - Yes
     -
     - No

   * - :class:`.OrthogonalMatchingPursuit`
     - Yes
     - No
     -
     - Yes

   * - :class:`.Ridge`
     - Yes
     - Yes
     -
     - Yes

   * - :class:`.RidgeCV`
     - Yes
     - Yes
     -
     - Yes

   * - :class:`.SGDRegressor`
     - Yes
     - Yes
     -
     - Yes

   * - :class:`.SVR`
     - Yes
     - Yes
     -
     - No


Clustering
-----------

.. list-table::
   :header-rows: 1

   * - Estimator
     - Dense
     - Sparse
     - Scalability

   * - :class:`.AffinityPropagation`
     - Yes
     - No
     -

   * - :class:`.DBSCAN`
     - Yes
     - No
     -

   * - :class:`.KMeans`
     - Yes
     - Yes
     -

   * - :class:`.MiniBatchKMeans`
     - Yes
     - Yes
     -

   * - :class:`.MeanShift`
     - Yes
     - No
     -

   * - :class:`.SpectralClustering`
     - Yes
     - No
     -

   * - :class:`.Ward`
     - Yes
     - No
     -
