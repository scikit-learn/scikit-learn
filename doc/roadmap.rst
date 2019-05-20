.. _roadmap:

Roadmap
=======

Purpose of this document
------------------------
This document list general directions that core contributors are interested
to see developed in scikit-learn. The fact that an item is listed here is in
no way a promise that it will happen, as resources are limited. Rather, it
is an indication that help is welcomed on this topic.

Statement of purpose: Scikit-learn in 2018
------------------------------------------
Eleven years after the inception of Scikit-learn, much has changed in the
world of machine learning. Key changes include:

* Computational tools: The exploitation of GPUs, distributed programming
  frameworks like Scala/Spark, etc.
* High-level Python libraries for experimentation, processing and data
  management: Jupyter notebook, Cython, Pandas, Dask, Numba...
* Changes in the focus of machine learning research: artificial intelligence
  applications (where input structure is key) with deep learning,
  representation learning, reinforcement learning, domain transfer, etc.

A more subtle change over the last decade is that, due to changing interests
in ML, PhD students in machine learning are more likely to contribute to
PyTorch, Dask, etc. than to Scikit-learn, so our contributor pool is very
different to a decade ago.

Scikit-learn remains very popular in practice for trying out canonical
machine learning techniques, particularly for applications in experimental
science and in data science. A lot of what we provide is now very mature.
But it can be costly to maintain, and we cannot therefore include arbitrary
new implementations. Yet Scikit-learn is also essential in defining an API
framework for the development of interoperable machine learning components
external to the core library.

**Thus our main goals in this era are to**:

* continue maintaining a high-quality, well-documented collection of canonical
  tools for data processing and machine learning within the current scope
  (i.e. rectangular data largely invariant to column and row order;
  predicting targets with simple structure)
* improve the ease for users to develop and publish external components
* improve inter-operability with modern data science tools (e.g. Pandas, Dask)
  and infrastructures (e.g. distributed processing)

Many of the more fine-grained goals can be found under the `API tag
<https://github.com/scikit-learn/scikit-learn/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc+label%3AAPI>`_
on the issue tracker.

Architectural / general goals
-----------------------------
The list is numbered not as an indication of the order of priority, but to
make referring to specific points easier. Please add new entries only at the
bottom.

#. Everything in Scikit-learn should conform to our API contract

   * `Pipeline <pipeline.Pipeline>` and `FeatureUnion` modify their input
     parameters in fit. Fixing this requires making sure we have a good
     grasp of their use cases to make sure all current functionality is
     maintained. :issue:`8157` :issue:`7382`

#. Improved handling of Pandas DataFrames and SparseDataFrames

   * document current handling
   * column reordering issue :issue:`7242`
   * avoiding unnecessary conversion to ndarray
   * returning DataFrames from transformers :issue:`5523`
   * getting DataFrames from dataset loaders
   * Sparse currently not considered

#. Improved handling of categorical features

   * Tree-based models should be able to handle both continuous and categorical
     features :issue:`4899`
   * In dataset loaders
   * As generic transformers to be used with ColumnTransforms (e.g. ordinal
     encoding supervised by correlation with target variable)

#. Improved handling of missing data

   * Making sure meta-estimators are lenient towards missing data
   * Non-trivial imputers
   * Learners directly handling missing data
   * An amputation sample generator to make parts of a dataset go missing
   * Handling mixtures of categorical and continuous variables

#. Passing around information that is not (X, y): Sample properties

   * We need to be able to pass sample weights to scorers in cross validation.
   * We should have standard/generalised ways of passing sample-wise properties
     around in meta-estimators. :issue:`4497` :issue:`7646`

#. Passing around information that is not (X, y): Feature properties

   * Feature names or descriptions should ideally be available to fit for, e.g.
     . :issue:`6425` :issue:`6424`
   * Per-feature handling (e.g. "is this a nominal / ordinal / English language
     text?") should also not need to be provided to estimator constructors,
     ideally, but should be available as metadata alongside X. :issue:`8480`

#. Passing around information that is not (X, y): Target information

   * We have problems getting the full set of classes to all components when
     the data is split/sampled. :issue:`6231` :issue:`8100`
   * We have no way to handle a mixture of categorical and continuous targets.

#. Make it easier for external users to write Scikit-learn-compatible
   components

   * More flexible estimator checks that do not select by estimator name
     :issue:`6599` :issue:`6715`
   * Example of how to develop a meta-estimator
   * More self-sufficient running of scikit-learn-contrib or a similar resource

#. Support resampling and sample reduction

   * Allow subsampling of majority classes (in a pipeline?) :issue:`3855`
   * Implement random forests with resampling :issue:`8732`

#. Better interfaces for interactive development

   * __repr__ and HTML visualisations of estimators :issue:`6323`
   * Include plotting tools, not just as examples. :issue:`9173`

#. Improved tools for model diagnostics and basic inference

   * alternative feature importances implementations (e.g. methods or wrappers)
   * better ways to handle validation sets when fitting
   * better ways to find thresholds / create decision rules :issue:`8614`

#. Better tools for selecting hyperparameters with transductive estimators

   * Grid search and cross validation are not applicable to most clustering
     tasks. Stability-based selection is more relevant.

#. Improved tracking of fitting

   * Verbose is not very friendly and should use a standard logging library
     :issue:`6929`
   * Callbacks or a similar system would facilitate logging and early stopping

#. Distributed parallelism

   * Joblib can now plug onto several backends, some of them can distribute the
     computation across computers
   * However, we want to stay high level in scikit-learn

#. A way forward for more out of core

   * Dask enables easy out-of-core computation. While the dask model probably
     cannot be adaptable to all machine-learning algorithms, most machine
     learning is on smaller data than ETL, hence we can maybe adapt to very
     large scale while supporting only a fraction of the patterns.

#. Better support for manual and automatic pipeline building

   * Easier way to construct complex pipelines and valid search spaces
     :issue:`7608` :issue:`5082` :issue:`8243`
   * provide search ranges for common estimators??
   * cf. `searchgrid <https://searchgrid.readthedocs.io/en/latest/>`_

#. Support for working with pre-trained models

   * Estimator "freezing". In particular, right now it's impossible to clone a
     `CalibratedClassifierCV` with prefit. :issue:`8370`. :issue:`6451`

#. Backwards-compatible de/serialization of some estimators

   * Currently serialization (with pickle) breaks across versions. While we may
     not be able to get around other limitations of pickle re security etc, it
     would be great to offer cross-version safety from version 1.0. Note: Gael
     and Olivier think that this can cause heavy maintenance burden and we
     should manage the trade-offs. A possible alternative is presented in the
     following point.

#. Documentation and tooling for model lifecycle management

   * Document good practices for model deployments and lifecycle: before
     deploying a model: snapshot the code versions (numpy, scipy, scikit-learn,
     custom code repo), the training script and an alias on how to retrieve
     historical training data + snapshot a copy of a small validation set +
     snapshot of the predictions (predicted probabilities for classifiers)
     on that validation set.
   * Document and tools to make it easy to manage upgrade of scikit-learn
     versions:

     * Try to load the old pickle, if it works, use the validation set
       prediction snapshot to detect that the serialized model still behave
       the same;
     * If joblib.load / pickle.load not work, use the versioned control
       training script + historical training set to retrain the model and use
       the validation set prediction snapshot to assert that it is possible to
       recover the previous predictive performance: if this is not the case
       there is probably a bug in scikit-learn that needs to be reported.

#. (Optional) Improve scikit-learn common tests suite to make sure that (at
   least for frequently used) models have stable predictions across-versions
   (to be discussed);

   * Extend documentation to mention how to deploy models in Python-free
     environments for instance  `ONNX <https://github.com/onnx/onnxmltools>`_.
     and use the above best practices to assess predictive consistency between
     scikit-learn and ONNX prediction functions on validation set.
   * Document good practices to detect temporal distribution drift for deployed
     model and good practices for re-training on fresh data without causing
     catastrophic predictive performance regressions.

#. More didactic documentation

   * More and more options have been added to scikit-learn. As a result, the
     documentation is crowded which makes it hard for beginners to get the big
     picture. Some work could be done in prioritizing the information.

Subpackage-specific goals
-------------------------

:mod:`sklearn.cluster`

* kmeans variants for non-Euclidean distances, if we can show these have
  benefits beyond hierarchical clustering.

:mod:`sklearn.ensemble`

* a stacking implementation

:mod:`sklearn.model_selection`

* multi-metric scoring is slow :issue:`9326`
* perhaps we want to be able to get back more than multiple metrics
* the handling of random states in CV splitters is a poor design and
  contradicts the validation of similar parameters in estimators.
* exploit warm-starting and path algorithms so the benefits of `EstimatorCV`
  objects can be accessed via `GridSearchCV` and used in Pipelines.
  :issue:`1626`
* Cross-validation should be able to be replaced by OOB estimates whenever a
  cross-validation iterator is used.
* Redundant computations in pipelines should be avoided (related to point
  above) cf `daskml
  <https://dask-ml.readthedocs.io/en/latest/hyper-parameter-search.html#avoid-repeated-work>`_

:mod:`sklearn.neighbors`

* Ability to substitute a custom/approximate/precomputed nearest neighbors
  implementation for ours in all/most contexts that nearest neighbors are used
  for learning. :issue:`10463`

:mod:`sklearn.pipeline`

* Performance issues with `Pipeline.memory`
* see "Everything in Scikit-learn should conform to our API contract" above
