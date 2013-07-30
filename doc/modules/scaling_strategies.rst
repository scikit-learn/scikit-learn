.. _scaling_strategies:

==================
Scaling Strategies
==================

For some applications the amount of examples, features (or both) and/or the 
speed at which they need to be processed are challenging for traditional 
approaches. In these cases Scikit-learn has a number of options you can 
consider to make your system scale. 

Scaling with instances using out-of-core learning
=================================================

Out-of-core (or "external memory") learning is a technique used to learn from
data that cannot fit in a computer's main memory (RAM). 

Here is sketch of a system designed to achieve this goal:

  1. a way to stream instances
  2. a way to extract features from instances
  3. an online algorithm

Streaming instances
------------------
Basically, \1. could be either a reader that yields instances from files on a
hard drive or a program that listens to a socket and parses instances from a 
network stream. This is the case of the *all-reduce* computing paradigm where
all mappers of a computing cluster send their data to one learner as 
popularized by projects like 
`Vowpal Wabbit <https://github.com/JohnLangford/vowpal_wabbit/wiki>`_. However, 
details on how to achieve this are beyond the scope of this documentation.

Extracting features
-------------------
\2. could be any relevant way to extract features among the 
different :ref:`feature extraction <feature_extraction>` methods supported by
Scikit-learn. However, when working with data that needs vectorization and 
where the set of features or values is not known in advance one should take 
explicit care by using a stateless feature extractor. A good example is text
classification where unknown terms are likely to be found during training.
Currently the preferred way to do this is to use the so-called 
:ref:`hashing trick<feature_hashing>` as implemented by 
:class:`sklearn.feature_extraction.FeatureHasher`.

Online learning
---------------
Finally, for \3. we have a number of options inside Scikit-learn. Although all
learning algorithms can not learn online (i.e. without seeing all the instances
at once), all estimators implementing the `partial_fit` API are candidates.
Actually, the ability to learn iteratively from a mini-batch of instances is key
to out-of-core learning as it guarantees that at any given time there will be
only a small amount of instances in the main memory. Choosing a good size for
the mini-batch that balances relevancy and memory footprint could involve some
tuning.

Here is a list of online estimators for different tasks:

  - Classification
      + :class:`sklearn.naive_bayes.MultinomialNB`
      + :class:`sklearn.naive_bayes.BernoulliNB`
      + :class:`sklearn.linear_model.Perceptron`
      + :class:`sklearn.linear_model.SGDClassifier`
      + :class:`sklearn.linear_model.PassiveAggressiveClassifier`
  - Regression
      + :class:`sklearn.linear_model.SGDRegressor`
      + :class:`sklearn.linear_model.PassiveAggressiveRegressor`
  - Clustering
      + :class:`sklearn.cluster.MiniBatchKMeans`
  - Decomposition / feature Extraction
      + :class:`sklearn.decomposition.MiniBatchDictionaryLearning`
      + :class:`sklearn.cluster.MiniBatchKMeans`

An somewhat important thing to note is that although a stateless feature 
extraction routine may be able to to cope with new/unseen attributes, the 
online learner itself may be unable to cope with new/unseen targets. E.g. for
classification you need to pass all the possible classes to the first 
`partial_fit` call using the `classes=` parameter.

Examples
--------
Finally, we have a full-fledged example of
:ref:`example_applications_plot_out_of_core_classification.py` comparing the
performance of different algorithms with the number of processed examples.




