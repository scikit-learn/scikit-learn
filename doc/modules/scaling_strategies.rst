.. _scaling_strategies:

==================
Scaling Strategies
==================

For some applications the amount of examples, features (or both) and/or the 
speed at which they need to be processed are challenging for traditional 
approaches. In these cases Scikit-learn has a number of options you can 
consider to make your system scale. 

Scaling with #examples using out-of-core learning
=================================================

Out-of-core (or "external memory") learning is a technique used to learn from
data that cannot fit in a computer's main memory (RAM). 

Here is sketch of a system designed to achieve this goal:
1. a way to stream examples
2. a way to extract features from examples
3. an online algorithm

Basically, \1. could be either a reader that `yield`s examples from files on a
hard drive or a program that listens to a socket and parses examples from a 
network stream. This is the case of the *all-reduce* computing paradigm where
all mappers of a computing cluster send their data to one learner as 
popularized by projects as `Vowpal wabbit
<https://github.com/JohnLangford/vowpal_wabbit/wiki>`_. However, details on 
how to achieve this are beyond the scope of this documentation.

\2. could be any relevant way to extract features among the 
different :ref:`feature extraction <feature_extraction>` methods supported by
Scikit-learn. However, when working with data that needs vectorization and 
where the set of features or values is not known in advance one should take 
explicit care by using a stateless feature extractor. A good example is text
classification where unknown terms are likely to be found during training.
Currently the preferred way to do this is to use the so-called 
:ref:`hashing trick<feature_hashing>` as implemented by 
:class:`sklearn.feature_extraction.FeatureHasher`.

Finally, for \3. we have a number of options inside Scikit-learn. Although all
learning algorithms can not learn online (i.e. without seeing all the examples
at once), all estimators implementing the `partial_fit` API are candidates.

Here is list of online estimators for different tasks:
- Classification
  + :class:`Perceptron`
  + :class:`MultinomialNB`
  + :class:`BernoulliNB`
  + :class:`SGDClassifier`
  + :class:`PassiveAggressiveClassifier`
- Regression
  + :class:`SGDRegressor`
  + :class:`PassiveAggressiveRegressor`
- Clustering
  + :class:`MinibatchKMeans`
- Decomposition / feature Extraction
  + :class:`MinibatchDictionaryLearning`
  + :class:`MinibatchKMeans`

An somewhat important thing to note is that although a stateless feature 
extraction routine may be able to to cope with new/unseen attributes, the 
online learner itself is often unable to cope with new/unseen classes. E.g. for
classification you need to pass all the possible classes to the first 
`partial_fit` call using the `classes=` parameter.







