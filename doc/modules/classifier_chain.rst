.. _classifier_chain:

===========================
Classifier Chain
===========================

.. currentmodule:: sklearn.classifier_chain

Classifier chains are a way of combining a number of binary classifiers
into a single multi-label model that is capable of exploiting
correlations among labels.

The fitting procedure for a classifier chain is as follows. For a multi-label
classification problem with N classes, N binary classifiers are arranged
into a "chain". Each classifier in the chain is fit using all of the
features in the training data plus the true labels of the classes whose
models appear earlier in the chain.

The first model in the chain is fit on just the features in the training
data. The second model in the chain is fit on the features in the training
data plus an additional binary feature indicating the presence of the first
label. The features for the third model are augmented by two additional
features indicating the presence of the first two labels and so on. When
predicting the true labels will not be available. Instead the predictions of
 each model are passed on to the subsequent models in the chain to be used
 as features. The additional features provided to each model by earlier
 models in the chain makes it possible to exploit correlations among the
 labels.

 Clearly the order of the chain is important. The first model in the chain
 has no information about the other labels while the last model in the chain
  has features indicating the presence of all of the other labels. In
  general one does not know the optimal ordering of the models in the chain
  so typically many randomly ordered chains are fit and their predictions
  are averaged together.


.. topic:: References:
    Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank,
        "Classifier Chains for Multi-label Classification", 2009.
