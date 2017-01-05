.. _classifier_chain:

================
Classifier Chain
================

.. currentmodule:: sklearn.classifier_chain

Classifier chains are a way of combining a number of binary classifiers
into a single multi-label model that is capable of exploiting correlations
among labels.

For a multi-label classification problem with N classes, N binary
classifiers are assigned an integer between 0 and N-1. These integers
define the order of models in the chain. Each classifier is then fit on the
available training data plus the true labels of the classes whose
models were assigned a low order number.

When predicting the true labels will not be available. Instead the
predictions of each model are passed on to the subsequent models in the
chain to be used as features.

Clearly the order of the chain is important. The first model in the chain
has no information about the other labels while the last model in the chain
has features indicating the presence of all of the other labels. In general
one does not know the optimal ordering of the models in the chain so
typically many randomly ordered chains are fit and their predictions are
averaged together.

.. topic:: References:
    Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank,
        "Classifier Chains for Multi-label Classification", 2009.
