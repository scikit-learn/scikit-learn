.. include:: includes/big_toc_css.rst

.. _inspection:

Inspection
----------

Predictive performance is often the main goal of developing machine learning
models. Yet summarising performance with an evaluation metric is often
insufficient: it assumes that the evaluation metric and test dataset
perfectly reflect the target domain, which is rarely true. In certain domains,
a model needs a certain level of interpretability before it can be deployed.
A model that is exhibiting performance issues needs to be debugged for one to 
understand the model's underlying issue. The 
:mod:`sklearn.inspection` module provides tools to help understand the 
predictions from a model and what affects them. This can be used to 
evaluate assumptions and biases of a model, design a better model, or
to diagnose issues with model performance.

.. toctree::

    modules/partial_dependence
    modules/permutation_importance
