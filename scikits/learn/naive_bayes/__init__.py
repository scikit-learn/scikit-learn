"""
Naive Bayes models
==================

Naive Bayes algorithms are a set of supervised learning methods based on
applying Bayes' theorem with strong (naive) independence assumptions. 

See http://scikit-learn.sourceforge.net/modules/naive_bayes.html for
complete documentation.
"""

from .naive_bayes import GNB, MultinomialNB

from . import sparse

