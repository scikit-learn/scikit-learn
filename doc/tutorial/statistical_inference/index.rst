.. _stat_learn_tut_index:

==========================================================================
A tutorial on statistical-learning for scientific data processing
==========================================================================

.. topic:: Statistical learning 

    `Machine learning <http://en.wikipedia.org/wiki/Machine_learning>`_ is 
    a technique with a growing importance, as the
    size of the datasets experimental sciences are facing is rapidly
    growing. Problems it tackles range from building a prediction function
    linking different observations, to classifying observations, or
    learning the structure in an unlabeled dataset. 
    
    This tutorial will explore `statistical learning`, that is the use of
    machine learning techniques with the goal of `statistical inference 
    <http://en.wikipedia.org/wiki/Statistical_inference>`_:
    drawing conclusions on the data at hand.

    ``sklearn`` is a Python module integrating classic machine
    learning algorithms in the tightly-knit world of scientific Python
    packages (`numpy <http://www.scipy.org>`_, `scipy
    <http://www.scipy.org>`_, `matplotlib
    <http://matplotlib.sourceforge.net/>`_).

.. include:: ../../includes/big_toc_css.rst

.. warning::

    In scikit-learn release 0.9, the import path has changed from
    `scikits.learn` to `sklearn`. To import with cross-version 
    compatibility, use::

        try:
            from sklearn import something
        except ImportError:
            from scikits.learn import something


.. toctree::
   :maxdepth: 2

   settings
   supervised_learning
   model_selection
   unsupervised_learning
   putting_together
   finding_help

