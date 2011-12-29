
..
    We are putting the title as a raw HTML so that it doesn't appear in
    the contents

.. raw:: html

    <h1>scikit-learn: machine learning in Python</h1>
    <style type="text/css">
    p {
        margin: 7px 0 7px 0 ;
    }
    span.linkdescr a {
        color:  #3E4349 ;
    }
    </style>

.. only:: html

    .. |banner1| image:: auto_examples/svm/images/plot_oneclass_1.png
       :height: 150
       :target: auto_examples/svm/plot_oneclass.html

    .. |banner2| image:: auto_examples/cluster/images/plot_ward_structured_vs_unstructured_2.png
       :height: 150
       :target: auto_examples/cluster/plot_ward_structured_vs_unstructured.html

    .. |banner3| image:: auto_examples/gaussian_process/images/plot_gp_regression_1.png
       :height: 150
       :target: auto_examples/gaussian_process/plot_gp_regression.html

    .. |banner4| image:: auto_examples/cluster/images/plot_lena_ward_segmentation_1.png
       :height: 150
       :target: auto_examples/cluster/plot_lena_ward_segmentation.html

    .. |center-div| raw:: html

        <div style="text-align: center; margin: -7px 0 -13px 0;">

    .. |end-div| raw:: html

        </div>


    |center-div| |banner1| |banner2| |banner3| |banner4| |end-div|


.. topic:: Easy-to-use and general-purpose machine learning in Python

    ``scikit-learn`` is a Python module integrating classic machine
    learning algorithms in the tightly-knit world of scientific Python
    packages (`numpy <http://numpy.scipy.org>`_, `scipy
    <http://www.scipy.org>`_, `matplotlib
    <http://matplotlib.sourceforge.net/>`_).
    It aims to provide simple and efficient solutions to learning
    problems that are accessible to everybody and reusable in various
    contexts: **machine-learning as a versatile tool for science and
    engineering**.


.. raw:: html

  <table class="contentstable" style="width: 100% ; margin-top: -8px">
    <tr valign="top"><td width="28%">
      <p class="biglink"><a class="biglink" href="supervised_learning.html">
                Supervised learning</a><br/>
         <span class="linkdescr">
                <a href="modules/svm.html">Support vector machines</a>,
                <a href="modules/linear_model.html">linear models</a>,
                <a href="modules/naive_bayes.html">naives Bayes</a>,
                <a href="modules/gaussian_process.html">Gaussian process</a>...
         </span></p>
    </td><td align="center" width="32%">
      <p class="biglink"><a class="biglink" href="unsupervised_learning.html">
        Unsupervised learning</a><br/>
         <span class="linkdescr">
                <a href="modules/clustering.html">Clustering</a>,
                <a href="modules/mixture.html">Gaussian mixture models</a>,
                <a href="modules/manifold.html">manifold learning</a>,
                <a href="modules/decomposition.html">matrix factorization</a>,
                <a href="modules/covariance.html">covariance</a>...
         </span></p>
    </td><td align="right" width="30%">
      <p class="biglink"><a class="biglink" href="index.html#user-guide">
        And much more</a><br/>
         <span class="linkdescr">
                <a href="model_selection.html">Model selection</a>,
                <a href="datasets/index.html">datasets</a>,
                <a href="modules/feature_extraction.html">feature extraction...</a>
                <strong>See below</strong>.</span></p>
    </td></tr>
  </table>

**License:** Open source, commercially usable: **BSD license** (3 clause)

.. include:: includes/big_toc_css.rst

Documentation for scikit-learn **version** |release|. For other versions and
printable format, see :ref:`documentation_resources`.

User Guide
==========

.. toctree::
   :maxdepth: 2

   user_guide.rst

Example Gallery
===============

.. toctree::
   :maxdepth: 2

   auto_examples/index


Development
===========
.. toctree::
   :maxdepth: 2

   developers/index
   developers/performance
   developers/utilities
   about
