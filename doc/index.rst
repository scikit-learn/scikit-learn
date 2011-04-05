
..  
    We are putting the title as a raw HTML so that it doesn't appear in
    the contents

.. raw:: html

    <h1>scikits.learn: machine learning in Python</h1>

.. only:: html

    .. |banner1| image:: auto_examples/cluster/images/plot_affinity_propagation_1.png
       :height: 150
       :target: auto_examples/cluster/plot_affinity_propagation.html

    .. |banner2| image:: auto_examples/gaussian_process/images/plot_gp_regression_1.png
       :height: 150
       :target: auto_examples/gaussian_process/plot_gp_regression.html

    .. |banner3| image:: auto_examples/svm/images/plot_oneclass_1.png
       :height: 150
       :target: auto_examples/svm/plot_oneclass.html

    .. |banner4| image:: auto_examples/cluster/images/plot_lena_ward_segmentation_1.png
       :height: 150
       :target: auto_examples/cluster/plot_lena_ward_segmentation.html

    .. |center-div| raw:: html

        <div style="text-align: center; margin: 0px 0 -5px 0;">

    .. |end-div| raw:: html

        </div>


    |center-div| |banner1| |banner2| |banner3| |banner4| |end-div| 


.. topic:: Easy-to-use and general-purpose machine learning in Python

    ``scikits.learn`` is a Python module integrating classic machine
    learning algorithms in the tightly-knit world of scientific Python
    packages (`numpy <http://www.scipy.org>`_, `scipy
    <http://www.scipy.org>`_, `matplotlib
    <http://matplotlib.sourceforge.net/>`_).
    
    It aims to provide simple and efficient solutions to learning
    problems that are accessible to everybody and reusable in various
    contexts: **machine-learning as a versatile tool for science and
    engineering**.
    


:Features:
  * **Solid**: :ref:`supervised-learning`: :ref:`svm`, :ref:`linear_model`.

  * **Work in progress**: :ref:`unsupervised-learning`:
    :ref:`clustering`, :ref:`mixture`, manifold learning, :ref:`ICA
    <ICA>`, :ref:`gaussian_process`

  * **Planed**: Gaussian graphical models, matrix factorization

:License:
  Open source, commercially usable: **BSD license** (3 clause)


.. include:: big_toc_css.rst

.. note:: This document describes scikits.learn |release|. For other
   versions and printable format, see :ref:`documentation_resources`.

User Guide
==========

.. toctree::
   :maxdepth: 2

   contents

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
   developers/neighbors
   performance
   about
