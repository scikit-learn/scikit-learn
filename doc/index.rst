
===========================================
scikits.learn: machine learning in python
===========================================

.. |banner1| image:: auto_examples/cluster/images/plot_affinity_propagation.png
   :height: 150
   :target: auto_examples/cluster/plot_affinity_propagation.html


.. |banner2| image:: auto_examples/glm/images/plot_lasso_lars.png
   :height: 150
   :target: auto_examples/svm/plot_custom_kernel.html

.. |banner3| image:: auto_examples/svm/images/plot_oneclass.png
   :height: 150
   :target: auto_examples/svm/plot_oneclass.html

.. |banner4| image:: auto_examples/cluster/images/plot_lena_segmentation.png
   :height: 150
   :target: auto_examples/cluster/plot_lena_segmentation.html

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
  * **Solid**: :ref:`supervised-learning`: classification, regression

  * **Work in progress**: :ref:`unsupervised-learning`: :ref:`clustering`, 
    :ref:`gmm`, manifold learning, ICA

  * **Planed**: Gaussian graphical models, matrix factorization

:License:
  Open source, commercially usable: **BSD license** (3 clause)


.. raw:: html

   <div class="example_digits">

:ref:`A simple Example: recognizing hand-written digits <example_plot_digits_classification.py>` ::

    import pylab as pl

    from scikits.learn import datasets, svm
    digits = datasets.load_digits()
    for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
        pl.subplot(2, 4, index+1)
        pl.imshow(image, cmap=pl.cm.gray_r)
        pl.title('Training: %i' % label)
    
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))
    
    classifier = svm.SVC()
    classifier.fit(data[:n_samples/2], digits.target[:n_samples/2])
    
    for index, image in enumerate(digits.images[n_samples/2:n_samples/2+4]):
        pl.subplot(2, 4, index+5)
        pl.imshow(image, cmap=pl.cm.gray_r)
        pl.title('Prediction: %i' % classifier.predict(image.ravel()))

.. image:: images/plot_digits_classification.png
   :height: 140
   :target: auto_examples/plot_digits_classification.html

.. raw:: html

    </div>


User guide
======================

.. warning:: 

   This documentation is relative to the development version,
   documentation for the stable version can be found `here
   <http://scikit-learn.sourceforge.net/old_doc/>`__

.. toctree::
   :maxdepth: 3

   install

Tutorial
========

.. toctree::
   :maxdepth: 2

   tutorial

Reference
=========

.. toctree::
   :maxdepth: 2

   supervised_learning
   unsupervised_learning
   model_selection
   cross_validation
   modules/classes

Gallery
=======

.. toctree::
   :maxdepth: 2

   auto_examples/index

Developement
============
.. toctree::
   :maxdepth: 2

   developers/index
   performance
