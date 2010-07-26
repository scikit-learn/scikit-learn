.. raw:: html

  <style type="text/css">
    li.toctree-l1 {
        padding: 0.5em 0 1em 0 ;
        list-style-type: none;
        font-size: 150% ;
        }

    li.toctree-l2 {
        font-size: 70% ;
        list-style-type: square;
        }

    li.toctree-l3 {
        font-size: 85% ;
        list-style-type: circle;
        }

    div.bodywrapper h1 {
        text-align: center;
        font-size: 300% ;
    }
  
  </style>



===========================================
Scikits.learn: machine learning in Python
===========================================

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
    
.. image:: auto_examples/images/plot_digits_classification.png
    :align: right
    :scale: 50


.. raw:: html

    <small>

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
    
.. raw:: html

    </small>


**Features:**

 * **Solid**: supervised learning: classification, regression

 * **Work in progress**: unsupervised learning: clustering, mixture modeling,
   manifold learning

 * **Planed**: Gaussian graphical models, matrix factorization, ICA

Download
========

Click `here <https://sourceforge.net/projects/scikit-learn/files/>`__
to download latest release. Previous releases can also be found in
that directory.

Mailing List
============

Visit `this page
<https://lists.sourceforge.net/lists/listinfo/scikit-learn-general>`_
to subscribe to the mailing list and keep informed about scikit-learn
development


User guide
======================

.. warning:: 

   This documentation is relative to the development version,
   documentation for the stable version can be found `here
   <http://scikit-learn.sourceforge.net/old_doc/>`__

.. toctree::
   :maxdepth: 3

   install
   tutorial
   supervised_learning
   unsupervised_learning
   auto_examples/index
   contribute
   .. API
