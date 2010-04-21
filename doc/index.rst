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



===========================
Scikits.learn Documentation
===========================

Introduction
============
``scikits.learn`` is a Python module for machine learning. It aims to
implement classic machine learning algorithms while remaining simple
and efficient.

It implements several machine learning algorithms, including Support
Vector Machines, Gaussian Mixture Models, Neural Networks, Nearest
Neighbors, Generalized Linear Models, etc.

.. image:: auto_examples/images/plot_digits_classification.png
    :align: right
    :scale: 50


.. raw:: html

    <small>

**Example**::

    import pylab as pl

    from scikits.learn import datasets, svm
    digits = datasets.load_digits()
    for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
        pl.subplot(2, 4, index+1)
        pl.imshow(image, cmap=pl.cm.gray_r)
        pl.title('Training: %i' % label)
    
    n_features = len(digits.images)
    data = digits.images.reshape((n_features, -1))
    
    classifier = svm.SVC()
    classifier.fit(data[:n_features/2], digits.target[:n_features/2])
    
    for index, image in enumerate(digits.images[n_features/2:n_features/2+4]):
        pl.subplot(2, 4, index+5)
        pl.imshow(image, cmap=pl.cm.gray_r)
        pl.title('Prediction: %i' % classifier.predict(image.ravel()))
    
.. raw:: html

    </small>

    
User guide: contents
======================

.. toctree::
   :maxdepth: 2

   install
   tutorial
   module/index
   auto_examples/index
   contribute
   .. API
