=========================
Putting it all together
=========================

..  Imports
    >>> import numpy as np

Pipelining
============

We have seen that some estimators can transform data and that some estimators
can predict variables. We can also create combined estimators:

.. image:: ../../auto_examples/images/plot_digits_pipe_1.png
   :target: ../../auto_examples/plot_digits_pipe.html
   :scale: 65
   :align: right

.. literalinclude:: ../../auto_examples/plot_digits_pipe.py
    :lines: 24-67




Face recognition with eigenfaces
=================================

The dataset used in this example is a preprocessed excerpt of the
"Labeled Faces in the Wild", also known as LFW_:

  http://vis-www.cs.umass.edu/lfw/lfw-funneled.tgz (233MB)

.. _LFW: http://vis-www.cs.umass.edu/lfw/

.. literalinclude:: ../../auto_examples/applications/face_recognition.py

.. |prediction| image:: ../../images/plot_face_recognition_1.png
   :scale: 50
 
.. |eigenfaces| image:: ../../images/plot_face_recognition_2.png
   :scale: 50

.. list-table::
   :class: centered

   *

     - |prediction|

     - |eigenfaces|

   * 

     - **Prediction**

     - **Eigenfaces**

Expected results for the top 5 most represented people in the dataset::

                     precision    recall  f1-score   support

  Gerhard_Schroeder       0.91      0.75      0.82        28
    Donald_Rumsfeld       0.84      0.82      0.83        33
         Tony_Blair       0.65      0.82      0.73        34
       Colin_Powell       0.78      0.88      0.83        58
      George_W_Bush       0.93      0.86      0.90       129

        avg / total       0.86      0.84      0.85       282


Open problem: Stock Market Structure
=====================================

Can we predict the variation in stock prices for Google over a given time frame?

:ref:`stock_market`
