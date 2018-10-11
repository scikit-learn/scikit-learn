.. _fast_kernel:

================================================================
Fast Kernel Machine (EigenPro) for Regression and Classification
================================================================

.. currentmodule:: sklearn.fast_kernel

Fast Kernel Machine is a very efficient implementation of kernel regression/classification
using *EigenPro iteration* [MB17]_,
an optimization method based on preconditioned stochastic gradient descent.
It essentially implements a "ridgeless" kernel regression.
Regularization, when necessary, can be achieved by early stopping.

Optimization parameters, such as step size, batch size, and the size of the preconditioning
block are chosen automatically and optimally. (They can also be set up manually.)
This results in a simple and user-friendly interface.

The figure below compares the Fast Kernel Classifier (EigenPro) and the Support Vector Classifier (:class:`SVC`) 
on MNIST digits classification task.
We see that EigenPro and SVC give competitive and similar accuracy on test set.
Notably, on the full MNIST training and testing using EigenPro are
approximately 100 times faster than that using SVC.

.. |mnist_acc| image:: ../images/fast_kernel_mnist_accuracy.png
    :target: ../auto_examples/fast_kernel/plot_mnist.html
    :scale: 50

.. |mnist_time| image:: ../images/fast_kernel_mnist_time.png
    :target: ../auto_examples/fast_kernel/plot_mnist.html
    :scale: 50

.. centered:: |mnist_acc| |mnist_time|


The next figure compares EigenPro and SVC on a binary classification problem with synthetic features.
Again, EigenPro demonstrates nearly 100 times acceleration on training and testing without loss of accuracy.

.. |synth_acc| image:: ../images/fast_kernel_classification_accuracy.png
    :target: ../auto_examples/fast_kernel/plot_classification.html
    :scale: 50

.. |synth_time| image:: ../images/fast_kernel_classification_times.png
    :target: ../auto_examples/fast_kernel/plot_classification.html
    :scale: 50

.. centered:: |synth_acc| |synth_time|


.. topic:: References:

    .. [MB17] Siyuan Ma and Mikhail Belkin,
       `"Diving into the shallows: a computational perspective on large-scale shallow learning"
       <https://arxiv.org/abs/1703.10622>`_,
       Advances in Neural Information Processing Systems, 2017.
