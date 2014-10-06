.. _neural_network:

==================================
Neural network models (supervised)
==================================

.. currentmodule:: sklearn.neural_network


.. _multilayer_perceptron:
 
Multi-layer Perceptron
======================

**Multi-layer Perceptron (MLP)** is a supervised learning algorithm that learns
a function :math:`f(\cdot): R \rightarrow R` by training on a dataset. Given a 
set of features :math:`X = {x_1, x_2, ..., x_m}` and a target :math:`y`, it 
can learn a non-linear, complex functions for either classification
or regression. It is different from logistic regression, in that between the input 
and the output layer, there can be one or more hidden layers. Figure 1 shows 
a one hidden layer MLP.

.. figure:: ../images/multilayerperceptron_network.png
   :align: center
   :scale: 60%

   **Figure 1 : One hidden layer MLP.**

The leftmost layer, known as the input layer, consists of a set of neurons 
:math:`\{x_i | x_1, x_2, ..., x_m\}` representing the input features. Each hidden
layer transforms the values from the previous layer by a weighted linear summation
:math:`w_1x_1 + w_2x_2 + ... + w_mx_m`, followed by a non-linear activation function
:math:`g(\cdot):R \rightarrow R` - like the hyperbolic tan function. The output layer
receives the values from the last hidden layer and transforms them into output values.

The module contains the public attributes ``layers_coef_`` and ``layers_intercept_``.
``layers_coef_`` is a list of weight matrices, where weight matrix at index
:math:`i` represents the weights between layer :math:`i` and layer 
:math:`i+1`. ``layers_intercept_`` is a list of bias vectors, where the vector
at index :math:`i` represents the bias values added to layer :math:`i+1`.

The advantages of Multi-layer Perceptron are:

    + Capability to learn complex/non-linear models.

    + Capability to learn models in real-time (on-line learning) 
      using ``partial_fit``.
      
      
The disadvantages of Multi-layer Perceptron (MLP) include:

    + MLP with hidden layers have a non-convex loss function where there exists 
      more than one local minimum. Therefore different random weight 
      initializations can lead to different validation accuracy.

    + MLP suffers from the Backpropagation diffusion problem; layers far from 
      the output update with decreasing momentum, leading to slow convergence.

    + MLP requires tuning a number of hyperparameters such as the number of 
      hidden neurons, layers, and iterations.

    + MLP is sensitive to feature scaling.

Please see :ref:`Tips on Practical Use <mlp_tips>` section that addresses 
some of these disadvantages.


Classification
==============

Class :class:`MultilayerPerceptronClassifier` implements  
a multi layer perceptron (MLP) algorithm that trains using Backpropagation. 

MLP trains on two arrays: array X of size (n_samples, n_features), which holds 
the training samples represented as floating point feature vectors; and array 
y of size (n_samples,), which holds the target values (class labels) for the 
training samples::

    >>> from sklearn.neural_network import MultilayerPerceptronClassifier
    >>> X = [[0., 0.], [1., 1.]]
    >>> y = [0, 1]
    >>> clf = MultilayerPerceptronClassifier(hidden_layer_sizes=(5, 2), random_state=1)
    >>> clf.fit(X, y)
    MultilayerPerceptronClassifier(activation='relu', algorithm='l-bfgs',
                    alpha=1e-05, batch_size=200, hidden_layer_sizes=(5, 2),
                    learning_rate='constant', learning_rate_init=0.5,
                    max_iter=200, power_t=0.5, random_state=1, shuffle=False,
                    tol=1e-05, verbose=False, warm_start=False)

After fitting (training), the model can predict labels for new samples::

    >>> clf.predict([[2., 2.], [-1., -2.]]) 
    array([1, 0])

MLP can fit a non-linear model to the training data. ``clf.layers_coef_`` 
contains the weight matrices that constitute the model parameters::

    >>> [coef.shape for coef in clf.layers_coef_]
    [(2, 5), (5, 2), (2, 1)]

To get the raw values before applying the output activation function, run the
following command,

use :meth:`MultilayerPerceptronClassifier.decision_function`::

    >>> clf.decision_function([[2., 2.], [1., 2.]])  # doctest: +ELLIPSIS
    array([ 11.55...,  11.55...])

Currently, :class:`MultilayerPerceptronClassifier` supports only the 
Cross-Entropy loss function, which allows probability estimates by running the 
``predict_proba`` method.

MLP trains using backpropagation. For classification, it minimizes the 
Cross-Entropy loss function, giving a vector of probability estimates 
:math:`P(y|x)` per sample :math:`x`:: 

    >>> clf.predict_proba([[2., 2.], [1., 2.]])  # doctest: +ELLIPSIS
    array([[  9.5...e-06,   9.99...e-01],
           [  9.5...e-06,   9.99...e-01]])

:class:`MultilayerPerceptronClassifier` supports multi-class classification by 
applying `Softmax <http://en.wikipedia.org/wiki/Softmax_activation_function>`_
as the output function. 

Further, the algorithm supports :ref:`multi-label classification <multiclass>` 
in which a sample can belong to more than one class. For each class, the output 
of :meth:`MultilayerPerceptronClassifier.decision_function` passes through the 
logistic function. Values larger or equal to `0.5` are rounded to `1`, 
otherwise to `0`. For a predicted output of a sample, the indices where the
value is `1` represents the assigned classes of that samples::

    >>> X = [[0., 0.], [1., 1.]]
    >>> y = [[0, 1], [1]]
    >>> clf = MultilayerPerceptronClassifier(hidden_layer_sizes=(15,), random_state=1)
    >>> clf.fit(X, y)
    MultilayerPerceptronClassifier(activation='relu', algorithm='l-bfgs',
                    alpha=1e-05, batch_size=200, hidden_layer_sizes=(15,),
                    learning_rate='constant', learning_rate_init=0.5,
                    max_iter=200, power_t=0.5, random_state=1, shuffle=False,
                    tol=1e-05, verbose=False, warm_start=False)
    >>> clf.predict([1., 2.])
    [(1,)]
    >>> clf.predict([0., 0.])
    [(0, 1)]

See the examples below and the doc string of 
:meth:`MultilayerPerceptronClassifier.fit` for further information.

.. topic:: Examples:

 * :ref:`example_plot_mlp_alpha.py`
 * :ref:`example_plot_mlp_nonlinear.py`


Regression
==========

Class :class:`MultilayerPerceptronRegressor` implements  
a multi layer perceptron (MLP) that trains using backpropagation with no 
activation function in the output layer. Therefore, it uses the square error as 
the loss function, and the output is a set of continuous values.

:class:`MultilayerPerceptronRegressor` also supports multi-output regression, in 
which a sample can have more than one target.


Algorithms
==========

MLP trains using either `Stochastic Gradient Descent <http://en.wikipedia.org/wiki/Stochastic_gradient_descent>`_ 
or `L-BFGS <http://en.wikipedia.org/wiki/Limited-memory_BFGS>`_. 
Stochastic Gradient Descent (SGD)  computes the gradient of the loss function 
with respect to a parameter that needs adaptation, i.e.

.. math::

    w \leftarrow w - \eta (\alpha \frac{\partial R(w)}{\partial w}
    + \frac{\partial Loss}{\partial w})
    
where :math:`\eta` is the learning rate which controls the step-size in
the parameter space search.  :math:`Loss` is the loss function used for the network.

With SGD, training supports online and mini-batch learning.

More details can be seen in the documentation of 
`SGD <http://scikit-learn.org/stable/modules/sgd.html>`_ 

L-BFGS is a fast learning algorithm that approximates the Hessian matrix which 
is the second-order partial derivative of a function. Further it approximates 
the inverse of the Hessian matrix to perform parameter update. 
The implementation uses the Scipy version of 
`L-BFGS <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_..

If the selected algorithm is 'L-BFGS', training does not support online nor 
mini-batch learning.


Complexity
==========

Suppose there are :math:`n` training samples, :math:`m` features, :math:`k` 
hidden layers, each containing :math:`h` neurons - for simplicity, and :math:`o`
output neurons.  The time complexity of backpropogation is 
:math:`O(n\cdot m \cdot h^k \cdot o \cdot i)`, where :math:`i` is the number 
of iterations. Since backpropagation has a high time complexity, it is advisable
to start with smaller number of hidden neurons and few hidden layers for 
training.


Mathematical formulation
========================

Given a set of training examples :math:`(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)` 
where :math:`x_i \in \mathbf{R}^n` and :math:`y_i \in \{0, 1\}`, a one hidden 
layer mlp learns the score function :math:`f(x) = W_2^T g(W_1^T x + b_1) + b_2` 
where :math:`W_1, W_2 \in \mathbf{R}^m` and :math:`b_1, b_2 \in \mathbf{R}` are 
model parameters. :math:`W_1, W_2` represent the weights of the input layer and 
hidden layer, resepctively; and :math:`b_1, b_2` represent the bias vectors added 
to the hidden layer and the output layer, respectively. 
:math:`g(\cdot) : R \rightarrow R` is the activation function, set by default as 
the hyperbolic tan. It is given as,

.. math::
      g(x)= \frac{e^x-e^{-x}}{e^x+e^{-x}}

For binary classification, :math:`f(x)` passes through the logistic function
:math:`g(x)= 1/(1+e^{-x})` to get  output values between zero and one. A threshold,
set to 0.5, would assign samples of outputs larger or equal 0.5 to the positive 
class, and the rest to the negative class.

If there are more than two classes, :math:`f(x)` would instead pass through
the softmax function, which is written as,

.. math::
      \text{Softmax} = \frac{\exp(W_i^Tx)}{\sum_{l=1}^k\exp(W_l^Tx)} 

where :math:`W_i` contains the weights incident to the output neuron representing
class :math:`i`, and :math:`K` is the number of classes. The result is a vector 
containing the probabilities that sample :math:`x` belong to each class. The 
output is the class with the highest probability.

In regression, the output remains as :math:`f(x)`; therefore, there is no
output activation function.

MLP uses different loss functions depending on the problem type. The loss 
function for classification is Cross-Entropy, which is given as,

.. math::

    Loss(x,y,W) = -y \ln {f(x)}+(1-y)\ln{(1-f(x))} + \alpha ||W||_2^2

where :math: `\alpha ||W||_2^2` is an L2-norm regularization term (aka penalty)
that penalizes model complexity; and :math:`\alpha > 0` is a non-negative 
hyperparameter that controls the magnitude of the penalty term.

For regression, MLP uses the Square Error loss function; written as,

.. math::

    Loss(x,y,W) = \frac{1}{2}||f(x) - y ||_2^2 + \alpha ||W||_2^2


Starting from initial random weights, multi layer perceptron (MLP) minimizes 
the loss function by repeatedly updating these weights. After computing
the loss function value, a backward pass propagates it from the output layer 
to the initial layer, providing each weight parameter with an update value 
meant to decrease the loss function.

In gradient descent, the weight gradient :math:`\nabla W_{loss}` with respect
to the loss function is computed and added to the corresponding :math:`W`.
More formally, this is expressed as,

.. math::

    W^{i+1} = W^i + \epsilon \nabla W^i_{loss}


where :math:`i` is the iteration step, and :math:`\epsilon` is the learning rate 
with a value ranging from 0 and 1. 

The algorithm stops when it reaches the set number of iterations; or 
when the loss value is below a certain, small number.



.. _mlp_tips:

Tips on Practical Use
=====================

  * Multi-layer Perceptron is sensitive to feature scaling, so it
    is highly recommended to scale your data. For example, scale each
    attribute on the input vector X to [0, 1] or [-1, +1], or standardize
    it to have mean 0 and variance 1. Note that you must apply the *same* 
    scaling to the test vector for meaningful results. 
    You can use :class:`StandardScaler`:: for standarization.

      from sklearn.preprocessing import StandardScaler
      scaler = StandardScaler()
      scaler.fit(X_train)  # Don't cheat - fit only on training data
      X_train = scaler.transform(X_train)
      X_test = scaler.transform(X_test)  # apply same transformation to test data

  * Finding a reasonable regularization term :math:`\alpha` is
    best done using :class:`GridSearchCV`, usually in the
    range ``10.0 ** -np.arange(1, 7)``.

  * Empirically, we observed that `L-BFGS` converges faster and
    with better solutions than `SGD`. Therefore, if mini-batch
    and online learning are unnecessary, it is best advised
    to set :meth:`MultilayerPerceptronClassifier.algorithm` as
    'l-bfgs'.


.. topic:: References:

    * `"Learning representations by back-propagating errors."
      <http://www.iro.umontreal.ca/~pift6266/A06/refs/backprop_old.pdf>`_
      Rumelhart, David E., Geoffrey E. Hinton, and Ronald J. Williams.

    * `"Stochastic Gradient Descent" <http://leon.bottou.org/projects/sgd>`_ L. Bottou - Website, 2010.
    
    * `"Backpropagation" <http://ufldl.stanford.edu/wiki/index.php/Backpropagation_Algorithm>`_
      Andrew Ng, Jiquan Ngiam, Chuan Yu Foo, Yifan Mai, Caroline Suen - Website, 2011.

    * `"Efficient BackProp" <yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf>`_
      Y. LeCun, L. Bottou, G. Orr, K. Müller - In Neural Networks: Tricks
      of the Trade 1998.
         
