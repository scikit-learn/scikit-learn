.. _neural_network:

==================================
Neural network models (supervised)
==================================

.. currentmodule:: sklearn.neural_network


.. _multilayer_perceptron:
 
Multi-layer Perceptron
======================

**Multi-layer Perceptron (MLP)** is a supervised algorithm that tries to learn 
the relationship between the input features :math:`X` and the target :math:`y`. 
Here we implemented a single-layer feedforward network as shown in Figure 1.

.. figure:: ../images/multilayerperceptron_network.png
   :align: center

   **Figure 1 : Single-layer feedforward network**

Across the module, :math:`W_1` is the weight matrix between the input layer and
the hidden ``coef_hidden_``, :math:`W_2` is the weight matrix between the 
hidden layer and output ``coef_output_`` :math:`b_1` is the intercept vector 
between the input layer and the hidden ``intercept_hidden_``, and :math:`b_2` 
is the intercept vector between the hidden layer and output ``intercept_output_``.

It can be regarded as a sequence of connected Logistic Regression models. 
Recalling Figure 1, The input layer and the hidden layer constitute one 
Logistic Regression model. The values of the hidden layer form the input 
of the next Logistic Regression model composing of that layer and the output 
layer. This algorithm, however, is limited to one single hidden layer.

The advantages of Multi-layer Perceptron are:

    + Capability to learn complex/non-linear models.

    + Flexibility (parameters can be tuned to address different problems).
      
    + Scalability (the readability of the code provides users opportunities to 
      develop more sophisticated networks).
      
The disadvantages of Multi-layer Perceptron include:

    + Hidden layers in MLP make the loss function non-convex (has more than 
      one local minimum), and therefore solutions may not converge to the 
      global minimum.

    + MLP suffers from the Backpropagation diffusion problem; layers far from 
      the output update with decreasing momentum, leading to slow convergence.

    + MLP requires tuning a number of hyperparameters such as the number of 
      hidden neurons and iterations.

    + MLP is sensitive to feature scaling.

Please see the :ref:`Tips on Practical Use <mlp_tips>` section that addresses 
some of the aforementioned disadvantages.


Classification
==============

The class :class:`MultilayerPerceptronClassifier` implements single-layer 
feedforward network that trains using Backpropagation. 

Like all classifiers, MLP trains on two arrays: an array X
of size [n_samples, n_features] holding the training samples, and an
array y of size [n_samples] holding the target values (class labels)
for the training samples::

    >>> from sklearn.neural_network import MultilayerPerceptronClassifier
    >>> X = [[0., 0.], [1., 1.]]
    >>> y = [0, 1]
    >>> clf = MultilayerPerceptronClassifier(n_hidden=2)
    >>> clf.fit(X, y)
    MultilayerPerceptronClassifier(activation='logistic', algorithm='l-bfgs',
                alpha=1e-05, batch_size=200, eta0=0.5,
                learning_rate='constant', max_iter=200, n_hidden=2,
                power_t=0.25, random_state=None, shuffle=False, tol=1e-05,
                verbose=False, warm_start=False)

After fitting (training), the model can predict labels for new samples::

    >>> clf.predict([[2., 2.], [1., 2.]]) 
    >>> array([1, 1])

MLP can fit a non-linear model to the training data. The members 
``clf.coef_hidden_`` and ``clf.coef_output_`` (denoted as :math:`W_1` and 
:math:`W_2`  in Figure 1, respectively) constitute the model parameters::

    >>> clf.coef_hidden_ 
    array([[-3.03439668,  3.10077632],
           [-3.52570416,  3.11540912]])
    >>> clf.coef_output_
    array([[-9.85365392],
           [ 9.94168421]])


The members ``clf.intercept_hidden_`` and ``clf.intercept_output_`` hold the 
intercepts (denoted as :math:`b_1` and :math:`b_2` in Figure 1, respectively):


To get the raw values before the application of an output function (which is 
applied at the last/right-most layer of the classifier)
use :meth:`MultilayerPerceptronClassifier.decision_function`::

    >>> clf.decision_function([[2., 2.], [1., 2.]])
    array([9.80650883, 9.77826778])



Currently, :class:`MultilayerPerceptronClassifier` supports only the 
Cross-Entropy loss function, which allows probability estimates by running the 
``predict_proba`` method.

MLP trains using Backpropagation. For classification, it minimizes the 
Cross-Entropy loss function, giving a vector of probability estimates 
:math:`P(y|x)` per sample :math:`x`:: 

    >>> clf.predict_proba([[2., 2.], [1., 2.]])
    array([[5.50888114e-05, 9.99944911e-01],
           [5.66666648e-05, 9.99943333e-01]])

:class:`MultilayerPerceptronClassifier` supports multi-class classification by 
applying `Softmax <http://en.wikipedia.org/wiki/Softmax_activation_function>`_
as the output function. 

Further, the algorithm supports :ref:`multi-label classification <multiclass>` 
in which a sample can belong to more than one class. For each class, the output 
of :meth:`MultilayerPerceptronClassifier.decision_function` passes through the 
logistic function. Values larger or equal to `0.5` are rounded to `1`, 
otherwise `0`. Classes with value `1` are returned in the prediction::

    >>> X = [[0., 0.], [1., 1.]]
    >>> y = [[0, 1], [1]]
    >>> clf = MultilayerPerceptronClassifier(n_hidden=2)
    >>> clf.fit(X, y)
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

The class :class:`MultilayerPerceptronRegressor` implements a single-layer 
feedforward network that trains using Backpropagation. It only differs from 
:class:`MultilayerPerceptronClassifier`, in that there is no output gate 
function and the loss function is the Square Error. The output is the result 
returned by math:`MultilayerPerceptronRegressor.decision_function`


:class:`MultilayerPerceptronRegressor` supports multi-output regression, in 
which a sample has more than one target.


Algorithms
==========

MLP trains either by `Stochastic Gradient Descent <http://en.wikipedia.org/wiki/Stochastic_gradient_descent>`_ 
or `L-BFGS <http://en.wikipedia.org/wiki/Limited-memory_BFGS>`_. 
Stochastic Gradient Descent (SGD)  computes the gradient of the loss function 
with respect to a parameter that needs adaptation, i.e.

.. math::

    w \leftarrow w - \eta (\alpha \frac{\partial R(w)}{\partial w}
    + \frac{\partial Loss}{\partial w})
    
where :math:`\eta` is the learning rate which controls the step-size in
the parameter space.  :math:`Loss` is the loss function used for the network.

With SGD, training supports online and mini-batch learning.

More details can be seen in the documentation of 
`SGD <http://scikit-learn.org/stable/modules/sgd.html>`_ 

L-BFGS is a powerful method for parameter search. It computes the Hessian 
matrix which is the second-order partial derivative of a function. Further it 
approximates the inverse of the Hessian matrix to perform
parameter update. The implementation uses the Scipy's version of 
`L-BFGS <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html>`_..

If the selected algorithm is 'L-BFGS', training would not support online nor 
mini-batch learning.


Complexity
==========

Single-layer feedforward network has a time complexity that is linear
to the number of training examples. If X is a matrix of size (n, p), and 
the number of hidden neurons is m, training has a cost of 
:math:`O(k n \bar p m)`, where k is the number of iterations (epochs) and 
:math:`\bar p` is the average number of non-zero attributes per sample.

Mathematical formulation
========================

Given a set of training examples :math:`(x_1, y_1), \ldots, (x_n, y_n)` where
:math:`x_i \in \mathbf{R}^n` and :math:`y_i \in \{0, 1\}`, our goal is to
learn a scoring function :math:`\hat{y} = g(W_2^T f(W_1^T x + b_1) + b_2)` 
with model parameters :math:`W_1, W_2 \in \mathbf{R}^m` and intercepts 
:math:`b_1, b_2 \in \mathbf{R}`. `f()` and `g()` are
activation functions. In binary classification, `g()` is the logistic function,

.. math::
      \text{logistic} = 1/(1+e^{-x})

so predictions larger than or equal 0.5 are set to 1, otherwise to 0. The 
output activation for multi-classification is the Softmax,

.. math::
      \text{Softmax} = \frac{\exp(W_j^Tx)}{\sum_{l=1}^k\exp(W_l^Tx)} 

where :math:`W_i` is  the incoming weights to output neuron :math:`i`, and `K`
 is the number of classes. For a sample :math:`x`, ``Softmax`` returns the 
 class with the highest probability.


For classification, it finds the model parameters by minimizing the regularized
Cross-Entropy loss function,

.. math::

    \text{Loss} = -y \ln {\hat{y}}+(1-y)\ln{(1-\hat{y})} + \alpha ||W||_2^2

where :math: `\alpha ||W||_2^2` is an L2-norm regularization term (aka penalty)
that penalizes model complexity; :math: `W` comprises :math: `W_1` and :math: 
`W_2`;  and :math:`\alpha > 0` is a non-negative hyperparameter that controls
the magnitude of the penalty term.

For regression, it finds the model parameters by minimizing the regularized 
Square Error loss function,

.. math::

    \text{Loss} = \frac{1}{2}||W_2 f(W_1 X + b_1) + b_2 - X ||_2^2 + \alpha ||W||_2^2


.. _mlp_tips:

Tips on Practical Use
=====================

  * Multi-layer Perceptron is sensitive to feature scaling, so it
    is highly recommended to scale your data. For example, scale each
    attribute on the input vector X to [0, 1] or [-1, +1], or standardize
    it to have mean 0 and variance 1. Note that the *same* scaling
    must be applied to the test vector to obtain meaningful
    results. This can be easily done using :class:`StandardScaler`::

      from sklearn.preprocessing import StandardScaler
      scaler = StandardScaler()
      scaler.fit(X_train)  # Don't cheat - fit only on training data
      X_train = scaler.transform(X_train)
      X_test = scaler.transform(X_test)  # apply same transformation to test data

    If your attributes have an intrinsic scale (e.g. word frequencies or
    indicator features) scaling is not needed.

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
         
