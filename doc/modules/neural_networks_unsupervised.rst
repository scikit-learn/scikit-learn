.. _neural_network:

==================================
Neural network models (supervised)
==================================

.. currentmodule:: sklearn.neural_network

Randomly weighted neural networks
=================================

Randomly weighted neural networks (RW-NN) is a supervised learning algorithm
that trains a single-hidden layer feedforward network with the help of randomization. 
It computes :math:`\w1 \in R^{d \times h}`, :math:`\w2 \in R^{h \times o}`, and 
:math:`\b \in R^{d}` to solve the following equation:

.. math::

   g(Xw1 + b)w2 = y

where :math:`g(\cdot): R \rightarrow R` is the activation function; :math:`w1 \in R^{d \times k}`
is the weight parameter vector between the input layer of the network and
the hidden layer;  :math:`w2 \in R^{k \times o}` is the weight parameter vector between the hidden 
layer of the network and the output layer;  :math:`b \in R^{d}` is the intercept vector
for the hidden layer. Figure 1 shows an example of such network.

.. figure:: ../auto_examples/neural_networks/images/plot_slfn_001.png
   :target: ../auto_examples/neural_networks/plot_slfn.html
   :align: center
   :scale: 100%

The algorithm takes the following steps:

  *  Generate the matrices :math:`w1 \in R^{d \times k}` and :math:`b \in R^d` with random values using the uniform distribution;
  *  compute :math:`H = g(Xw1 + b)`; and then
  *  solve for :math:`w2` using the ridge implementation as :math:`(H^T H + (1 / C) * I)^{-1} H^T y` - where 
     `C` is the regularization term.

:math:`k` is the number of hidden neurons. Larger :math:`k` allows for higher capacity to learn complex functions. 
It allows the neural network to have randomly combined features in the hidden layer that characterize the training dataset.
This technique is shallow (cannot learn highly complex functions) since the errors resulting from solving :math:`w2` using ridge
are not propagated to the previous layer for better approximation.

For classification, one can use a pipeline comprising the :class:`RandomBasisFunction` and :class:`RidgeClassifier` as
shown in the following example::

    >>> from sklearn.neural_network import RandomBasisFunction
    >>> from sklearn.linear_model import RidgeClassifier
    >>> from sklearn.pipeline import make_pipeline

    >>> X = [[0, 0], [1, 1]]
    >>> y = [0, 1]

    >>> reg = make_pipeline(RandomBasisFunction(random_state=1), RidgeClassifier(alpha=0))
    >>> reg.fit(X, y)
    Pipeline(steps=[('randombasisfunction', RandomBasisFunction(activation='tanh', intercept=True, n_outputs=10,
              random_state=1, weight_scale='auto')), ('ridgeclassifier', RidgeClassifier(alpha=0, class_weight=None, copy_X=True, fit_intercept=True,
            max_iter=None, normalize=False, solver='auto', tol=0.001))])

    >>> reg.predict(X)
    array([0, 1])

For regression, one can use a pipeline comprising the :class:`RandomBasisFunction` and :class:`Ridge` as
shown in the following example::

    >>> from sklearn.neural_network import RandomBasisFunction
    >>> from sklearn.linear_model import Ridge
    >>> from sklearn.pipeline import make_pipeline

    >>> X = [[0, 0], [1, 1]]
    >>> y = [0.5, 0.2]

    >>> reg = make_pipeline(RandomBasisFunction(random_state=1), Ridge(alpha=0))
    >>> reg.fit(X, y)
    Pipeline(steps=[('randombasisfunction', RandomBasisFunction(activation='tanh', intercept=True, n_outputs=10,
              random_state=1, weight_scale='auto')), ('ridge', Ridge(alpha=0, copy_X=True, fit_intercept=True, max_iter=None,
       normalize=False, solver='auto', tol=0.001))])

    >>> reg.predict(X)
    array([ 0.5,  0.2])
     
The references below show examples of how tuning some of the hyper-parameters of the pipeline affect the resulting
decision function::

  * :ref:`example_neural_networks_plot_random_neural_network.py`

  * :ref:`example_neural_networks_plot_random_nn_overfitting.py`

.. topic:: References:

  * Schmidt, Wouter F., Martin A. Kraaijveld, and Robert PW Duin. 
    "Feedforward neural networks with random weights." Pattern Recognition, 
    1992. Vol. II. Conference B: Pattern Recognition Methodology and Systems, 
    Proceedings., 11th IAPR International Conference on. IEEE, 1992.

====================================
Neural network models (unsupervised)
====================================

.. currentmodule:: sklearn.neural_network

.. _random_basis_function:

Random basis function
=====================

Random basis function :math: `f(X): R \rightarrow R` that projects matrix
:math: `X` into another feature space where the number of features is less, equal
or higher than the original feature space. The output matrix :math: `H` is
computed as follows:

.. math::

   H = g(Xw + b)

where :math: `g(\cdot): R \rightarrow R` is the activation function, :math: `w`
is the weight parameter vector, and :math: `b` is the intercept vector.

:math: `w \in R^{d \times k}`, and :math: `b \in R^{d}` are generated based
on the uniform distribution scaled between two values, set by the user.


The example code below illustrates using this function::

    >>> from sklearn.neural_network import RandomBasisFunction
    >>> X = [[0, 0], [1, 1]]
    >>> fe = RandomBasisFunction(random_state=1, n_outputs=2)
    >>> fe.fit(X)
    RandomBasisFunction(activation='tanh', intercept=True, n_outputs=2,
              random_state=1, weight_scale='auto')
    >>> fe.transform(X)
    array([[-0.69896184, -0.76098975],
           [-0.97981807, -0.73662692]])

This function can be useful for training some neural network structures as 
described in the next section.


.. _rbm:

Restricted Boltzmann machines
=============================

Restricted Boltzmann machines (RBM) are unsupervised nonlinear feature learners
based on a probabilistic model. The features extracted by an RBM or a hierarchy
of RBMs often give good results when fed into a linear classifier such as a
linear SVM or a perceptron.

The model makes assumptions regarding the distribution of inputs. At the moment,
scikit-learn only provides :class:`BernoulliRBM`, which assumes the inputs are
either binary values or values between 0 and 1, each encoding the probability
that the specific feature would be turned on.

The RBM tries to maximize the likelihood of the data using a particular
graphical model. The parameter learning algorithm used (:ref:`Stochastic
Maximum Likelihood <sml>`) prevents the representations from straying far
from the input data, which makes them capture interesting regularities, but
makes the model less useful for small datasets, and usually not useful for
density estimation.

The method gained popularity for initializing deep neural networks with the
weights of independent RBMs. This method is known as unsupervised pre-training.

.. figure:: ../auto_examples/neural_networks/images/plot_rbm_logistic_classification_001.png
   :target: ../auto_examples/neural_networks/plot_rbm_logistic_classification.html
   :align: center
   :scale: 100%

.. topic:: Examples:

   * :ref:`example_neural_networks_plot_rbm_logistic_classification.py`


Graphical model and parametrization
-----------------------------------

The graphical model of an RBM is a fully-connected bipartite graph.

.. image:: ../images/rbm_graph.png
   :align: center

The nodes are random variables whose states depend on the state of the other
nodes they are connected to. The model is therefore parameterized by the
weights of the connections, as well as one intercept (bias) term for each
visible and hidden unit, ommited from the image for simplicity.

The energy function measures the quality of a joint assignment:

.. math:: 

   E(\mathbf{v}, \mathbf{h}) = \sum_i \sum_j w_{ij}v_ih_j + \sum_i b_iv_i
     + \sum_j c_jh_j

In the formula above, :math:`\mathbf{b}` and :math:`\mathbf{c}` are the
intercept vectors for the visible and hidden layers, respectively. The
joint probability of the model is defined in terms of the energy:

.. math::

   P(\mathbf{v}, \mathbf{h}) = \frac{e^{-E(\mathbf{v}, \mathbf{h})}}{Z}


The word *restricted* refers to the bipartite structure of the model, which
prohibits direct interaction between hidden units, or between visible units.
This means that the following conditional independencies are assumed:

.. math::

   h_i \bot h_j | \mathbf{v} \\
   v_i \bot v_j | \mathbf{h}

The bipartite structure allows for the use of efficient block Gibbs sampling for
inference.

Bernoulli Restricted Boltzmann machines
---------------------------------------

In the :class:`BernoulliRBM`, all units are binary stochastic units. This
means that the input data should either be binary, or real-valued between 0 and
1 signifying the probability that the visible unit would turn on or off. This
is a good model for character recognition, where the interest is on which
pixels are active and which aren't. For images of natural scenes it no longer
fits because of background, depth and the tendency of neighbouring pixels to
take the same values.

The conditional probability distribution of each unit is given by the
logistic sigmoid activation function of the input it receives:

.. math::

   P(v_i=1|\mathbf{h}) = \sigma(\sum_j w_{ij}h_j + b_i) \\
   P(h_i=1|\mathbf{v}) = \sigma(\sum_i w_{ij}v_i + c_j)

where :math:`\sigma` is the logistic sigmoid function:

.. math::

   \sigma(x) = \frac{1}{1 + e^{-x}}

.. _sml:

Stochastic Maximum Likelihood learning
--------------------------------------

The training algorithm implemented in :class:`BernoulliRBM` is known as
Stochastic Maximum Likelihood (SML) or Persistent Contrastive Divergence
(PCD). Optimizing maximum likelihood directly is infeasible because of
the form of the data likelihood:

.. math::

   \log P(v) = \log \sum_h e^{-E(v, h)} - \log \sum_{x, y} e^{-E(x, y)}

For simplicity the equation above is written for a single training example.
The gradient with respect to the weights is formed of two terms corresponding to
the ones above. They are usually known as the positive gradient and the negative
gradient, because of their respective signs.  In this implementation, the
gradients are estimated over mini-batches of samples.

In maximizing the log-likelihood, the positive gradient makes the model prefer
hidden states that are compatible with the observed training data. Because of
the bipartite structure of RBMs, it can be computed efficiently. The
negative gradient, however, is intractable. Its goal is to lower the energy of
joint states that the model prefers, therefore making it stay true to the data.
It can be approximated by Markov chain Monte Carlo using block Gibbs sampling by
iteratively sampling each of :math:`v` and :math:`h` given the other, until the
chain mixes. Samples generated in this way are sometimes refered as fantasy
particles. This is inefficient and it is difficult to determine whether the
Markov chain mixes.

The Contrastive Divergence method suggests to stop the chain after a small
number of iterations, :math:`k`, usually even 1. This method is fast and has
low variance, but the samples are far from the model distribution.

Persistent Contrastive Divergence addresses this. Instead of starting a new
chain each time the gradient is needed, and performing only one Gibbs sampling
step, in PCD we keep a number of chains (fantasy particles) that are updated
:math:`k` Gibbs steps after each weight update. This allows the particles to
explore the space more thoroughly.

.. topic:: References:

    * `"A fast learning algorithm for deep belief nets"
      <http://www.cs.toronto.edu/~hinton/absps/fastnc.pdf>`_
      G. Hinton, S. Osindero, Y.-W. Teh, 2006

    * `"Training Restricted Boltzmann Machines using Approximations to
      the Likelihood Gradient"
      <http://www.cs.toronto.edu/~tijmen/pcd/pcd.pdf>`_
      T. Tieleman, 2008
