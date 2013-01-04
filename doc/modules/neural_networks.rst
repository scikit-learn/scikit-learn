.. _neural_networks:

====================================
Neural network models (unsupervised)
====================================

.. currentmodule:: sklearn.neural_networks


.. _rbm:

Restricted Boltzmann Machines
=============================

Restricted Boltzmann Machines (RBM) are unsupervised nonlinear feature learners
based on a probabilistic model. The algorithm used in their training is an
approximation that lowers the quality of the probability density estimation
learned, but focuses on learning useful hidden representations.  Because of
this, the features extracted by an RBM give good results when fed into a linear
classifier such as a linear SVM or perceptron.

The learning algorithm prevents the representations from straying far from the
input data, which makes them capture interesting regularities, but makes the
model less useful for small datasets.

Deep neural networks that are notoriously difficult to train from scratch can
be simplified by initializing each layer's weights with the weights of an RBM.

The makes assumptions regarding the distribution of inputs.  At the moment,
scikit-learn only provides :class:`BernoulliRBM`, which assumes the input are
either binary values or values between 0 and 1, each encoding the probability
that the specific feature would be turned on.

Graphical model and parametrization
-----------------------------------

The graphical model of an RBM is a fully-connected bipartite graph.

.. image:: ../images/rbm_graph.png
   :align: center

The nodes are random variables whose states depend on the state of the other
nodes they are connected to.  The model is therefore parameterized by the
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
This makes inference much easier by assuming the following conditional
independencies:

.. math::

   h_i \bot h_j | \mathbf{v} \\
   v_i \bot v_j | \mathbf{h}


Bernoulli Restricted Boltzmann Machines
---------------------------------------

In the :class:`BernoulliRBM`, all units are binary stochastic units.  This
means that the input data should either be binary, or real-valued between 0 and
1 signifying the probability that the visible unit would turn on or off.  This
is a good model for character recognition, where the interest is on which
pixels are active and which aren't.  For images of natural scenes it no longer
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

