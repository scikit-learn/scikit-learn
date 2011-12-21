.. _polynomial_classifier:

=====================
Polynomial Classifier
=====================

.. currentmodule:: sklearn.polynomial_classifier

The Polynomial Classifier solves the problem of assigning a feature vector :math:`c\in\mathbb{R}^d` to one of :math:`k` classes :math:`\boldsymbol\Omega=\{\Omega_1,\ldots,\Omega_k\}`. Instead of using probability densities for representing the feature space, like e.g. Bayes Classifier, the Polynomial Classifier relies on a linear discriminant function :math:`d_\lambda(\mathbf c)` for each of the :math:`k` classes:

.. math::

   d_\lambda(\mathbf{c}) = \mathbf{a}^T_\lambda \boldsymbol\varphi(\mathbf c), \qquad \lambda =1,\ldots, k

where :math:`\mathbf a_\lambda` is a parameter vector and :math:`\boldsymbol\varphi(\mathbf c)\in \mathbb{R}^m` transforms the feature :math:`\mathbf c=(c_1,\ldots,c_d)^T` to a polynomial space with degree :math:`p`,

.. math::

   \boldsymbol\varphi(\mathbf c) = (\varphi_1(\mathbf c), \varphi_2(\mathbf c), \ldots, \varphi_m(\mathbf c))^T

e.g. in the linear case (:math:`p=1`):

.. math::

   \boldsymbol\varphi(\mathbf c) = (1, c_1, c_2, \ldots, c_d)^T

or the quadratic form (:math:`p=2`):

.. math::

   \boldsymbol\varphi(\mathbf c) = (1, c_1, c_2, \ldots, c_d, c_1 c_1, c_2 c_1, \ldots, c_d c_d)^T

With increasing degree :math:`p` the size of :math:`\boldsymbol\varphi(\mathbf c)` (= number of polynomials :math:`m`) rises according to :math:`m=\begin{pmatrix} d + p \\ p \end{pmatrix}=\frac{(d+p)!}{d!\, p!}`. Therefore, values for the degree normally ranges in :math:`p\in\{1,2,3\}` to keep :math:`m` as low as possible. It is also a good idea to reduce the dimensionality of the feature :math:`\mathbf c` vector before training and using the Polynomial Classifier, e.g. by applying the Karhunen-Loeve transform (KLT, PCA).

Classification
--------------

A feature :math:`\mathbf c` is assigned to a class :math:`\Omega_j`, :math:`j\in 1,\ldots,k`, by taking the maximum answer of the discriminant function:

.. math::

   \mathbf c \in \Omega_j\colon j = \arg\max_\lambda d_\lambda(\mathbf c)

Instead of separately evaluating the discriminant function for each class, a more compact (matrix-like) representation can be found:

.. math::

   \mathbf d(\mathbf c) = \mathbf A^T \boldsymbol\varphi(\mathbf c)

where :math:`\mathbf A = (\mathbf a_1,\ldots, \mathbf a_k)`, :math:`\mathbf A\in \mathbb{R}^{m\times k}`, is the parameter matrix.

Learning the Parameter Matrix
-----------------------------

The parameter matrix :math:`\mathbf A` is not known in advance and, therefore, is estimated from a set of given examples :math:`\{(\mathbf c^{(i)}, \boldsymbol\delta^{(i)})\}_{i=1,\ldots,N}`. The "label" :math:`\boldsymbol\delta^{(i)}\in \mathbb{R}^k` is a vector filled by zeros but the j-th entry with :math:`\mathbf c^{(i)}\in \Omega_j` is set to one:

.. math::

   \delta_j^{(i)} = 
   \begin{cases}
   1 & \text{if } \mathbf c^{(i)} \in \Omega_j \\
   0 & \text{otherwise}
   \end{cases}\qquad
   j=1,\ldots,k

The residual is defined as the expectation over all quadratic errors:

.. math::

   \varepsilon(\mathbf A) = E\{ ( \boldsymbol\delta - \mathbf A^T \boldsymbol\varphi(\mathbf c) )^2 \}

with an analytical solution for :math:`\mathbf A`:

.. math::

   \mathbf A^\star = \left( \frac{1}{N} \boldsymbol\Phi \boldsymbol\Phi^T \right)^{-1} \left( \frac{1}{N} \boldsymbol\Phi \boldsymbol\Delta^T \right)

with

.. math::

   \boldsymbol\Phi = \begin{pmatrix} \boldsymbol\varphi(\mathbf c^{(1)}), \ldots, \boldsymbol\varphi(\mathbf c^{(N)}) \end{pmatrix} 

.. math::

   \boldsymbol\Delta = \begin{pmatrix} \boldsymbol\delta^{(1)}, \ldots, \boldsymbol\delta^{(N)} \end{pmatrix}

The direct computation of :math:`\mathbf A^\star` might cause problems because of inverting a singular matrix. Niemann (2003) proposes the application of the Gauss-Jordan algorithm with pivot selection.

.. topic:: Examples:

 * :ref:`examples/polynomial_classifier.py`

.. topic:: References:

 * J. Schürmann (1996).
   Pattern Classification: A Unified View of Statistical and Neural Approaches. 
   John Wiley & Sons, Inc.

 * H. Niemann (2003). 
   Klassifikation von Mustern, University Erlangen-Nürnberg.
   <http://www5.informatik.uni-erlangen.de/fileadmin/Persons/NiemannHeinrich/klassifikation-von-mustern/m00links.html>

