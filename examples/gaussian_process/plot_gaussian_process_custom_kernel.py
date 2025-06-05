"""
==========================================================
Custom kernel for Gaussian Process Regression (GPR)
==========================================================
"""

# A new custom kernel can also be created once it subclasses
# :class:`sklearn.gaussian_process.kernels.Kernel` and once the abstract
# methods :meth:`__call__`, :meth:`diag`, and :meth:`is_stationary`
# are implemented.
#

# Here is a step-by-step construction of the Sigmoid Kernel
# .. math::
# k(x_i, x_j)= \tanh(\alpha x_i \cdot x_j + \sigma)
# First, :math:`\alpha` and :math:`\sigma` are read in as hyperparameters::

from sklearn.gaussian_process.kernels import Kernel, Hyperparameter

class MySigmoidKernel(Kernel):
      def __init__(self, alpha=0.05,
                  alpha_bounds=(1e-5, 1e5),
                  sigma=5,
                  sigma_bounds=(1e-5, 1e5)):
         self.alpha = alpha
         self.alpha_bounds = alpha_bounds
         self.sigma = sigma
         self.sigma_bounds = sigma_bounds

      @property
      def hyperparameter_scaling(self):
         return Hyperparameter("alpha", "numeric", self.alpha_bounds)

      @property
      def hyperparameter_shifting(self):
         return Hyperparameter("sigma", "numeric", self.sigma_bounds)



# Note that the custom kernel must be positive-definite for the choice of the
# hyperparameters.

# Now, from the formula, this sigmoid kernel is non-stationary.
# If a custom kernel is instead stationary, it can subclass
# :func:`~sklearn.gaussian_process.kernels.StationaryKernelMixin` and skip this step;
# otherwise, the method :meth:`is_stationary` has to be implemented as follows::

      def is_stationary(self):
         return False

# Next, the diagonal must be calculated. Note that if a custom
# kernel is normalized, i.e. :math:`k(x,x)=1`, then it can subclass
# :func:`~sklearn.gaussian_process.kernels.NormalizedKernelMixin`, in which case
# this step can be skipped::

      def diag(self, X):
         return np.tanh(self.alpha * np.einsum('ij,ij->i', X, X) + self.sigma)

# Finally, for the :meth:`__call__` method, the gradient needs to be
# calculated with respect to the log-transformed hyperparameters.
# In this case, there are two hyperparameters:


# :math:`\alpha` and :math:`\sigma` with respect to which the partial derivatives
# need to be taken (note that if either of them is fixed, then an empty matrix
# is returned as the partial derivative). After the partial derivatives are
# taken, the method :func:`numpy.dstack` combines them into a variable
# `K_gradient`::



      def __call__(self, X, Y=None, eval_gradient=False):
         X = np.atleast_2d(X)
         if Y is None:
               K = np.tanh(self.sigma + self.alpha * np.inner(X, X))
         else:
               if eval_gradient:
                  raise ValueError(
                     "Gradient can only be evaluated when Y is None.")
               K = np.tanh(self.sigma + self.alpha * np.inner(X, Y))

         if eval_gradient:
               # gradient with respect to sigma

            if not self.hyperparameter_shifting.fixed:
                  sigma_gradient = np.empty((K.shape[0], K.shape[1], 1))
                  sigma_gradient[...,0] = self.sigma * (1 - K**2)
            else:
                sigma_gradient = np.empty((X.shape[0], X.shape[1], 0))

            # gradient with respect to log-transformed alpha
            if not self.hyperparameter_scaling.fixed:
                alpha_gradient = np.empty((K.shape[0], K.shape[1], 1))
                alpha_gradient[...,0] = self.alpha * np.inner(X,X) * (1-K**2)
            else:
                alpha_gradient = np.empty((K.shape[0], K.shape[1], 0))

            K_gradient = np.dstack((alpha_gradient, sigma_gradient))
            return K, K_gradient

         else:
               return K

# Once these methods have been implemented, the custom kernel is ready to be used.

# References
# ----------
