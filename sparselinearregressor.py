#!/usr/bin/env python
# coding: utf-8

# # Implementing the SparseLinearRegressor
# 
# * **Issue number**:30139 
# * **Title of the Issue**: The `input_tags.sparse` flag is often incorrect

# In[46]:


import numpy as np
from scipy import sparse as sp
from sklearn.base import RegressorMixin, MultiOutputMixin
from sklearn.linear_model._base import LinearModel
from sklearn.utils.validation import check_X_y
from scipy.sparse.linalg import lsqr
from scipy import optimize
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed

class SparseLinearRegressor(MultiOutputMixin, RegressorMixin, LinearModel):
    """
    Custom implementation of a linear regression model that supports
    fitting on sparse data and can enforce non-negativity constraints.

    Parameters
    ----------
    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set to False,
        no intercept will be used in calculations (i.e., data is expected to be centered).
        
    copy_X : bool, default=True
        If True, a copy of X will be made; else, it may overwrite X.
        
    n_jobs : int or None, default=None
        The number of jobs to use for the computation. None means 1 unless
        in a joblib.parallel_backend context.
        
    positive : bool, default=False
        If True, the coefficients are constrained to be positive.
    """
    
    def __init__(self, *, fit_intercept=True, copy_X=True, n_jobs=None, positive=False):
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.n_jobs = n_jobs
        self.positive = positive

    def _get_tags(self):
        """
        Retrieve the tags associated with the estimator.

        Returns
        -------
        dict
            Tags related to the model's behavior.
        """
        tags = super()._get_tags()
        tags['sparse'] = not self.positive  # Adjust tag based on positivity constraint
        return tags

    def fit(self, X, y, sample_weight=None):
        """
        Fit the model to the provided data.

        Parameters
        ----------
        X : array-like or sparse matrix of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Individual weights for each sample.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        n_jobs_ = self.n_jobs
        accept_sparse = False if self.positive else ["csr", "csc", "coo"]

        # Validating input data
        X, y = check_X_y(X, y, accept_sparse=accept_sparse, multi_output=True, y_numeric=True)

        # Preprocess data (center and scale)
        X, y, X_offset, y_offset, X_scale = _preprocess_data(X, y, fit_intercept=self.fit_intercept, copy=self.copy_X)

        if self.positive:
            # Use non-negative least squares for fitting
            if y.ndim < 2:
                self.coef_ = optimize.nnls(X, y)[0]
            else:
                outs = Parallel(n_jobs=n_jobs_)(delayed(optimize.nnls)(X, y[:, j]) for j in range(y.shape[1]))
                self.coef_ = np.vstack([out[0] for out in outs])
        elif sp.issparse(X):
            # For sparse inputs, solve using lsqr
            X_offset_scale = X_offset / X_scale
            X_centered = sp.linalg.LinearOperator(
                shape=X.shape,
                matvec=lambda b: X.dot(b) - b.dot(X_offset_scale),
                rmatvec=lambda b: X.T.dot(b) - X_offset_scale * b.sum()
            )
            if y.ndim < 2:
                self.coef_ = lsqr(X_centered, y)[0]
            else:
                outs = Parallel(n_jobs=n_jobs_)(delayed(lsqr)(X_centered, y[:, j].ravel()) for j in range(y.shape[1]))
                self.coef_ = np.vstack([out[0] for out in outs])
        else:
            # For dense inputs, use least squares
            self.coef_, _, self.rank_, self.singular_ = np.linalg.lstsq(X, y, rcond=None)
            self.coef_ = self.coef_.T

        # Ensuring coefficients are 1D if y is 1D
        if y.ndim == 1:
            self.coef_ = np.ravel(self.coef_)

        # Set intercept based on preprocessing
        self._set_intercept(X_offset, y_offset, X_scale)
        return self

# Utility functions needed for fitting
def _preprocess_data(X, y, fit_intercept=True, copy=True, sample_weight=None):
    """
    Preprocess the data by centering it and scaling if necessary.

    Parameters
    ----------
    X : array-like or sparse matrix of shape (n_samples, n_features)
        Training data.
    
    y : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Target values.
    
    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model.

    copy : bool, default=True
        If True, a copy of X will be made; else, it may overwrite X.

    sample_weight : array-like of shape (n_samples,), default=None
        Individual weights for each sample.

    Returns
    -------
    tuple
        Processed X, y, and offset/scale values.
    """
    if fit_intercept:
        # Centering data
        if sp.issparse(X):
            X_offset = np.average(X.toarray(), axis=0, weights=sample_weight)
        else:
            X_offset = np.average(X, axis=0, weights=sample_weight)
        X = X - X_offset

        y_offset = np.average(y, axis=0, weights=sample_weight)
        y = y - y_offset
    else:
        X_offset = np.zeros(X.shape[1], dtype=X.dtype)
        y_offset = 0. if y.ndim == 1 else np.zeros(y.shape[1], dtype=y.dtype)

    # Optionally scale data
    if fit_intercept:
        X_scale = np.ones(X.shape[1], dtype=X.dtype)
    else:
        X_scale = StandardScaler(with_mean=False).fit(X).scale_

    return X, y, X_offset, y_offset, X_scale

# Test Code
if __name__ == "__main__":
    # Sample sparse input data
    X_sparse = sp.csr_matrix([[1, 2], [3, 4], [5, 6]])
    y = np.array([1, 2, 3])

    # Testing with positive=False (should accept sparse)
    reg = SparseLinearRegressor(positive=False)
    tags = reg._get_tags()
    print("Sparse support with positive=False:", tags['sparse'])

    # Checking actual behavior with sparse input (positive=False)
    try:
        reg.fit(X_sparse, y)
        print("Fitting with sparse data (positive=False) succeeded.")
    except ValueError as e:
        print("Error during fit with sparse data (positive=False):", e)

    # Checking actual behavior with sparse input (positive=True)
    try:
        reg = SparseLinearRegressor(positive=True)
        reg.fit(X_sparse.toarray(), y)  # Convert sparse to dense
        print("Fitting with sparse data (positive=True) succeeded after conversion to dense.")
    except ValueError as e:
        print("Error during fit with sparse data (positive=True):", e)


# # Testing with some datas

# In[47]:


X_new = sp.csr_matrix([[2, 3], [4, 5]])

# For positive=True, convert to dense
reg_positive = SparseLinearRegressor(positive=True)
reg_positive.fit(X_sparse.toarray(), y)
predictions_positive = reg_positive.predict(X_new.toarray())  # Convert to dense
print("Predictions with positive=True:", predictions_positive)

# For positive=False, use sparse directly
reg_negative = SparseLinearRegressor(positive=False)
reg_negative.fit(X_sparse, y)
predictions_negative = reg_negative.predict(X_new)
print("Predictions with positive=False:", predictions_negative)


# In[48]:


import numpy as np
from scipy import sparse as sp

# Creating a larger sparse dataset
# For demonstration, let's create a sparse matrix of shape (1000, 50)
# where only 5% of the entries are non-zero.
n_samples = 1000
n_features = 50
density = 0.05  # 5% non-zero entries

# Generating a random sparse matrix
X_large_sparse = sp.rand(n_samples, n_features, density=density, format='csr')
y_large = np.random.rand(n_samples)  # Random target values

# New input data for predictions (also sparse)
X_new = sp.csr_matrix([[2, 0, 3] + [0] * (n_features - 3), 
                        [4, 0, 0] + [0] * (n_features - 3)])

# For positive=True, convert to dense
reg_positive = SparseLinearRegressor(positive=True)

# Fit the model with sparse input converted to dense
reg_positive.fit(X_large_sparse.toarray(), y_large)
predictions_positive = reg_positive.predict(X_new.toarray())  # Convert to dense for prediction
print("Predictions with positive=True:", predictions_positive)

# For positive=False, fit the model with sparse data directly
reg_negative = SparseLinearRegressor(positive=False)
reg_negative.fit(X_large_sparse, y_large)

# Using the sparse new input directly for prediction
predictions_negative = reg_negative.predict(X_new)
print("Predictions with positive=False:", predictions_negative)

