
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import issparse, csc_matrix

from ..base import TransformerMixin
from ..utils import atleast2d_or_csr, safe_mask
from ..externals import six


class FeatureSelectionMixin(six.with_metaclass(ABCMeta, TransformerMixin)):
    def get_support(self, indices=False):
        """
        Return a mask, or list, of the features/indices selected.
        """
        mask = self._get_support_mask()
        return mask if not indices else np.where(mask)[0]

    @abstractmethod
    def _get_support_mask(self):
        """
        Must return a boolean mask indicating which features are selected.
        """

    def transform(self, X):
        """Reduce X to the selected features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the features selected during the \
            elimination.
        """
        X = atleast2d_or_csr(X)
        mask = self._get_support_mask()
        if len(mask) != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")
        return atleast2d_or_csr(X)[:, safe_mask(X, mask)]

    def inverse_transform(self, X):
        """
        Reverse the transformation operation

        Returns `X` with columns of zeros inserted where features would have
        been removed by `transform`.
        """
        support_ = self.get_support()
        if issparse(X):
            X = X.tocsc()
            # insert additional entries in indptr:
            # e.g. if transform changed indptr from [0 2 6 7] to [0 2 3]
            # col_nonzeros here will be [2 0 1] so indptr becomes [0 2 2 3]
            col_nonzeros = self.inverse_transform(np.diff(X.indptr)).ravel()
            indptr = np.concatenate([[0], np.cumsum(col_nonzeros)])
            Xt = csc_matrix((X.data, X.indices, indptr),
                            shape=(X.shape[0], len(indptr) - 1), dtype=X.dtype)
            return Xt
        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support_.size), dtype=X.dtype)
        Xt[:, support_] = X
        return Xt

