import numpy as np
import scipy.sparse as sp_sparse
from .._config import get_config
from .validation import _num_samples


def _one_to_one(feature_names_in):
    return feature_names_in


class _DataAdapter:
    def __init__(self, X, needs_feature_names_in=True):
        self.needs_feature_names_in = needs_feature_names_in
        self.feature_names_in_ = None
        if get_config()['array_out'] == 'ndarray':
            return

        if self.needs_feature_names_in:
            self.feature_names_in_ = self._get_feature_names(X)

    def _get_feature_names(self, X):

        if hasattr(X, "columns"):
            return np.array(X.columns, dtype=object)

        elif hasattr(X, "dims") and isinstance(X.dims, tuple):
            # xarray DataArray
            if len(X.dims) != 2:
                raise ValueError("XArray.DataArray must be 2D")

            return np.array(X.coords[X.dims[1]], dtype=object)

    def check_X(self, X):
        """Get metadata for X that will be transformed"""
        self.n_transform_samples_ = _num_samples(X)
        if get_config()['array_out'] == 'ndarray':
            return self

        # TODO: check for column name consistency

        self.dims_ = None
        self.index_ = None

        if hasattr(X, "columns"):
            # dataframe
            self.index_ = X.index
        elif hasattr(X, "dims") and isinstance(X.dims, tuple):
            # xarray DataArray
            self.dims_ = dims = X.dims
            if len(dims) != 2:
                raise ValueError("XArray.DataArray must be 2D")
            self.index_ = X.coords[dims[0]]
        return self

    def transform(self, X, get_feature_names_out=_one_to_one):
        if not hasattr(self, 'n_transform_samples_'):
            raise ValueError("check_X must be called first")

        array_out = get_config()['array_out']
        if array_out == 'ndarray':
            return X

        if self.n_transform_samples_ != _num_samples(X):
            raise ValueError("The number of samples in fit must be the same "
                             "as transform")

        if self.needs_feature_names_in:
            if self.feature_names_in_ is None:
                feature_names_out = None
            else:
                feature_names_out = get_feature_names_out(
                    self.feature_names_in_)
        else:
            # feature names are not required
            feature_names_out = get_feature_names_out()

        if array_out == 'pandas':
            import pandas as pd
            if sp_sparse.issparse(X):
                return pd.DataFrame.sparse.from_spmatrix(
                    X, index=self.index_, columns=feature_names_out)

            # not sparse
            return pd.DataFrame(X, index=self.index_,
                                columns=feature_names_out)
        elif array_out == 'xarray':
            import xarray as xr

            dims = self.dims_
            index = self.index_
            if dims is None:
                dims = ('index', 'columns')
            if index is None:
                index = np.arange(self.n_transform_samples_)
            # sparse xarray
            if sp_sparse.issparse(X):
                # uses pydata/sparse
                import sparse as pydata_sparse
                X = pydata_sparse.COO.from_scipy_sparse(X)
            return xr.DataArray(X, dims=dims,
                                coords={dims[0]: index,
                                        dims[1]: feature_names_out})
        return X


class _ManyDataAdapter(_DataAdapter):
    def __init__(self, Xs):
        """Xs is a list of arrays or matrics"""
        if get_config()['array_out'] == 'ndarray':
            return
        self.needs_feature_names_in = True

        self.adapters = [
            _DataAdapter(X, needs_feature_names_in=True) for X in Xs]

        feature_name_ins = [adapter.feature_names_in_
                            for adapter in self.adapters]
        if any(names is None for names in feature_name_ins):
            self.feature_names_in_ = None
        else:
            self.feature_names_in_ = np.concatenate(feature_name_ins)

    def check_X(self, Xs):
        for X, adapter in zip(Xs, self.adapters):
            adapter.check_X(X)

        # TODO: Check this
        # assume the metadata is the same for all data
        first_adapter = self.adapters[0]
        self.index_ = first_adapter.index_
        self.dims_ = first_adapter.dims_
        self.n_transform_samples_ = first_adapter.n_transform_samples_

        return self
