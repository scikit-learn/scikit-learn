import numpy as np
import scipy.sparse as sp_sparse
from .._config import get_config
from .validation import _num_samples


def _one_to_one(feature_names_in):
    return feature_names_in


class _DataAdapter:
    def __init__(self, needs_feature_names_in=True):
        self.needs_feature_names_in = needs_feature_names_in

    def fit(self, X):
        self.n_samples_ = _num_samples(X)
        if get_config()['array_out'] == 'ndarray':
            return self

        self.dims_ = None
        self.index_ = None
        self.feature_names_in_ = None

        if hasattr(X, "columns"):
            # dataframe
            self.index_ = X.index

            if self.needs_feature_names_in:
                self.feature_names_in_ = np.array(X.columns, dtype=object)

        elif hasattr(X, "dims") and isinstance(X.dims, tuple):
            # xarray DataArray
            self.dims_ = dims = X.dims
            if len(dims) != 2:
                raise ValueError("XArray.DataArray must be 2D")
            self.index_ = X.coords[dims[0]]

            if self.needs_feature_names_in:
                self.feature_names_in_ = np.array(X.coords[dims[1]],
                                                  dtype=object)
        return self

    def transform(self, X, get_feature_names_out=_one_to_one):
        if self.n_samples_ != _num_samples(X):
            raise ValueError("The number of samples in fit must be the same "
                             "as transform")

        array_out = get_config()['array_out']
        if array_out == 'ndarray':
            return X

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
                index = np.arange(self.n_samples_)
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
    def fit(self, Xs):
        """Xs is a list of arrays or matrics"""
        if get_config()['array_out'] == 'ndarray':
            return self

        adapters = [
            _DataAdapter(needs_feature_names_in=True).fit(X) for X in Xs]
        # TODO: Should check this assumption
        # assume the dimensional metadata is the same for all Xs
        first_adapter = adapters[0]
        self.index_ = first_adapter.index_
        self.dims_ = first_adapter.dims_
        self.n_samples_ = first_adapter.n_samples_

        feature_name_ins = [adapter.feature_names_in_ for adapter in adapters]
        if any(names is None for names in feature_name_ins):
            self.feature_names_in_ = None
        else:
            self.feature_names_in_ = np.concatenate(feature_name_ins)

        return self
