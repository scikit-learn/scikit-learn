import numpy as np
import scipy.sparse as sp_sparse
from .._config import get_config
from .validation import _num_samples


def _one_to_one(feature_names_in):
    return feature_names_in


def _get_feature_names(X):

    if hasattr(X, "columns"):
        return np.array(X.columns, dtype=object)

    elif hasattr(X, "dims") and isinstance(X.dims, tuple):
        # xarray DataArray
        if len(X.dims) != 2:
            raise ValueError("XArray.DataArray must be 2D")

        return np.array(X.coords[X.dims[1]], dtype=object)


class _DataTransformer:

    def __init__(self, feature_names_in, n_samples, index, dims):
        self.feature_names_in = feature_names_in
        self.n_samples = n_samples
        self.index = index
        self.dims = dims

    def transform(self, X, get_feature_names_out=_one_to_one):
        array_out = get_config()['array_out']
        if array_out == 'default' or not self.n_samples:
            return X

        feature_names_in = self.feature_names_in
        if feature_names_in is None:
            feature_names_in = _get_feature_names(X)

        feature_names_out = get_feature_names_out(feature_names_in)

        # no names are found
        if feature_names_out is None:
            if array_out == 'pydata/sparse' and sp_sparse.issparse(X):
                # hack support for pydata sparse
                import sparse as pydata_sparse
                return pydata_sparse.COO.from_scipy_sparse(X)
            else:
                return X

        if array_out == 'pandas':
            import pandas as pd
            if sp_sparse.issparse(X):
                return pd.DataFrame.sparse.from_spmatrix(
                    X, index=self.index, columns=feature_names_out)

            # not sparse
            return pd.DataFrame(X, index=self.index,
                                columns=feature_names_out)
        elif array_out == 'xarray':
            import xarray as xr

            dims = self.dims
            index = self.index
            if dims is None:
                dims = ('index', 'columns')
            if index is None:
                index = np.arange(self.n_samples)
            # sparse xarray
            if sp_sparse.issparse(X):
                # uses pydata/sparse
                import sparse as pydata_sparse
                X = pydata_sparse.COO.from_scipy_sparse(X)
            return xr.DataArray(X, dims=dims,
                                coords={dims[0]: index,
                                        dims[1]: feature_names_out})
        return X


class _DataAdapter:
    def fit_get_transformer(self, X):
        return self.fit(X).get_transformer(X)

    def fit(self, X):
        # When X is None -> do not save any feature names
        if get_config()['array_out'] == 'default' or X is None:
            self.feature_names_in_ = None
            return self

        self.feature_names_in_ = _get_feature_names(X)
        return self

    def get_transformer(self, X):
        """Get metadata for X that will be transformed"""
        # TODO: It is possible to check for column name consistency here
        dims, index = None, None

        try:
            n_samples = _num_samples(X)
        except TypeError:
            return _DataTransformer(None, 0, index, dims)

        if get_config()['array_out'] == 'default':
            return _DataTransformer(self.feature_names_in_, 0, index, dims)

        if hasattr(X, "columns"):
            # dataframe
            index = X.index
        elif hasattr(X, "dims") and isinstance(X.dims, tuple):
            # xarray DataArray
            dims = X.dims
            if len(dims) != 2:
                raise ValueError("XArray.DataArray must be 2D")
            index = X.coords[dims[0]]
        return _DataTransformer(self.feature_names_in_,
                                n_samples, index, dims)


class _ManyDataAdapter(_DataAdapter):

    def fit(self, Xs):
        """Xs is a list of arrays or matrics"""
        self.adapters = [_DataAdapter().fit(X) for X in Xs]

        feature_name_ins = [adapter.feature_names_in_
                            for adapter in self.adapters]
        if any(names is None for names in feature_name_ins):
            self.feature_names_in_ = None
        else:
            self.feature_names_in_ = np.concatenate(feature_name_ins)
        return self

    def get_transformer(self, Xs):
        transformers = [
            adapter.get_transformer(X)
            for X, adapter in zip(Xs, self.adapters)
        ]

        # TODO: assume the metadata is the same for all data
        first_transformer = transformers[0]
        n_samples = first_transformer.n_samples
        index = first_transformer.index
        dims = first_transformer.dims

        return _DataTransformer(self.feature_names_in_,
                                n_samples, index, dims)
