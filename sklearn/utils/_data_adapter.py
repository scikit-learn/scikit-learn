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

    def __init__(self, df_adapter, n_samples, index, dims):
        self.df_adapter = df_adapter
        self.n_samples = n_samples
        self.index = index
        self.dims = dims

    def transform(self, X, get_feature_names_out=_one_to_one):
        array_out = get_config()['array_out']
        if array_out == 'default':
            return X

        if self.n_samples != _num_samples(X):
            raise ValueError("The number of samples in fit must be the same "
                             "as transform")

        if self.df_adapter.needs_feature_names_in:
            feature_names_in = self.df_adapter.feature_names_in_
            if feature_names_in is None:
                # try to get feature names from X
                feature_names_in = _get_feature_names(X)

            if feature_names_in is None:
                feature_names_out = None
            else:
                feature_names_out = get_feature_names_out(
                    self.df_adapter.feature_names_in_)
        else:
            # feature names are not required
            feature_names_out = get_feature_names_out()

        # no names found
        if feature_names_out is None:
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
    def __init__(self, needs_feature_names_in=True):
        self.needs_feature_names_in = needs_feature_names_in

    def fit_get_transformer(self, X):
        return self.fit(X).get_transformer(X)

    def fit(self, X):
        self.feature_names_in_ = None
        if get_config()['array_out'] == 'default':
            return self

        if self.needs_feature_names_in:
            self.feature_names_in_ = _get_feature_names(X)

        return self

    def get_transformer(self, X):
        """Get metadata for X that will be transformed"""
        dims, index = None, None
        n_samples = _num_samples(X)

        if get_config()['array_out'] == 'default':
            return _DataTransformer(self, n_samples, index, dims)

        if hasattr(X, "columns"):
            # dataframe
            index = X.index
        elif hasattr(X, "dims") and isinstance(X.dims, tuple):
            # xarray DataArray
            dims = X.dims
            if len(dims) != 2:
                raise ValueError("XArray.DataArray must be 2D")
            index = X.coords[dims[0]]
        return _DataTransformer(self, n_samples, index, dims)


class _ManyDataAdapter(_DataAdapter):
    def __init__(self):
        super().__init__(needs_feature_names_in=True)

    def fit(self, Xs):
        """Xs is a list of arrays or matrics"""
        self.adapters = [
            _DataAdapter(needs_feature_names_in=True).fit(X) for X in Xs]

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

        return _DataTransformer(self, n_samples, index, dims)
