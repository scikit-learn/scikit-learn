import numpy as np
import scipy.sparse as sp_sparse
from .._config import get_config


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


class _ArrayTransformer:

    def __init__(self, X, needs_feature_names_in=True):
        if (get_config()['array_out'] == 'default'
           or not needs_feature_names_in):
            self.feature_names_in = None
        else:
            self.feature_names_in = _get_feature_names(X)
        self.needs_feature_names_in = needs_feature_names_in
        self.dims = getattr(X, "dims", None)

    def transform(self, X, get_feature_names_out=_one_to_one):
        array_out = get_config()['array_out']
        if array_out == 'default':
            return X

        if self.needs_feature_names_in:
            feature_names_out = get_feature_names_out(self.feature_names_in)
        else:
            feature_names_out = get_feature_names_out()

        # no names are found
        if feature_names_out is None:
            if array_out == 'pydata/sparse' and sp_sparse.issparse(X):
                # hack support for pydata sparse
                import sparse as pydata_sparse
                return pydata_sparse.COO.from_scipy_sparse(X)
            else:
                return X

        if array_out.startswith('pandas'):
            import pandas as pd
            if sp_sparse.issparse(X):
                if array_out == 'pandas/pydata/sparse':
                    # extremely experimental based on
                    # https://github.com/TomAugspurger/pandas/tree/33182-sparse-block  # noqa
                    import sparse as pydata_sparse
                    X = pydata_sparse.COO.from_scipy_sparse(X)
                    return pd.DataFrame(X, columns=feature_names_out)

                # use standard pandas interface for sparse
                return pd.DataFrame.sparse.from_spmatrix(
                    X, columns=feature_names_out)

            # not sparse
            return pd.DataFrame(X, columns=feature_names_out)
        elif array_out == 'xarray':
            import xarray as xr

            dims = self.dims
            if dims is None:
                dims = ('index', 'columns')
            # sparse xarray
            if sp_sparse.issparse(X):
                # uses pydata/sparse
                import sparse as pydata_sparse
                X = pydata_sparse.COO.from_scipy_sparse(X)
            return xr.DataArray(X, dims=dims,
                                coords={dims[1]: feature_names_out})

        return X


class _ManyArrayTransformer(_ArrayTransformer):

    def __init__(self, Xs):
        # TODO: Check Xs has the same structure
        transformers = [_ArrayTransformer(X) for X in Xs]
        feature_name_ins = [trans.feature_names_in for trans in transformers]

        if any(names is None for names in feature_name_ins):
            self.feature_names_in = None
        else:
            self.feature_names_in = np.concatenate(feature_name_ins)

        # Assumes all Xs has the same structure
        first_transformer = transformers[0]
        self.dims = first_transformer.dims
