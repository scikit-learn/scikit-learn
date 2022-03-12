import itertools
from collections import OrderedDict
from collections.abc import Generator
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse

from ..externals._arff import ArffSparseDataType, ArffContainerType
from ..utils import (
    _chunk_generator,
    check_pandas_support,
    get_chunk_n_rows,
    is_scalar_nan,
)


def _split_sparse_columns(
    arff_data: ArffSparseDataType, include_columns: List
) -> ArffSparseDataType:
    """
    obtains several columns from sparse arff representation. Additionally, the
    column indices are re-labelled, given the columns that are not included.
    (e.g., when including [1, 2, 3], the columns will be relabelled to
    [0, 1, 2])

    Parameters
    ----------
    arff_data : tuple
        A tuple of three lists of equal size; first list indicating the value,
        second the x coordinate and the third the y coordinate.

    include_columns : list
        A list of columns to include.

    Returns
    -------
    arff_data_new : tuple
        Subset of arff data with only the include columns indicated by the
        include_columns argument.
    """
    arff_data_new: ArffSparseDataType = (list(), list(), list())
    reindexed_columns = {
        column_idx: array_idx for array_idx, column_idx in enumerate(include_columns)
    }
    for val, row_idx, col_idx in zip(arff_data[0], arff_data[1], arff_data[2]):
        if col_idx in include_columns:
            arff_data_new[0].append(val)
            arff_data_new[1].append(row_idx)
            arff_data_new[2].append(reindexed_columns[col_idx])
    return arff_data_new


def _sparse_data_to_array(
    arff_data: ArffSparseDataType, include_columns: List
) -> np.ndarray:
    # turns the sparse data back into an array (can't use toarray() function,
    # as this does only work on numeric data)
    num_obs = max(arff_data[1]) + 1
    y_shape = (num_obs, len(include_columns))
    reindexed_columns = {
        column_idx: array_idx for array_idx, column_idx in enumerate(include_columns)
    }
    # TODO: improve for efficiency
    y = np.empty(y_shape, dtype=np.float64)
    for val, row_idx, col_idx in zip(arff_data[0], arff_data[1], arff_data[2]):
        if col_idx in include_columns:
            y[row_idx, reindexed_columns[col_idx]] = val
    return y


def _feature_to_dtype(feature: Dict[str, str]):
    """Map feature to dtype for pandas DataFrame"""
    if feature["data_type"] == "string":
        return object
    elif feature["data_type"] == "nominal":
        return "category"
    # only numeric, integer, real are left
    elif feature["number_of_missing_values"] != "0" or feature["data_type"] in [
        "numeric",
        "real",
    ]:
        # cast to floats when there are any missing values
        return np.float64
    elif feature["data_type"] == "integer":
        return np.int64
    raise ValueError("Unsupported feature: {}".format(feature))


def _convert_arff_data(
    arff: ArffContainerType,
    col_slice_x: List[int],
    col_slice_y: List[int],
    shape: Optional[Tuple] = None,
) -> Tuple:
    """
    converts the arff object into the appropriate matrix type (np.array or
    scipy.sparse.csr_matrix) based on the 'data part' (i.e., in the
    liac-arff dict, the object from the 'data' key)

    Parameters
    ----------
    arff : dict
        As obtained from liac-arff object.

    col_slice_x : list
        The column indices that are sliced from the original array to return
        as X data

    col_slice_y : list
        The column indices that are sliced from the original array to return
        as y data

    Returns
    -------
    X : np.array or scipy.sparse.csr_matrix
    y : np.array
    """
    arff_data = arff["data"]
    if isinstance(arff_data, Generator):
        if shape is None:
            raise ValueError("shape must be provided when arr['data'] is a Generator")
        if shape[0] == -1:
            count = -1
        else:
            count = shape[0] * shape[1]
        data = np.fromiter(
            itertools.chain.from_iterable(arff_data), dtype="float64", count=count
        )
        data = data.reshape(*shape)
        X = data[:, col_slice_x]
        y = data[:, col_slice_y]
        return X, y
    elif isinstance(arff_data, tuple):
        arff_data_X = _split_sparse_columns(arff_data, col_slice_x)
        num_obs = max(arff_data[1]) + 1
        X_shape = (num_obs, len(col_slice_x))
        X = scipy.sparse.coo_matrix(
            (arff_data_X[0], (arff_data_X[1], arff_data_X[2])),
            shape=X_shape,
            dtype=np.float64,
        )
        X = X.tocsr()
        y = _sparse_data_to_array(arff_data, col_slice_y)
        return X, y
    else:
        # This should never happen
        raise ValueError("Unexpected Data Type obtained from arff.")


def _convert_arff_data_dataframe(
    arff: ArffContainerType, columns: List, features_dict: Dict[str, Any]
) -> Tuple:
    """Convert the ARFF object into a pandas DataFrame.

    Parameters
    ----------
    arff : dict
        As obtained from liac-arff object.

    columns : list
        Columns from dataframe to return.

    features_dict : dict
        Maps feature name to feature info from openml.

    Returns
    -------
    result : tuple
        tuple with the resulting dataframe
    """
    pd = check_pandas_support("fetch_openml with as_frame=True")

    attributes = OrderedDict(arff["attributes"])
    arff_columns = list(attributes)

    if not isinstance(arff["data"], Generator):
        raise ValueError(
            "arff['data'] must be a generator when converting to pd.DataFrame."
        )

    # calculate chunksize
    first_row = next(arff["data"])
    first_df = pd.DataFrame([first_row], columns=arff_columns)

    row_bytes = first_df.memory_usage(deep=True).sum()
    chunksize = get_chunk_n_rows(row_bytes)

    # read arff data with chunks
    columns_to_keep = [col for col in arff_columns if col in columns]
    dfs = []
    dfs.append(first_df[columns_to_keep])
    for data in _chunk_generator(arff["data"], chunksize):
        dfs.append(pd.DataFrame(data, columns=arff_columns)[columns_to_keep])
    df = pd.concat(dfs, ignore_index=True)

    for column in columns_to_keep:
        dtype = _feature_to_dtype(features_dict[column])
        if dtype == "category":
            cats_without_missing = [
                cat
                for cat in attributes[column]
                if cat is not None and not is_scalar_nan(cat)
            ]
            dtype = pd.api.types.CategoricalDtype(cats_without_missing)
        df[column] = df[column].astype(dtype, copy=False)
    return (df,)


def _liac_arff_parser(
    arff_container,
    output_arrays_type,
    features_dict,
    data_columns,
    target_columns,
    col_slice_x=None,
    col_slice_y=None,
    shape=None,
):
    if output_arrays_type == "pandas":
        nominal_attributes = None
        columns = data_columns + target_columns
        (frame,) = _convert_arff_data_dataframe(arff_container, columns, features_dict)
        X = frame[data_columns]
        if len(target_columns) >= 2:
            y = frame[target_columns]
        elif len(target_columns) == 1:
            y = frame[target_columns[0]]
        else:
            y = None
    else:
        frame = None
        X, y = _convert_arff_data(arff_container, col_slice_x, col_slice_y, shape)

        nominal_attributes = {
            k: v
            for k, v in arff_container["attributes"]
            if isinstance(v, list) and k in data_columns + target_columns
        }
        is_classification = {
            col_name in nominal_attributes for col_name in target_columns
        }
        if not is_classification:
            # No target
            pass
        elif all(is_classification):
            y = np.hstack(
                [
                    np.take(
                        np.asarray(nominal_attributes.pop(col_name), dtype="O"),
                        y[:, i : i + 1].astype(int, copy=False),
                    )
                    for i, col_name in enumerate(target_columns)
                ]
            )
        elif any(is_classification):
            raise ValueError(
                "Mix of nominal and non-nominal targets is not currently supported"
            )

        # reshape y back to 1-D array, if there is only 1 target column;
        # back to None if there are not target columns
        if y.shape[1] == 1:
            y = y.reshape((-1,))
        elif y.shape[1] == 0:
            y = None

    return X, y, frame, nominal_attributes
