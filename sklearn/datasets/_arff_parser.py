"""Implementation of ARFF parsers: via LIAC-ARFF and pandas."""
import itertools
from collections import OrderedDict
from collections.abc import Generator
from typing import List

import numpy as np
import scipy as sp

from ..externals import _arff
from ..externals._arff import ArffSparseDataType
from ..utils import (
    _chunk_generator,
    check_pandas_support,
    get_chunk_n_rows,
)


def _split_sparse_columns(
    arff_data: ArffSparseDataType, include_columns: List
) -> ArffSparseDataType:
    """Obtains several columns from sparse ARFF representation. Additionally,
    the column indices are re-labelled, given the columns that are not
    included. (e.g., when including [1, 2, 3], the columns will be relabelled
    to [0, 1, 2]).

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


def _cast_frame(frame, columns_info):
    """Cast the columns of a dataframe using the ARFF metadata.

    Parameters
    ----------
    frame : dataframe
        The dataframe to cast.

    columns_info : dict
        The ARFF metadata for the columns of the dataframe.

    Returns
    -------
    frame : dataframe
        The dataframe with the right casting.
    """
    dtypes = {}
    for name in frame.columns:
        column_dtype = columns_info[name]["data_type"]
        if column_dtype.lower() == "integer":
            dtypes[name] = "Int64"
        elif column_dtype.lower() == "nominal":
            dtypes[name] = "category"
        else:
            dtypes[name] = frame.dtypes[name]
    return frame.astype(dtypes)


def _post_process_frame(frame, feature_names, target_names):
    """Post process a dataframe to select the desired columns in `X` and `y`.

    Parameters
    ----------
    frame : dataframe
        The dataframe to split into `X` and `y`.

    feature_names : list of str
        The list of feature names to populate `X`.

    target_names : list of str
        The list of target names to populate `y`.

    Returns
    -------
    X : dataframe
        The dataframe containing the features.

    y : {series, dataframe} or None
        The series or dataframe containing the target.
    """
    X = frame[feature_names]
    if len(target_names) >= 2:
        y = frame[target_names]
    elif len(target_names) == 1:
        y = frame[target_names[0]]
    else:
        y = None
    return X, y


def _liac_arff_parser(
    gzip_file,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
    shape=None,
):
    """ARFF parser using the LIAC-ARFF library coded purely in Python.

    This parser is quite slow but consume a generator.

    Parameters
    ----------
    gzip_file : GzipFile instance
        The file compressed to be read.

    output_array_type : {"numpy", "sparse", "pandas"}
        The type of the arrays that will be returned. The possibilities ara:

        - `"numpy"`: both `X` and `y` will be NumPy arrays;
        - `"sparse"`: `X` will be sparse matrix and `y` will be a NumPy array;
        - `"pandas"`: `X` will be a pandas DataFrame and `y` will be either a
          pandas Series or DataFrame.

    columns_info : dict
        The information provided by OpenML regarding the columns of the ARFF
        file.

    feature_names_to_select : list of str
        A list of the feature names to be selected.

    target_names_to_select : list of str
        A list of the target names to be selected.

    Returns
    -------
    X : {ndarray, sparse matrix, dataframe}
        The data matrix.

    y : {ndarray, dataframe, series}
        The target.

    frame : dataframe or None
        A dataframe containing both `X` and `y`. `None` if
        `output_array_type != "pandas"`.

    nominal_attributes : list of str or None
        The names of the features that are categorical. `None` if
        `output_array_type == "pandas"`.
    """

    def _io_to_generator(gzip_file):
        for line in gzip_file:
            yield line.decode("utf-8")

    stream = _io_to_generator(gzip_file)

    # find which type (dense or sparse) ARFF type we will have to deal with
    return_type = _arff.COO if output_arrays_type == "sparse" else _arff.DENSE_GEN
    # we should not let LIAC-ARFF to encode the nominal attributes with NumPy
    # arrays to have only numerical values.
    encode_nominal = not (output_arrays_type == "pandas")
    arff_container = _arff.load(
        stream, return_type=return_type, encode_nominal=encode_nominal
    )
    columns_to_select = feature_names_to_select + target_names_to_select

    nominal_attributes = {
        name: categories
        for name, categories in arff_container["attributes"]
        if isinstance(categories, list) and name in columns_to_select
    }
    if output_arrays_type == "pandas":
        pd = check_pandas_support("fetch_openml with as_frame=True")

        columns_info = OrderedDict(arff_container["attributes"])
        columns_names = list(columns_info.keys())

        # calculate chunksize
        first_row = next(arff_container["data"])
        first_df = pd.DataFrame([first_row], columns=columns_names)

        row_bytes = first_df.memory_usage(deep=True).sum()
        chunksize = get_chunk_n_rows(row_bytes)

        # read arff data with chunks
        columns_to_keep = [col for col in columns_names if col in columns_to_select]
        dfs = [first_df[columns_to_keep]]
        for data in _chunk_generator(arff_container["data"], chunksize):
            dfs.append(pd.DataFrame(data, columns=columns_names)[columns_to_keep])
        frame = pd.concat(dfs, ignore_index=True)
        del dfs, first_df

        frame = _cast_frame(frame, columns_info_openml)
        X, y = _post_process_frame(
            frame, feature_names_to_select, target_names_to_select
        )
    else:
        arff_data = arff_container["data"]

        feature_indices_to_select = [
            int(columns_info_openml[col_name]["index"])
            for col_name in feature_names_to_select
        ]
        target_indices_to_select = [
            int(columns_info_openml[col_name]["index"])
            for col_name in target_names_to_select
        ]

        if isinstance(arff_data, Generator):
            if shape is None:
                raise ValueError(
                    "shape must be provided when arr['data'] is a Generator"
                )
            if shape[0] == -1:
                count = -1
            else:
                count = shape[0] * shape[1]
            data = np.fromiter(
                itertools.chain.from_iterable(arff_data),
                dtype="float64",
                count=count,
            )
            data = data.reshape(*shape)
            X = data[:, feature_indices_to_select]
            y = data[:, target_indices_to_select]
        elif isinstance(arff_data, tuple):
            arff_data_X = _split_sparse_columns(arff_data, feature_indices_to_select)
            num_obs = max(arff_data[1]) + 1
            X_shape = (num_obs, len(feature_indices_to_select))
            X = sp.sparse.coo_matrix(
                (arff_data_X[0], (arff_data_X[1], arff_data_X[2])),
                shape=X_shape,
                dtype=np.float64,
            )
            X = X.tocsr()
            y = _sparse_data_to_array(arff_data, target_indices_to_select)
        else:
            # This should never happen
            raise ValueError("Unexpected Data Type obtained from arff.")

        is_classification = {
            col_name in nominal_attributes for col_name in target_names_to_select
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
                    for i, col_name in enumerate(target_names_to_select)
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

    if output_arrays_type == "pandas":
        return X, y, frame, None
    return X, y, None, nominal_attributes


def _pandas_arff_parser(
    gzip_file,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
):
    """ARFF parser using `pandas.read_csv`.

    The metadata being fetch directly to OpenML, one can skip the metadata from
    the ARFF file and read the data as a CSV file.

    Parameters
    ----------
    gzip_file : GzipFile instance
        The file compressed to be read.

    output_array_type : {"numpy", "sparse", "pandas"}
        The type of the arrays that will be returned. The possibilities ara:

        - `"numpy"`: both `X` and `y` will be NumPy arrays;
        - `"sparse"`: `X` will be sparse matrix and `y` will be a NumPy array;
        - `"pandas"`: `X` will be a pandas DataFrame and `y` will be either a
          pandas Series or DataFrame.

    columns_info : dict
        The information provided by OpenML regarding the columns of the ARFF
        file.

    feature_names_to_select : list of str
        A list of the feature names to be selected.

    target_names_to_select : list of str
        A list of the target names to be selected.

    Returns
    -------
    X : {ndarray, sparse matrix, dataframe}
        The data matrix.

    y : {ndarray, dataframe, series}
        The target.

    frame : dataframe or None
        A dataframe containing both `X` and `y`. `None` if
        `output_array_type != "pandas"`.

    nominal_attributes : list of str or None
        The names of the features that are categorical. `None` if
        `output_array_type == "pandas"`.
    """
    import pandas as pd

    # read the file until the data section
    for line in gzip_file:
        if line.decode("utf-8").lower().startswith("@data"):
            break

    # ARFF represent missing values with "?"
    frame = pd.read_csv(
        gzip_file,
        header=None,
        na_values=["?"],  # missing values are represented by `?`
        comment="%",  # skip line starting by `%` since they are comments
    )
    frame.columns = [name for name in columns_info_openml]

    columns_to_select = feature_names_to_select + target_names_to_select
    columns_to_keep = [col for col in frame.columns if col in columns_to_select]
    frame = frame[columns_to_keep]

    frame = _cast_frame(frame, columns_info_openml)
    X, y = _post_process_frame(frame, feature_names_to_select, target_names_to_select)
    nominal_attributes = {
        col_name: frame[col_name].cat.categories.tolist()
        for col_name in frame.columns
        if hasattr(frame[col_name].dtype, "is_dtype")
        and frame[col_name].dtype.is_dtype("category")
    }

    if output_arrays_type == "pandas":
        return X, y, frame, None
    else:
        # FIXME: we should only use `.to_numpy` when supporting more recent version
        # of pandas.
        X = X.to_numpy() if hasattr(X, "to_numpy") else X.values
        if y is not None:
            y = y.to_numpy() if hasattr(y, "to_numpy") else np.asarray(y)
    return X, y, None, nominal_attributes


def load_arff_from_gzip_file(
    gzip_file,
    parser,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
    shape=None,
):
    """Load a compressed ARFF file using a given parser.

    Parameters
    ----------
    gzip_file : GzipFile instance
        The file compressed to be read.

    parser : {"liac-arff", "pandas"}
        The parser used to parse the ARFF file.

    output_array_type : {"numpy", "sparse", "pandas"}
        The type of the arrays that will be returned. The possibilities ara:

        - `"numpy"`: both `X` and `y` will be NumPy arrays;
        - `"sparse"`: `X` will be sparse matrix and `y` will be a NumPy array;
        - `"pandas"`: `X` will be a pandas DataFrame and `y` will be either a
          pandas Series or DataFrame.

    columns_info : dict
        The information provided by OpenML regarding the columns of the ARFF
        file.

    feature_names_to_select : list of str
        A list of the feature names to be selected.

    target_names_to_select : list of str
        A list of the target names to be selected.

    Returns
    -------
    X : {ndarray, sparse matrix, dataframe}
        The data matrix.

    y : {ndarray, dataframe, series}
        The target.

    frame : dataframe or None
        A dataframe containing both `X` and `y`. `None` if
        `output_array_type != "pandas"`.

    nominal_attributes : list of str or None
        The names of the features that are categorical. `None` if
        `output_array_type == "pandas"`.
    """
    if parser == "liac-arff":
        return _liac_arff_parser(
            gzip_file,
            output_arrays_type,
            columns_info_openml,
            feature_names_to_select,
            target_names_to_select,
            shape,
        )
    elif parser == "pandas":
        return _pandas_arff_parser(
            gzip_file,
            output_arrays_type,
            columns_info_openml,
            feature_names_to_select,
            target_names_to_select,
        )
    else:
        raise ValueError(
            f"Unknown parser: '{parser}'. Should be 'liac-arff' or 'pandas'."
        )
