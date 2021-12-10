"""This module is ARFF parser using LIAC-ARFF library implemented in pure Python."""
import itertools
from collections import OrderedDict
from collections.abc import Generator
from contextlib import closing
from io import BytesIO
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
    to [0, 1, 2])

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


def _cast_frame(frame, columns_info, infer_casting=False):
    """Cast the columns of a dataframe using the ARFF metadata.

    In addition, we provide a way to infer the casting to go beyond the
    limitations of the ARFF metadata.

    Parameters
    ----------
    frame : dataframe
        The dataframe to cast.

    columns_info : dict
        The ARFF metadata for the columns of the dataframe.

    infer_casting : bool
        Whether or not to infer the right casting. Indeed, ARFF does not
        provide the right information for two of the types:

        - `"numeric"` and `"real"` data type can also represent integer values.
          Therefore, when `infer_casting` is `True`, we will try to cast the
          columns to integer.
        - `"nominal"` columns are represented with categorical values that are
          only strings. Therefore, when `infer_casting` is `True`, we will try
          to cast these columns into integers or floats if possible.

    Returns
    -------
    frame : dataframe
        The dataframe with the right casting.
    """
    for name in columns_info:
        column_dtype = columns_info[name]["data_type"]
        if column_dtype.lower() == "integer" and infer_casting:
            frame[name] = frame[name].astype("Int64")
        elif column_dtype.lower() in ("numeric", "real") and infer_casting:
            # "numeric" and "real" could represent integer values.
            # LIAC-ARFF will parse those as float. We can try to cast them to int.
            frame[name] = frame[name].astype("Int64", errors="ignore")
        elif column_dtype.lower() == "nominal":
            # convert to categorical columns
            frame[name] = frame[name].astype("category")
            if infer_casting:
                # LIAC-ARFF will always parse categories as strings.
                # We can try to cast first as integer and then float. If that fails,
                # we can then keep them as strings.
                categories = frame[name].cat.categories
                rename_categories = False
                try:
                    categories = [int(category) for category in categories]
                    rename_categories = True
                except ValueError:
                    try:
                        categories = [float(category) for category in categories]
                        rename_categories = True
                    except ValueError:
                        pass
                if rename_categories:
                    frame[name] = frame[name].cat.rename_categories(categories)
    return frame


def _post_process_frame(frame, feature_names, target_names):
    """Post-process the dataframe to output `X, y, frame, nominal_attributes`."""
    X = frame[feature_names]
    if len(target_names) >= 2:
        y = frame[target_names]
    elif len(target_names) == 1:
        y = frame[target_names[0]]
    else:
        y = None
    nominal_attributes = None
    return X, y, frame, nominal_attributes


def _post_process_array(X, y, columns_info, feature_names, target_names):
    """Post-process the `X` and `y` to output `X, y, frame, nominal_attributes`."""
    columns_to_select = feature_names + target_names
    nominal_attributes = {
        name: columns_info[name]["nominal_value"]
        for name in columns_to_select
        if name in columns_to_select and columns_info[name]["data_type"] == "nominal"
    }
    frame = None

    print(X, X.shape, y, y.shape)

    is_classification = {col_name in nominal_attributes for col_name in target_names}
    if not is_classification:  # No target
        pass
    elif all(is_classification):
        y = np.hstack(
            [
                np.take(
                    np.asarray(nominal_attributes.pop(col_name), dtype="O"),
                    y[:, i : i + 1].astype(int, copy=False),
                )
                for i, col_name in enumerate(target_names)
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


def _liac_arff_parser(
    gzip_file,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
    infer_casting,
    shape=None,
):
    def streamer(gzip_file):
        for line in gzip_file:
            yield line.decode("utf-8")

    stream = streamer(gzip_file)

    return_type = _arff.COO if output_arrays_type == "sparse" else _arff.DENSE_GEN
    arff_container = _arff.load(stream, return_type=return_type)

    columns_to_select = feature_names_to_select + target_names_to_select
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

        frame = _cast_frame(frame, columns_info_openml, infer_casting)
        X, y, frame, nominal_attributes = _post_process_frame(
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

        X, y, frame, nominal_attributes = _post_process_array(
            X,
            y,
            columns_info_openml,
            feature_names_to_select,
            target_names_to_select,
        )

    return X, y, frame, nominal_attributes


def _pandas_arff_parser(
    gzip_file,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
    infer_casting,
):
    import pandas as pd

    # read the file until the data section
    for line in gzip_file:
        line_str = line.decode("utf-8")
        if line_str.lower().startswith("@data"):
            break

    # ARFF represent missing values with "?"
    frame = pd.read_csv(
        BytesIO(gzip_file.read()),
        header=None,
        na_values=["?"],
    )
    frame.columns = [name for name in columns_info_openml]
    frame = _cast_frame(frame, columns_info_openml, infer_casting)
    if output_arrays_type == "pandas":
        X, y, frame, nominal_attributes = _post_process_frame(
            frame, feature_names_to_select, target_names_to_select
        )
    else:
        X = frame[feature_names_to_select].to_numpy()
        y = frame[target_names_to_select].to_numpy()
        X, y, frame, nominal_attributes = _post_process_array(
            X, y, columns_info_openml, feature_names_to_select, target_names_to_select
        )
    return X, y, frame, nominal_attributes


def load_arff_from_gzip_file(
    gzip_file,
    parser,
    output_arrays_type,
    columns_info_openml,
    feature_names_to_select,
    target_names_to_select,
    infer_casting,
    shape=None,
):
    """Load a compressed ARFF file using a given parser.

    Parameters
    ----------
    gzip_file : GzipFile instance
        The file compressed to be read.

    parser : {"liac-arff", "pandas"}
        FIXME

    output_array_type : {"numpy", "sparse", "pandas"}
        The type of the arrays that will be returned. The possibilities ara:

        - `"numpy"`: both `X` and `y` will be NumPy arrays;
        - `"sparse"`: `X` will be sparse matrix and `y` will be a NumPy array;
        - `"pandas"`: `X` will be a pandas DataFrame and `y` will be either a
          pandas Series or DataFrame.

    columns_info : dict
        FIXME

    feature_names_to_select : list of str
        A list of the feature names to be selected.

    target_names_to_select : list of str
        A list of the target names to be selected.

    infer_casting : bool
        FIXME

    Returns
    -------
    frame : dataframe or None
        A dataframe containing both `X` and `y`. `None` if
        `output_array_type != "pandas"`.

    X : {ndarray, sparse matrix, dataframe}
        The data matrix.

    y : {ndarray, dataframe, series}
        The target.

    nominal_attributes : list of str or None
        The names of the features that are categorical. `None` if
        `output_array_type == "pandas"`.
    """
    with closing(gzip_file):

        if parser == "liac-arff":
            return _liac_arff_parser(
                gzip_file,
                output_arrays_type,
                columns_info_openml,
                feature_names_to_select,
                target_names_to_select,
                infer_casting,
                shape,
            )
        elif parser == "pandas":
            return _pandas_arff_parser(
                gzip_file,
                output_arrays_type,
                columns_info_openml,
                feature_names_to_select,
                target_names_to_select,
                infer_casting,
            )
