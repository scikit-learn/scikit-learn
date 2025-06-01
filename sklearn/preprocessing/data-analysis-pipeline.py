import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import logging
from pathlib import Path
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from typing import Optional, List, Tuple, Union
import sqlite3

logging.basicConfig(level=logging.INFO)


def read_data(filepath: str) -> pd.DataFrame:
    suffix = Path(filepath).suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(filepath)
    elif suffix == ".xlsx":
        return pd.read_excel(filepath)
    elif suffix == ".sql":
        conn = sqlite3.connect(filepath)
        # For example, read table named 'table_name' or run a query:
        df = pd.read_sql("SELECT * FROM table_name", conn)
        conn.close()
        return df
    
    elif suffix == ".json":
        return pd.read_json(filepath)
    elif suffix == ".parquet":
        return pd.read_parquet(filepath)
    else:
        raise ValueError(f"Unsupported file format: {suffix}")


def data_metrics_and_visualise(filepath: str, s: Optional[str] = None):
    df = read_data(filepath)
    numeric_df = df.select_dtypes(include=[np.number])
    mat = numeric_df.to_numpy().T

    if s == "all_numeric_data":
        for row in mat:
            plt.boxplot(row)
            plt.title("Boxplot")
            plt.show()

        for row in mat:
            plt.hist(row, bins=30)
            plt.title("Histogram")
            plt.show()
    df.info()
    return df.describe(), df.head(), df.corr(), df.median(), df.mode()


def preview(filepath: str):
    df = read_data(filepath)
    return df.head(), df.tail()


def clean_data(
    filepath: str,
    drop_cols: Optional[List[str]] = None,
    drop_null: bool = False,
    fill_null_with_zero: bool = False,
    column_to_clean: Optional[str] = None,
    address_col_to_standardize: Optional[str] = None,
    old_str: Optional[str] = None,
    new_str: Optional[str] = None,
) -> pd.DataFrame:
    df = read_data(filepath)
    df.drop_duplicates(inplace=True)

    if drop_cols:
        df.drop(columns=drop_cols, inplace=True)

    if drop_null:
        df.dropna(inplace=True)

    if fill_null_with_zero:
        df.fillna(0, inplace=True)

    if column_to_clean:
        df[column_to_clean] = df[column_to_clean].apply(
            lambda x: re.sub(r"[^a-zA-Z0-9]", "", str(x)) if pd.notnull(x) else x
        )

    if address_col_to_standardize and old_str is not None and new_str is not None:
        df[address_col_to_standardize] = df[address_col_to_standardize].apply(
            lambda x: str(x).replace(old_str, new_str) if pd.notnull(x) else x
        )

    return df


def transform_data(
    filepath: str,
    filepath_test: Optional[str] = None,
    label_encode_cols: Optional[List[str]] = None,
    one_hot_encode_cols: Optional[List[str]] = None,
    one_hot_encode_cols_sklearn: Optional[List[str]] = None,
    scale: bool = False,
    axis_concat: Optional[int] = None,
    concat_cols: Optional[Tuple[Union[str, int], Union[str, int]]] = None,
    anomaly_col_train: Optional[str] = None,
    anomaly_col_test: Optional[str] = None,
):
    df = read_data(filepath)
    df_test = read_data(filepath_test) if filepath_test else None

    if label_encode_cols:
        le = LabelEncoder()
        for col in label_encode_cols:
            df[col] = le.fit_transform(df[col])
            if df_test is not None:
                df_test[col] = le.transform(df_test[col])

    if one_hot_encode_cols_sklearn:
        oe = OneHotEncoder(sparse_output=False)
        for col in one_hot_encode_cols_sklearn:
            encoded = oe.fit_transform(df[[col]])
            df = df.drop(columns=[col])
            df[oe.get_feature_names_out([col])] = encoded

            if df_test is not None:
                encoded_test = oe.transform(df_test[[col]])
                df_test = df_test.drop(columns=[col])
                df_test[oe.get_feature_names_out([col])] = encoded_test

    if one_hot_encode_cols:
        df = pd.get_dummies(df, columns=one_hot_encode_cols)
        if df_test is not None:
            df_test = pd.get_dummies(df_test, columns=one_hot_encode_cols)
            df_test = df_test.reindex(columns=df.columns, fill_value=0)

    if scale:
        scaler = StandardScaler()
        df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

        if df_test is not None:
            df_test = pd.DataFrame(scaler.transform(df_test), columns=df_test.columns)

    if anomaly_col_train is not None:
        scaler = StandardScaler()
        df[anomaly_col_train] = scaler.fit_transform(df[[anomaly_col_train]])
        df = df[
            (df[anomaly_col_train] <= 3) & (df[anomaly_col_train] >= -3)
        ]

        if anomaly_col_test is not None and df_test is not None:
            df_test[anomaly_col_test] = scaler.transform(df_test[[anomaly_col_test]])
            df_test = df_test[
                (df_test[anomaly_col_test] <= 3) & (df_test[anomaly_col_test] >= -3)
            ]

    if concat_cols and axis_concat is not None:
        df = pd.concat(
            [df[[concat_cols[0]]], df[[concat_cols[1]]]],
            axis=axis_concat,
            ignore_index=True,
        )

    return (df, df_test) if df_test is not None else df
