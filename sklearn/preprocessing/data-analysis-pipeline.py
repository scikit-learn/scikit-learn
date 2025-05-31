import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

def read_data(filepath:str):
    if filepath.lower().endswith(".csv"):
        df = pd.read_csv(filepath)


    elif filepath.lower().endswith(".xlsx"):
        df = pd.read_excel(filepath)

    elif filepath.lower().endswith(".sql"):
        df = pd.read_sql(filepath)

    elif filepath.lower().endswith(".json"):
        df = pd.read_json(filepath)

    elif filepath.lower().endswith(".parquet"):
        df = pd.read_parquet(filepath)

    return df





def data_metricsandvisualise(filepath:str , s=None): ##for visualisation of every numeric data column
    df = read_data(filepath)
    mat = np.array(df)
    mat = mat.T
    if s=='all_numeric_data':
        for i in range(mat.shape[0]):
            plt.boxplot(mat[i])
            plt.show()
        for j in range(mat.shape[0]):
            plt.hist(mat[j])
            plt.show()

    return df.describe() , df.head() , df.info() , df.corr()  , df.median() , df.mode()




def preview(filepath:str):
    df = read_data(filepath)
    return df.head() , df.tail()


def clean_data(filepath:str , s=None,c=None , n=None , r=None ,t=None , h=None , g=None): ##r is the  column with non-alpha-numeric-characters we r cleaning and t is a column with street address
    df = read_data(filepath)
    df.drop_duplicates(inplace=True)
    if s=='drop_column':
        df.drop(columns=c , inplace=True)


    elif n=='drop_null':
        df.dropna(inplace=True)

    elif n=='fill_null_with_0':
        df.fillna(0 , inplace=True)

    elif r=='column_with_non-alphanumeric_characters': ##the goal is to remove the non-alphanumeric characters here
        df[r] = df[r].apply(lambda x: re.sub(r'[^a-zA-Z0-9]' , '' , str(x)) if pd.notnull(x) else x)

    elif t=="column_to_standardise":
        df[t] = df[t].apply(lambda x: str(x).replace(h, g) if pd.notnull(x) else x) ##h is string we are standardising g is the standardised string
    
    return df



def transform_data(filepath: str,filepath_test: str = None,label_encode_cols: list = None,one_hot_encode_cols: list = None ,one_hot_encode_cols_scikit_learn:list=None,scale: bool = False,axis_concat: int = None,concat_cols: tuple = None
,anomaly_col_train : str =None , anomaly_col_test : str = None):
    df = read_data(filepath)
    df_test = read_data(filepath_test) if filepath_test else None


    if label_encode_cols:
        le = LabelEncoder()
        for col in label_encode_cols:
            df[col] = le.fit_transform(df[col])
            if df_test is not None:
                df_test[col] = le.transform(df_test[col])

    if one_hot_encode_cols_scikit_learn:
        oe = OneHotEncoder(sparse=False)
        for col in one_hot_encode_cols_scikit_learn:
            df[col] = oe.fit_transform(df[col])
            if df_test is not None:
                df_test[col] = oe.transform(df_test[col])

    
    if one_hot_encode_cols:
        df = pd.get_dummies(df, columns=one_hot_encode_cols)
        if df_test is not None:
            df_test = pd.get_dummies(df_test, columns=one_hot_encode_cols)
            
            df_test = df_test.reindex(columns=df.columns, fill_value=0)

    
    if scale:
        scaler = StandardScaler()
        df_scaled = scaler.fit_transform(df)
        df = pd.DataFrame(df_scaled, columns=df.columns)

        if df_test is not None:
            df_test_scaled = scaler.transform(df_test)
            df_test = pd.DataFrame(df_test_scaled, columns=df_test.columns)
    if anomaly_col_train is not None:
        scaler = StandardScaler()
        df[anomaly_col_train] = scaler.fit_transform(df[[anomaly_col_train]])
        df = df[(df[anomaly_col_train] <= 3) & (df[anomaly_col_train] >=-3)]
        
        if anomaly_col_test is not None and df_test is not None:
            df_test[anomaly_col_test] = scaler.transform(df[[anomaly_col_test]])
            df_test = df_test[(df_test[anomaly_col_test] <= 3) & (df[anomaly_col_test] >=-3)] ###in desciptive stats this is the set threshold but if needed change it


    
    if concat_cols and axis_concat is not None:
        df = pd.concat([df[concat_cols[0]], df[concat_cols[1]]], axis=axis_concat, ignore_index=True)

    return (df, df_test) if df_test is not None else df


    

    

    


        
    








    


    

    

  


















