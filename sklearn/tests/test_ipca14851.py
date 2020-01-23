
import pandas as pd
import numpy as np

from sklearn.decomposition import IncrementalPCA


def test_ipca14851():
    #  PCA object
    ipca = IncrementalPCA(n_components=16)

    #  dummy data - row counts are from the run in the bug report
    df1 = pd.DataFrame(np.random.rand(56320, 16))
    df2 = pd.DataFrame(np.random.rand(58368, 16))

    #  Fit the batch (original code had a text column)
    X = df1.iloc[:, :].values
    ipca.partial_fit(X)
    assert ipca.n_samples_seen_ == 56320

    X = df2.iloc[:, :].values
    ipca.partial_fit(X)
    assert ipca.n_samples_seen_ == 56320 + 58368


