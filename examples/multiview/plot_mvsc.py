"""
===========================
Multiview Spectral Clustering
===========================

An example plot of multiview Spectral Clustering, using multiples views from
same data thanks data given in `UCI Machine Learning Repository
<https://archive.ics.uci.edu/ml/datasets/Multiple+Features>`_.
"""
import numpy as np
import pandas
from multiview.mvsc import MVSC


def readData(filename, data_type=0):
    """
    Given a txt file, returns a numpy matrix with the values, according
    to datatype specified in data_type parameters.

    Parameters
    ----------
    filename: string
        Path or name of the txt file.
    data_type: integer, default 0
        Specifies the matrix datatype. If data_type is 0, data loaded will be
        float type. If data_type is 1, matrix datatype will be float.

    Returns
    -------
    output: ndarray
        Matrix with data loaded.
    """
    if data_type != 0 and data_type != 1:
        raise ValueError('data_type must be either 0 or 1. Found value %d '
                         'instead.' % data_type)
    with open(filename) as txtfile:

        result = []
        myreader = txtfile.readlines()

        for row in myreader:
            if data_type == 0:
                result.append([float(x) for x in row.split()])
            elif data_type == 1:
                result.append([int(x) for x in row.split()])
    if data_type == 0:
        return np.array(result, dtype='float')
    else:
        return np.array(result, dtype='int')


#####################################################

fourier = readData("mfeat-fou.txt", 0)
profcorr = readData("mfeat-fac.txt", 1)
pixels = readData("mfeat-pix.txt", 1)
morpho = readData("mfeat-mor.txt", 0)

markers = ['o', '2', '<', '*', 'h', 'x', 'D', '|', '_', 'v']
mypalette = ['green', 'purple', 'pink', 'blue', 'black',
             'brown', 'yellow', 'orange', 'gray', 'red']

distance = [False] * 4

# Ten clusters
mvsc = MVSC(k=10)
clust = mvsc.fit_transform([fourier, profcorr, pixels, morpho], distance)
clustering = clust[0]
labels_count = [np.unique(
    clustering[200 * i: 200 * (i + 1)], return_counts=True)
    for i in range(10)]
clustered = np.zeros((10, 10), dtype=int)
for i in range(10):
    for index, j in enumerate(labels_count[i][0]):
        clustered[i, j] = labels_count[i][1][index]
classes = np.arange(10)
print("############################################")
print("#               TEN CLUSTERS               #")
print("############################################")
print(pandas.DataFrame(clustered, classes, classes))

# Ten clusters and neighbours 2
mvsc = MVSC(k=10, neighbours=2)
clust = mvsc.fit_transform([fourier, profcorr, pixels, morpho], distance)
clustering = clust[0]
labels_count = [np.unique(
    clustering[200 * i: 200 * (i + 1)], return_counts=True)
    for i in range(10)]
clustered = np.zeros((10, 10), dtype=int)
for i in range(10):
    for index, j in enumerate(labels_count[i][0]):
        clustered[i, j] = labels_count[i][1][index]
classes = np.arange(10)
print("############################################")
print("#       TEN CLUSTERS AND NEIGHBOURS 2      #")
print("############################################")
print(pandas.DataFrame(clustered, classes, classes))
