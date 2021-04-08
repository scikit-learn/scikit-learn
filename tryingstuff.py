#%%
from sklearn.cluster import SpectralBiclustering
import numpy as np
X = np.array([[1, 1], [2, 1], [1, 0],[4, 7], [3, 5], [3, 6]])
clustering = SpectralBiclustering(n_clusters=2, random_state=0).fit(X)
clustering.row_labels_

clustering.column_labels_

clustering



# %%
