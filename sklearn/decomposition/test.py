import numpy as np
from .pca import PCA
import pandas as pd

X = np.array([[-1, -1,3,4,-1, -1,3,4], [-2, -1,5,-1, -1,3,4,2], [-3, -2,1,-1, -1,3,4,1],
[1, 1,4,-1, -1,3,4,2], [2, 1,0,-1, -1,3,4,2], [3, 2,10,-1, -1,3,4,10]])

ipca = PCA(n_components = 7, svd_solver= "arpack")

ipca.fit(X)
result = ipca.transform(X)

print result.shape
print ipca.n_components_
