import som
import pandas as pd
import vsom
from   sklearn import datasets

# import iris datasets
iris = datasets.load_iris()

# initialize label, if there is no label set labels = None
labels = iris.target

# initialize data
data = pd.DataFrame(iris.data[:, :4])
data.columns = iris.feature_names

# som written entirely in python, som_f written by fortran,
# which is much faster than python
algorithm = "som" 

# Build a map
m = som.build(data,labels,algorithm=algorithm) 

# display the neuron at position (4,4)
som.neuron(m,3,3)			