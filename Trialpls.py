import numpy as np
from sklearn.cross_decomposition import PLSRegression
#from scikit-learn import 

X = np.genfromtxt('train_x.csv',delimiter=';', skip_header=0, skip_footer=0, dtype= float)
Y = np.genfromtxt('train_y.csv',delimiter=';', skip_header=0, skip_footer=0, dtype= float)
