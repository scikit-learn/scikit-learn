from scikits.learn.svm2 import libsvm
import numpy as np

default_kwargs = {
	'svm_type': 0, # C_SVC
	'kernel_type': 0, # RBF
	'degree': 3,
	'gamma': 0.0,
	'coef0': 0.0,
	'nu': 0.5,
	'C': 1.0,
	'eps': 0.001,
	'p': 0.1,
	'shrinking': 1,
	'probability': 0,
	'nr_weight': 0,
	'weight_label': np.array(0),
	'weight': np.array(0)}

X = np.array( ((0,0), (-1, 1), (1, -1), (1,1), (2,0), (0, 2)) , dtype=np.float)
Y = np.array( (-1, -1, -1, 1, 1, 1 ) )

print libsvm.svm_train_wrap(X, Y, **default_kwargs)
