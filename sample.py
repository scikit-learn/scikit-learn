import numpy as np
from sklearn import linear_model

est_no_intercept = linear_model.Lasso(fit_intercept=False)
est_no_intercept.fit([ [0] ], [1])
#print(est_no_intercept.coef_.shape)
assert est_no_intercept.coef_.shape  == (1,)