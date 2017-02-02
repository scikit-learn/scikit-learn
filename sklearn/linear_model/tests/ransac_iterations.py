from sklearn import linear_model
from sklearn.linear_model.ransac import _dynamic_max_trials
import numpy as np

X =np.linspace(0,100,100)[:,None]
print X.shape
Y = 2*X
Y[80:] = np.random.rand(20,1)
Y[:80] += np.random.rand(80,1)

for i in range(10):
    clt_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression(),max_trials=100)
    clt_ransac.fit(X, Y)
    yhat_ransac = clt_ransac.predict(X)
    inlier_mask = clt_ransac.inlier_mask_
    print '*************************'
    print 'actual     trials:' , clt_ransac.n_trials_
    print 'should be dynamic:', _dynamic_max_trials(np.sum(inlier_mask), 100,
                        2, clt_ransac.stop_probability)
    print np.count_nonzero(inlier_mask)
