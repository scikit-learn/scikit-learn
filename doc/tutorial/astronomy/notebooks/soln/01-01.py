# 01-01.py
clf_0 = gmm.GMM(1, 'diag')
i0 = (y_train == 0)
clf_0.fit(X_train[i0])

clf_1 = gmm.GMM(1, 'diag')
i1 = (y_train == 1)
clf_1.fit(X_train[i1])
