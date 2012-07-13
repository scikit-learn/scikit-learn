# 01-03.py
logL = np.zeros((2, Ncrossval))
logL[0] = clf_0.score(X_crossval) + np.log(prior0)
logL[1] = clf_1.score(X_crossval) + np.log(prior1)
