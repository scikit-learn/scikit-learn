# 01-04.py
def GMMBayes(X_test, n_components, covariance_type):
    clf_0 = gmm.GMM(n_components, covariance_type, random_state=0)
    i0 = (y_train == 0)
    clf_0.fit(X_train[i0])

    clf_1 = gmm.GMM(n_components, covariance_type, random_state=0)
    i1 = (y_train == 1)
    clf_1.fit(X_train[i1])

    logL = np.zeros((2, X_test.shape[0]))
    logL[0] = clf_0.score(X_test) + np.log(prior0)
    logL[1] = clf_1.score(X_test) + np.log(prior1)

    y_pred = np.argmax(logL, 0)

    return y_pred
