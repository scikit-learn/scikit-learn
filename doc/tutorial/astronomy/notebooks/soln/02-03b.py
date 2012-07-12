#02-03b.py

#------------------------------------------------------------
# compute and plot the outlier fraction
# as a function of number of samples
n_samples_array = np.linspace(50, Ntrain, 20).astype(int)
train_error_2 = np.zeros(n_samples_array.shape)
cv_error_2 = np.zeros(n_samples_array.shape)

for i, n_samples in enumerate(n_samples_array):
    # print progress update
    print ' %i / %i' % (n_samples, Ntrain)

    clf = DecisionTreeRegressor(max_depth=max_depth)
    clf.fit(X_train[:n_samples], y_train[:n_samples])

    y_train_pred = clf.predict(X_train[:n_samples])
    y_cv_pred = clf.predict(X_cv)

    train_error_2[i] = compute_outlier_fraction(y_train_pred,
                                       y_train[:n_samples])
    cv_error_2[i] = compute_outlier_fraction(y_cv_pred, y_cv)
    
pl.figure()
pl.plot(n_samples_array, cv_error_2, label='cross-val error')
pl.plot(n_samples_array, train_error_2, label='training error')

pl.legend(loc=0)
pl.xlabel('number of samples')
pl.ylabel('error')

pl.title('max_depth = %s' % max_depth)
