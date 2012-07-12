#02-03a.py

#------------------------------------------------------------
# first compute and plot the outlier fraction as a function
# of max_depth
max_depth_array = np.arange(1, 21)
train_error = np.zeros(len(max_depth_array))
cv_error = np.zeros(len(max_depth_array))

for i, max_depth in enumerate(max_depth_array):
    # print progress update
    print '%i / %i' % (max_depth, max_depth_array[-1])

    clf = DecisionTreeRegressor(max_depth=max_depth)
    clf.fit(X_train, y_train)

    y_train_pred = clf.predict(X_train)
    y_cv_pred = clf.predict(X_cv)

    train_error[i] = compute_outlier_fraction(y_train_pred, y_train)
    cv_error[i] = compute_outlier_fraction(y_cv_pred, y_cv)

pl.figure()
pl.plot(max_depth_array, cv_error, label='cross-val error')
pl.plot(max_depth_array, train_error, label='training error')

pl.legend(loc=0)
pl.xlabel('max depth')
pl.ylabel('error')

# select the value of max_depth which led to the best results
max_depth = max_depth_array[np.argmin(cv_error)]
print "max_depth = %i" % max_depth
