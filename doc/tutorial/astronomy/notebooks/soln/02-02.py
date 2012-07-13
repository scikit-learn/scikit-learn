#02-02.py
for i, n_samples in enumerate(n_samples_array):
    # print progress update
    print ' %i / %i' % (n_samples, Ntrain)

    clf = DecisionTreeRegressor(max_depth=max_depth)
    clf.fit(X_train[:n_samples], y_train[:n_samples])

    y_train_pred = clf.predict(X_train[:n_samples])
    y_cv_pred = clf.predict(X_cv)

    train_error_2[i] = compute_rms_error(y_train_pred,
                                       y_train[:n_samples])
    cv_error_2[i] = compute_rms_error(y_cv_pred, y_cv)
