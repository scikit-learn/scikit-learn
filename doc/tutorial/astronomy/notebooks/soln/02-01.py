# 02-01.py
for i, max_depth in enumerate(max_depth_array):
    # print progress update
    print '%i / %i' % (max_depth, max_depth_array[-1])

    clf = DecisionTreeRegressor(max_depth=max_depth)
    clf.fit(X_train, y_train)

    y_train_pred = clf.predict(X_train)
    y_cv_pred = clf.predict(X_cv)

    train_error[i] = compute_rms_error(y_train_pred, y_train)
    cv_error[i] = compute_rms_error(y_cv_pred, y_cv)
