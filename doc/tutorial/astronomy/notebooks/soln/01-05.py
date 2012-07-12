# 01-05.py
y_pred_gmm = GMMBayes(X_test, 5, 'full')
y_pred_gnb = gnb.predict(X_test)
