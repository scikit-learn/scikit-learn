from sklearn.metrics import classification_report
y_true = [0, 0, 2, 0, 0]
y_pred = [0, 2, 2, 0, 0]
target_names = ['class 0', 'class 1', 'class 2']
print(classification_report(y_true, y_pred, target_names=target_names))