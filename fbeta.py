import numpy as np

from sklearn.metrics import fbeta_score

y_true = [0, 1, 2, 0, 1, 2]
y_pred_empty = [0, 0, 0, 0, 0, 0]
fbeta_score(y_true, y_pred_empty, average="macro", zero_division=np.nan, beta=0.5)
