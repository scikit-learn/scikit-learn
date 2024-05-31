import numpy as np
from sklearn.metrics import confusion_matrix

class TauGeneralized:
    def __init__(self, y_true, y_pred, labels=None, sample_weight=None, normalize=True):
        self.y_true = y_true
        self.y_pred = y_pred
        self.labels = labels
        self.sample_weight = sample_weight
        self.normalize = normalize
        self.cm = confusion_matrix(y_true, y_pred, labels=labels, sample_weight=sample_weight, normalize='all' if normalize else None)
        self.model_point, self.perfect_point, self.random_point = self.__measure()

    def __measure(self):
        cm = self.cm
        n_classes = cm.shape[0]
        perfect_point = np.ones(n_classes)
        random_point = np.full(n_classes, 0.5)

        model_point = cm.diagonal() / cm.sum(axis=1) if self.normalize else cm.diagonal()
        return model_point, perfect_point, random_point

    def __get_dist_from_random(self):
        return np.linalg.norm(self.model_point - self.random_point)

    def __get_dist_from_perfect(self):
        return np.linalg.norm(self.model_point - self.perfect_point)

    def get_tau(self):
        dist_upper_bound = np.sqrt(len(self.model_point))  # Assuming space dimensionality from model_point length
        dist_from_perfect = self.__get_dist_from_perfect()
        return 1 - (dist_from_perfect / dist_upper_bound) if self.normalize else dist_from_perfect
