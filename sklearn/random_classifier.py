
import numpy as np

from .base import BaseEstimator, ClassifierMixin
from .preprocessing import LabelEncoder
from .utils import check_random_state
from .utils.validation import safe_asarray

class RandomClassifier(BaseEstimator, ClassifierMixin):

    def __init__(self, sampling="most_frequent", random_state=0):
        self.sampling = sampling
        self.random_state = random_state

    def fit(self, X, y):
        self.label_encoder_ = LabelEncoder()
        y = self.label_encoder_.fit_transform(y)
        self.class_distrib_ = np.bincount(y) / float(y.shape[0])

    def predict(self, X):
        X = safe_asarray(X)
        n_samples = X.shape[0]
        rs = check_random_state(self.random_state)

        if self.sampling == "most_frequent":
            ret = np.ones(n_samples, dtype=int) * self.class_distrib_.argmax()
        elif self.sampling == "stratified":
            ret = rs.multinomial(1, self.class_distrib_,
                                 size=n_samples).argmax(axis=1)
        elif self.sampling == "uniform":
            ret = rs.randint(len(self.class_distrib_), size=n_samples)
        else:
            raise ValueError("Unknown sampling type.")

        return self.label_encoder_.inverse_transform(ret)
