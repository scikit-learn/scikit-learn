from ._encoders import _BaseEncoder


class ImpactEncoder(_BaseEncoder):
    def __init__(self, categories='auto'):
        self.categories = categories

    def fit(self, X, y):
        return self

    def transform(self, X, y=None):
        return X
