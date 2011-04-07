from .base import BaseEnsemble

class Committee(BaseEnsemble):

    def fit(self, *args, **kwargs):

        raise NotImplementedError("You must fit each estimator separately before adding them to the committee")

    def predict(self, X):
        """
        Average the predictions of all models in committee
        """ 
        norm = 0.
        prediction = np.zeros(X.shape[0])
        for weight, estimator in self:
            prediction += weight * estimator.predict(X)
            norm += weight
        if norm != 0:
            prediction /= norm
        return prediction
