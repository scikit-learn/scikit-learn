from .base import BaseEnsemble

class Committee(BaseEnsemble):

    def fit(self, *args, **kwargs):

        raise NotImplementedError("You must fit each estimator separately before adding them to the committee")

    def predict(self, X):
        """
        Average the predictions of all models in committee
        """ 
        prediction = np.zeros(X.shape[0])
        for alpha, estimator in self:
            prediction += estimator.predict(X)
        prediction /= len(self)
        return prediction
