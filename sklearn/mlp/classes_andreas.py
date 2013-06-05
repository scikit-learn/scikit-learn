from ..base import ClassifierMixin
from .base_andreas import BaseMLP
from ..preprocessing import LabelBinarizer


class MLPClassifierA(BaseMLP, ClassifierMixin):
    """ Multilayer Perceptron Classifier.

    Uses a neural network with one hidden layer.


    Parameters
    ----------


    Attributes
    ----------

    Notes
    -----


    References
    ----------"""
    def __init__(self, n_hidden=200, lr=0.1, l2decay=0, loss='cross_entropy',
            output_layer='softmax', batch_size=100, verbose=0):
        super(MLPClassifierA, self).__init__(n_hidden, lr, l2decay, loss,
                output_layer, batch_size, verbose)

    def fit(self, X, y, max_epochs=10, shuffle_data=False):
        self.lb = LabelBinarizer()
        one_hot_labels = self.lb.fit_transform(y)
        super(MLPClassifierA, self).fit(
                X, one_hot_labels, max_epochs,
                shuffle_data)
        return self

    def predict(self, X):
        prediction = super(MLPClassifierA, self).predict(X)
        return self.lb.inverse_transform(prediction)
