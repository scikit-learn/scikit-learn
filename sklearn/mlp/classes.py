from ..base import ClassifierMixin
from .base import BaseMLP
from ..preprocessing import LabelBinarizer


class MLPClassifier(BaseMLP, ClassifierMixin):
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
    def __init__(self,
                 n_hidden,
                 lr,
                 batch_size,
                 loss_function='squared',
                 output_function='tanh',
                 learning_method='backprop',
                 verbose=0):
        super(MLPClassifier, self).__init__(n_hidden,
                                            lr,
                                            batch_size,
                                            loss_function,
                                            output_function,
                                            learning_method,
                                            verbose=0)

    def fit(self, X, y, max_epochs=10, shuffle_data=False):
        self.lb = LabelBinarizer()
        one_hot_labels = self.lb.fit_transform(y)
        super(MLPClassifier, self).fit(
                X, one_hot_labels, max_epochs,
                shuffle_data)
        return self

    def predict(self, X):
        prediction = super(MLPClassifier, self).predict(X)
        return self.lb.inverse_transform(prediction)
