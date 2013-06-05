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
                 lr_moment,
                 batch_size,
                 loss_function='cross-entropy',
                 output_function='softmax',
                 learning_method='backprop',
                 shuffle_data=False,
                 verbose=0):
        super(MLPClassifier, self).__init__(n_hidden,
                                            lr,
                                            lr_moment,
                                            batch_size,
                                            loss_function,
                                            output_function,
                                            learning_method,
                                            shuffle_data,
                                            verbose=0)

    def fit(self, X, y, max_epochs=10):
        self.lb = LabelBinarizer()
        one_hot_labels = self.lb.fit_transform(y)
        super(MLPClassifier, self).fit(
                X, one_hot_labels, max_epochs)
        return self

    def predict(self, X):
        prediction = super(MLPClassifier, self).predict(X)
        return self.lb.inverse_transform(prediction)
