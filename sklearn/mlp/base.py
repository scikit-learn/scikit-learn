import numpy as np

from ..base import BaseEstimator
from .mlp_fast import sgd, predict, SquaredLoss, SoftMax, MultiCrossEntropyLoss, CrossEntropyLoss, LogSig
from .mlp_fast import Tanh


class BaseMLP(BaseEstimator):

    def __init__(self,
                 n_hidden,
                 lr,
                 batch_size,
                 loss_function='squared',
                 output_function='tanh',
                 learning_method='backprop',
                 shuffle_data=False,
                 verbose=0):

        self.n_hidden = n_hidden
        self.learning_method = learning_method
        self.lr = lr
        self.batch_size = batch_size
        self.shuffle_data = shuffle_data

        if loss_function == 'squared':
            self.loss_function = SquaredLoss()
        elif loss_function == 'cross-entropy' and output_function == 'softmax':
            self.loss_function = MultiCrossEntropyLoss()
        else:
            # TODO
            raise ValueError('Loss function must be one of ...')

        if output_function == 'tanh':
            self.output_function = Tanh()
        elif output_function == 'softmax':
            self.output_function = SoftMax()
        else:
            # TODO
            raise ValueError('Output function must be one of ...')

    def fit(self, X, y, max_epochs, verbose=0):

        n_samples, n_features = X.shape
        n_outs = y.shape[1]

        # Only one class
        if n_outs == 1 and isinstance(self.loss_function, MultiCrossEntropyLoss):
            self.loss_function = CrossEntropyLoss()
            self.output_function = LogSig()

        self.weights_hidden = np.zeros((n_features, self.n_hidden))
        self.bias_hidden = np.zeros(self.n_hidden)
        self.weights_output = np.zeros((self.n_hidden, n_outs))
        self.bias_output = np.zeros(n_outs)

        sgd(X, y,
            loss=self.loss_function,
            output=self.output_function,
            hidden=Tanh(),
            weights_hidden=self.weights_hidden,
            weights_output=self.weights_output,
            bias_hidden=self.bias_hidden,
            bias_output=self.bias_output,
            lr=self.lr,
            n_hidden=self.n_hidden,
            max_epochs=max_epochs,
            batch_size=self.batch_size,
            shuffle_data=self.shuffle_data)

        return self

    def predict(self, X):
        return predict(X,
                       self.weights_hidden,
                       self.bias_hidden,
                       self.weights_output,
                       self.bias_output,
                       self.output_function,
                       Tanh())
