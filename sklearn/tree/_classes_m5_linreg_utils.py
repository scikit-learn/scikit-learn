from abc import abstractmethod, ABCMeta

import numpy as np

from ..linear_model import LinearRegression
from ..linear_model._base import LinearModel
from ..preprocessing import StandardScaler
from ..utils.extmath import safe_sparse_dot


def linreg_model_to_text(model, feature_names=None, target_name=None,
                         precision=3, line_breaks=False):
    """
    Converts a linear regression model to a text representation.

    :param model:
    :param feature_names:
    :param target_name:
    :param precision:
    :param line_breaks: if True, each term in the sum shows in a different line
    :return:
    """
    bits = []

    # Template for numbers: we want scientific notation with a given precision
    nb_tpl = "%%0.%se" % precision

    # Handle multi-dimensional y (its second dim must be size 1, though)
    if len(model.coef_.shape) > 1:
        assert model.coef_.shape[0] == 1
        assert len(model.coef_.shape) == 2
        coefs = np.ravel(model.coef_)  # convert to 1D
        assert len(model.intercept_) == 1
        intercept = model.intercept_.item()  # extract scalar
    else:
        coefs = model.coef_           # a 1D array
        intercept = model.intercept_  # a scalar already

    # First all coefs * drivers
    for i, c in enumerate(coefs):
        var_name = ("X[%s]" % i) if feature_names is None else feature_names[i]

        if i == 0:
            # first term
            if c < 1:
                # use scientific notation
                product_text = (nb_tpl + " * %s") % (c, var_name)
            else:
                # use standard notation with precision
                c = np.round(c, precision)
                product_text = "%s * %s" % (c, var_name)
        else:
            # all the other terms: the sign should appear
            lb = '\n' if line_breaks else ""
            coef_abs = np.abs(c)
            coef_sign = '+' if np.sign(c) > 0 else '-'
            if coef_abs < 1:
                # use scientific notation
                product_text = (("%s%s " + nb_tpl + " * %s")
                                % (lb, coef_sign, coef_abs, var_name))
            else:
                # use standard notation with precision
                coef_abs = np.round(coef_abs, precision)
                product_text = ("%s%s %s * %s"
                                % (lb, coef_sign, coef_abs, var_name))

        bits.append(product_text)

    # Finally the intercept
    if len(bits) == 0:
        # intercept is the only term in the sum
        if intercept < 1:
            # use scientific notation only for small numbers (otherwise 12
            # would read 1.2e1 ... not friendly)
            constant_text = nb_tpl % intercept
        else:
            # use standard notation with precision
            i = np.round(intercept, precision)
            constant_text = "%s" % i
    else:
        # there are other terms in the sum: the sign should appear
        lb = '\n' if line_breaks else ""
        coef_abs = np.abs(intercept)
        coef_sign = '+' if np.sign(intercept) > 0 else '-'
        if coef_abs < 1:
            # use scientific notation
            constant_text = ("%s%s " + nb_tpl) % (lb, coef_sign, coef_abs)
        else:
            # use standard notation with precision
            coef_abs = np.round(coef_abs, precision)
            constant_text = "%s%s %s" % (lb, coef_sign, coef_abs)

    bits.append(constant_text)

    txt = " ".join(bits)
    if target_name is not None:
        txt = target_name + " = " + txt

    return txt


class DeNormalizableMixIn(metaclass=ABCMeta):
    """
    An abstract class that models able to de-normalize should implement.
    """
    __slots__ = ()

    @abstractmethod
    def denormalize(self,
                    x_scaler: StandardScaler = None,
                    y_scaler: StandardScaler = None
                    ):
        """
        Denormalizes the model, knowing that it was fit with the given
        x_scaler and y_scaler
        """


class DeNormalizableLinearModelMixIn(DeNormalizableMixIn, LinearModel):
    """
    A mix-in class to add 'denormalization' capability to a linear model
    """
    def denormalize(self,
                    x_scaler: StandardScaler = None,
                    y_scaler: StandardScaler = None
                    ):
        """
        De-normalizes the linear regression model.
        Before this function is executed,
            (y-y_mean)/y_scale = self.coef_.T <dot> (x-x_mean)/x_scale + self.intercept_
        so
            (y-y_mean)/y_scale = (self.coef_/x_scale).T <dot> x + (self.intercept_ - self.coef_.T <dot> x_mean/x_scale)
        that is
            (y-y_mean)/y_scale = new_coef.T <dot> x + new_intercept
        where
         * new_coef = (self.coef_/x_scale)
         * new_intercept = (self.intercept_ - (self.intercept_ - self.coef_.T <dot> x_mean/x_scale)

        Then going back to y
            y = (new_coef * y_scale).T <dot> x + (new_intercept * y_scale + y_mean)

        :param self:
        :param x_scaler:
        :param y_scaler:
        :return:
        """
        # First save old coefficients
        self.normalized_coef_ = self.coef_
        self.normalized_intercept_ = self.intercept_

        # denormalize coefficients to take into account the x normalization
        if x_scaler is not None:
            new_coef = self.coef_ / x_scaler.scale_
            new_intercept = (
                    self.intercept_ -
                    safe_sparse_dot(x_scaler.mean_, new_coef.T,
                                    dense_output=True)
            )

            self.coef_ = new_coef
            self.intercept_ = new_intercept

        # denormalize them further to take into account the y normalization
        if y_scaler is not None:
            new_coef = self.coef_ * y_scaler.scale_
            new_intercept = y_scaler.inverse_transform(
                np.atleast_1d(self.intercept_)
            )
            if np.isscalar(self.intercept_):
                new_intercept = new_intercept[0]
            self.coef_ = new_coef
            self.intercept_ = new_intercept


class DeNormalizableLinearRegression(LinearRegression,
                                     DeNormalizableLinearModelMixIn):
    """
    A Denormalizable linear regression. The old normalized coefficients are
    kept in a new field named `feature_importances_`
    """
    @property
    def feature_importances_(self):
        if hasattr(self, '_feature_importances_'):
            return self._feature_importances_
        else:
            return self.coef_

    def denormalize(self,
                    x_scaler: StandardScaler = None,
                    y_scaler: StandardScaler = None):
        """
        Denormalizes the model, and saves a copy of the old normalized
        coefficients in self._feature_importances.

        :param x_scaler:
        :param y_scaler:
        :return:
        """
        self._feature_importances_ = self.coef_
        super(DeNormalizableLinearRegression, self).denormalize(x_scaler,
                                                                y_scaler)


# For all other it should work too
# class DeNormalizableLasso(Lasso, DeNormalizableLinearModelMixIn):
#     pass
