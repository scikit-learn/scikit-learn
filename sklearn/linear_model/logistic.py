import numpy as np

from ..base import ClassifierMixin
from ..linear_model.base import CoefSelectTransformerMixin
from ..svm.base import BaseLibLinear
from ..svm import liblinear


class LogisticRegression(BaseLibLinear, ClassifierMixin,
                         CoefSelectTransformerMixin):
    """Logistic Regression classifier.

    In the multiclass case, the training algorithm uses a one-vs.-all (OvA)
    scheme, rather than the "true" multinomial LR (aka maximum entropy/MaxEnt).
    This class implements L1 and L2 regularized logistic regression using the
    `liblinear` library.

    Parameters
    ----------

    penalty : string, 'l1' or 'l2'
        Used to specify the norm used in the penalization

    dual : boolean
        Dual or primal formulation. Dual formulation is only
        implemented for l2 penalty. Prefer dual=False when
        n_samples > n_features.

    C : float
        Specifies the strength of the regularization. The smaller it is
        the bigger in the regularization.

    fit_intercept : bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function

    intercept_scaling : float, default: 1
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    tol: float, optional
         tolerance for stopping criteria

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independant
        of the number of samples.

    Attributes
    ----------

    `coef_` : array, shape = [n_classes-1, n_features]
        Coefficient of the features in the decision function.

        `coef_` is readonly property derived from `raw_coef_` that \
        follows the internal memory layout of liblinear.

    `intercept_` : array, shape = [n_classes-1]
        intercept (a.k.a. bias) added to the decision function.
        It is available only when parameter intercept is set to True

    See also
    --------
    LinearSVC

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    **References**:

    LIBLINEAR -- A Library for Large Linear Classification
    http://www.csie.ntu.edu.tw/~cjlin/liblinear/
    """

    def __init__(self, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1,
                 scale_C=False):

        super(LogisticRegression, self).__init__(penalty=penalty,
            dual=dual, loss='lr', tol=tol, C=C,
            fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
            scale_C=scale_C)

    def predict_proba(self, X):
        """Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.
        """
        X = np.asarray(X, dtype=np.float64, order='C')
        probas = liblinear.predict_prob_wrap(X, self.raw_coef_,
                                      self._get_solver_type(),
                                      self.tol, self.C,
                                      self.class_weight_label,
                                      self.class_weight, self.label_,
                                      self._get_bias())
        return probas[:, np.argsort(self.label_)]

    def predict_log_proba(self, X):
        """Log of Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.
        """
        return np.log(self.predict_proba(X))
