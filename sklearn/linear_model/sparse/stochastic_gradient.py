from ...utils import deprecated
from ..stochastic_gradient import SGDClassifier as SGDClassifier_
from ..stochastic_gradient import SGDRegressor as SGDRegressor_


@deprecated("""to be removed in v0.12;
use sklearn.linear_model.SGDClassifier directly""")
class SGDClassifier(SGDClassifier_):
    pass


@deprecated("""to be removed in v0.12;
use sklearn.linear_model.SGDRegressor directly""")
class SGDRegressor(SGDRegressor_):
    pass
