from ..coordinate_descent import Lasso
from ..coordinate_descent import ElasticNet
from ...utils import deprecated


@deprecated("""to be removed in v0.13;
use sklearn.linear_model.ElasticNet instead""")
class ElasticNet(ElasticNet):
    pass


@deprecated("""to be removed in v0.13;
use sklearn.linear_model.Lasso instead""")
class Lasso(Lasso):
    pass
