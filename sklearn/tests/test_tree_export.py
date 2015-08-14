# Author: Frank Zalkow
# License: BSD 3 clause

from sklearn.utils.testing import assert_in

def test_friedman_mse_in_graphviz():

    import re
    from sklearn.tree import DecisionTreeRegressor, export_graphviz
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.externals.six import StringIO

    clf = DecisionTreeRegressor(criterion="friedman_mse")
    clf.fit([[0,0], [1,1]], [0,1])
    dot_data = StringIO()
    export_graphviz(clf, out_file=dot_data)

    clf = GradientBoostingClassifier(n_estimators=2)
    clf.fit([[0,0], [1,1]], [0,1])
    for estimator in clf.estimators_:
        export_graphviz(estimator[0], out_file=dot_data)

    for finding in re.finditer("\[.+?\] ;", dot_data.getvalue()):
        assert_in("friedman_mse", finding.group())

        