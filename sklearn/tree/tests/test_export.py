"""
Testing for export functions of decision trees (sklearn.tree.export).
"""

from re import finditer

from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import export_graphviz
from sklearn.externals.six import StringIO
from sklearn.utils.testing import assert_in

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
y2 = [[-1, 1], [-1, 1], [-1, 1], [1, 2], [1, 2], [1, 3]]
w = [1, 1, 1, .5, .5, .5]


def test_graphviz_toy():
    # Check correctness of export_graphviz
    clf = DecisionTreeClassifier(max_depth=3,
                                 min_samples_split=1,
                                 criterion="gini",
                                 random_state=2)
    clf.fit(X, y)

    # Test export code
    out = StringIO()
    export_graphviz(clf, out_file=out)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box] ;\n' \
                '0 [label="X[0] <= 0.0\\ngini = 0.5\\nsamples = 6\\n' \
                'value = [3, 3]"] ;\n' \
                '1 [label="gini = 0.0\\nsamples = 3\\nvalue = [3, 0]"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="True"] ;\n' \
                '2 [label="gini = 0.0\\nsamples = 3\\nvalue = [0, 3]"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="False"] ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test with feature_names
    out = StringIO()
    export_graphviz(clf, out_file=out, feature_names=["feature0", "feature1"])
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box] ;\n' \
                '0 [label="feature0 <= 0.0\\ngini = 0.5\\nsamples = 6\\n' \
                'value = [3, 3]"] ;\n' \
                '1 [label="gini = 0.0\\nsamples = 3\\nvalue = [3, 0]"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="True"] ;\n' \
                '2 [label="gini = 0.0\\nsamples = 3\\nvalue = [0, 3]"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="False"] ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test with class_names
    out = StringIO()
    export_graphviz(clf, out_file=out, class_names=["yes", "no"])
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box] ;\n' \
                '0 [label="X[0] <= 0.0\\ngini = 0.5\\nsamples = 6\\n' \
                'value = [3, 3]\\nclass = yes"] ;\n' \
                '1 [label="gini = 0.0\\nsamples = 3\\nvalue = [3, 0]\\n' \
                'class = yes"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="True"] ;\n' \
                '2 [label="gini = 0.0\\nsamples = 3\\nvalue = [0, 3]\\n' \
                'class = no"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="False"] ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test plot_options
    out = StringIO()
    export_graphviz(clf, out_file=out, filled=True, impurity=False,
                    proportion=True, special_characters=True, rounded=True)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box, style="filled, rounded", color="black", ' \
                'fontname=helvetica] ;\n' \
                'edge [fontname=helvetica] ;\n' \
                '0 [label=<X<SUB>0</SUB> &le; 0.0<br/>samples = 100.0%<br/>' \
                'value = [0.5, 0.5]>, fillcolor="#e5813900"] ;\n' \
                '1 [label=<samples = 50.0%<br/>value = [1.0, 0.0]>, ' \
                'fillcolor="#e58139ff"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="True"] ;\n' \
                '2 [label=<samples = 50.0%<br/>value = [0.0, 1.0]>, ' \
                'fillcolor="#399de5ff"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="False"] ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test max_depth
    out = StringIO()
    export_graphviz(clf, out_file=out, max_depth=0, class_names=True)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box] ;\n' \
                '0 [label="X[0] <= 0.0\\ngini = 0.5\\nsamples = 6\\n' \
                'value = [3, 3]\\nclass = y[0]"] ;\n' \
                '1 [label="(...)"] ;\n' \
                '0 -> 1 ;\n' \
                '2 [label="(...)"] ;\n' \
                '0 -> 2 ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test max_depth with plot_options
    out = StringIO()
    export_graphviz(clf, out_file=out, max_depth=0, filled=True,
                    node_ids=True)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box, style="filled", color="black"] ;\n' \
                '0 [label="node #0\\nX[0] <= 0.0\\ngini = 0.5\\n' \
                'samples = 6\\nvalue = [3, 3]", fillcolor="#e5813900"] ;\n' \
                '1 [label="(...)", fillcolor="#C0C0C0"] ;\n' \
                '0 -> 1 ;\n' \
                '2 [label="(...)", fillcolor="#C0C0C0"] ;\n' \
                '0 -> 2 ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test multi-output with weighted samples
    clf = DecisionTreeClassifier(max_depth=2,
                                 min_samples_split=1,
                                 criterion="gini",
                                 random_state=2)
    clf = clf.fit(X, y2, sample_weight=w)

    out = StringIO()
    export_graphviz(clf, out_file=out, filled=True, impurity=False)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box, style="filled", color="black"] ;\n' \
                '0 [label="X[0] <= 0.0\\nsamples = 6\\n' \
                'value = [[3.0, 1.5, 0.0]\\n' \
                '[3.0, 1.0, 0.5]]", fillcolor="#e5813900"] ;\n' \
                '1 [label="samples = 3\\nvalue = [[3, 0, 0]\\n' \
                '[3, 0, 0]]", fillcolor="#e58139ff"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="True"] ;\n' \
                '2 [label="X[0] <= 1.5\\nsamples = 3\\n' \
                'value = [[0.0, 1.5, 0.0]\\n' \
                '[0.0, 1.0, 0.5]]", fillcolor="#e5813986"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="False"] ;\n' \
                '3 [label="samples = 2\\nvalue = [[0, 1, 0]\\n' \
                '[0, 1, 0]]", fillcolor="#e58139ff"] ;\n' \
                '2 -> 3 ;\n' \
                '4 [label="samples = 1\\nvalue = [[0.0, 0.5, 0.0]\\n' \
                '[0.0, 0.0, 0.5]]", fillcolor="#e58139ff"] ;\n' \
                '2 -> 4 ;\n' \
                '}'

    assert_equal(contents1, contents2)

    # Test regression output with plot_options
    clf = DecisionTreeRegressor(max_depth=3,
                                min_samples_split=1,
                                criterion="mse",
                                random_state=2)
    clf.fit(X, y)

    out = StringIO()
    export_graphviz(clf, out_file=out, filled=True, leaves_parallel=True,
                    rotate=True, rounded=True)
    contents1 = out.getvalue()
    contents2 = 'digraph Tree {\n' \
                'node [shape=box, style="filled, rounded", color="black", ' \
                'fontname=helvetica] ;\n' \
                'graph [ranksep=equally, splines=polyline] ;\n' \
                'edge [fontname=helvetica] ;\n' \
                'rankdir=LR ;\n' \
                '0 [label="X[0] <= 0.0\\nmse = 1.0\\nsamples = 6\\n' \
                'value = 0.0", fillcolor="#e5813980"] ;\n' \
                '1 [label="mse = 0.0\\nsamples = 3\\nvalue = -1.0", ' \
                'fillcolor="#e5813900"] ;\n' \
                '0 -> 1 [labeldistance=2.5, labelangle=-45, ' \
                'headlabel="True"] ;\n' \
                '2 [label="mse = 0.0\\nsamples = 3\\nvalue = 1.0", ' \
                'fillcolor="#e58139ff"] ;\n' \
                '0 -> 2 [labeldistance=2.5, labelangle=45, ' \
                'headlabel="False"] ;\n' \
                '{rank=same ; 0} ;\n' \
                '{rank=same ; 1; 2} ;\n' \
                '}'

    assert_equal(contents1, contents2)


def test_graphviz_errors():
    # Check for errors of export_graphviz
    clf = DecisionTreeClassifier(max_depth=3, min_samples_split=1)
    clf.fit(X, y)

    # Check feature_names error
    out = StringIO()
    assert_raises(IndexError, export_graphviz, clf, out, feature_names=[])

    # Check class_names error
    out = StringIO()
    assert_raises(IndexError, export_graphviz, clf, out, class_names=[])


def test_friedman_mse_in_graphviz():
    clf = DecisionTreeRegressor(criterion="friedman_mse", random_state=0)
    clf.fit(X, y)
    dot_data = StringIO()
    export_graphviz(clf, out_file=dot_data)

    clf = GradientBoostingClassifier(n_estimators=2, random_state=0)
    clf.fit(X, y)
    for estimator in clf.estimators_:
        export_graphviz(estimator[0], out_file=dot_data)

    for finding in finditer("\[.*?samples.*?\]", dot_data.getvalue()):
        assert_in("friedman_mse", finding.group())
