"""
Testing for export functions of decision trees (sklearn.tree.export).
"""

from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_graphviz
from sklearn.externals.six import StringIO

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]


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
    contents2 = "digraph Tree {\n" \
                "0 [label=\"X[0] <= 0.0000\\ngini = 0.5\\n" \
                "samples = 6\", shape=\"box\"] ;\n" \
                "1 [label=\"gini = 0.0000\\nsamples = 3\\n" \
                "value = [ 3.  0.]\", shape=\"box\"] ;\n" \
                "0 -> 1 ;\n" \
                "2 [label=\"gini = 0.0000\\nsamples = 3\\n" \
                "value = [ 0.  3.]\", shape=\"box\"] ;\n" \
                "0 -> 2 ;\n" \
                "}"

    assert_equal(contents1, contents2)

    # Test with feature_names
    out = StringIO()
    export_graphviz(clf, out_file=out, feature_names=["feature0", "feature1"])
    contents1 = out.getvalue()
    contents2 = "digraph Tree {\n" \
                "0 [label=\"feature0 <= 0.0000\\ngini = 0.5\\n" \
                "samples = 6\", shape=\"box\"] ;\n" \
                "1 [label=\"gini = 0.0000\\nsamples = 3\\n" \
                "value = [ 3.  0.]\", shape=\"box\"] ;\n" \
                "0 -> 1 ;\n" \
                "2 [label=\"gini = 0.0000\\nsamples = 3\\n" \
                "value = [ 0.  3.]\", shape=\"box\"] ;\n" \
                "0 -> 2 ;\n" \
                "}"

    assert_equal(contents1, contents2)

    # Test max_depth
    out = StringIO()
    export_graphviz(clf, out_file=out, max_depth=0)
    contents1 = out.getvalue()
    contents2 = "digraph Tree {\n" \
                "0 [label=\"X[0] <= 0.0000\\ngini = 0.5\\n" \
                "samples = 6\", shape=\"box\"] ;\n" \
                "1 [label=\"(...)\", shape=\"box\"] ;\n" \
                "0 -> 1 ;\n" \
                "2 [label=\"(...)\", shape=\"box\"] ;\n" \
                "0 -> 2 ;\n" \
                "}"

    assert_equal(contents1, contents2)


def test_graphviz_errors():
    # Check for errors of export_graphviz
    clf = DecisionTreeClassifier(max_depth=3, min_samples_split=1)
    clf.fit(X, y)

    out = StringIO()
    assert_raises(IndexError, export_graphviz, clf, out, feature_names=[])


if __name__ == "__main__":
    import nose
    nose.runmodule()
