from contextlib import closing
from io import StringIO

import pytest

from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.pipeline import Pipeline
from sklearn.pipeline import FeatureUnion
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import VotingClassifier
from sklearn.feature_selection import SelectPercentile
from sklearn.preprocessing import OneHotEncoder
from sklearn.svm import LinearSVC
from sklearn.svm import LinearSVR
from sklearn.tree import DecisionTreeClassifier
from sklearn.multiclass import OneVsOneClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import StackingRegressor
from sklearn._display_estimator import _write_label_html
from sklearn._display_estimator import _get_visual_block
from sklearn._display_estimator import _estimator_repr_html


@pytest.mark.parametrize("checked", [True, False])
def test_write_label_html(checked):
    name = "LogisticRegression"
    tool_tip = "hello-world"

    with closing(StringIO()) as out:
        _write_label_html(out, name, tool_tip, checked=checked)
        html_label = out.getvalue()
        assert 'LogisticRegression</label>' in html_label
        assert html_label.startswith('<div class="sk-label-container">')
        assert '<pre>hello-world</pre>' in html_label
        if checked:
            assert 'checked>' in html_label


@pytest.mark.parametrize('est', ['passthrough', 'drop', None])
def test_get_visual_block_single_str_none(est):
    est_html_info = _get_visual_block(est)
    assert est_html_info.kind == 'single'
    assert est_html_info.estimators == est
    assert est_html_info.names == str(est)
    assert est_html_info.name_details == str(est)


def test_get_visual_block_single_estimator():
    est = LogisticRegression(C=10.0)
    est_html_info = _get_visual_block(est)
    assert est_html_info.kind == 'single'
    assert est_html_info.estimators == est
    assert est_html_info.names == est.__class__.__name__
    assert est_html_info.name_details == str(est)


def test_get_visual_block_pipeline():
    pipe = Pipeline([
        ('imputer', SimpleImputer()),
        ('do_nothing', 'passthrough'),
        ('do_nothing_more', None),
        ('classifier', LogisticRegression())
    ])
    est_html_info = _get_visual_block(pipe)
    assert est_html_info.kind == 'serial'
    assert est_html_info.estimators == tuple(step[1] for step in pipe.steps)
    assert est_html_info.names == ['imputer: SimpleImputer',
                                   'do_nothing: passthrough',
                                   'do_nothing_more: passthrough',
                                   'classifier: LogisticRegression']
    assert est_html_info.name_details is None


def test_get_visual_block_feature_union():
    f_union = FeatureUnion([
        ('pca', PCA()), ('svd', TruncatedSVD())
    ])
    est_html_info = _get_visual_block(f_union)
    assert est_html_info.kind == 'parallel'
    assert est_html_info.names == ('pca', 'svd')
    assert est_html_info.estimators == tuple(
        trans[1] for trans in f_union.transformer_list)
    assert est_html_info.name_details == (None, None)


def test_get_visual_block_voting():
    clf = VotingClassifier([
        ('log_reg', LogisticRegression()),
        ('mlp', MLPClassifier())
    ])
    est_html_info = _get_visual_block(clf)
    assert est_html_info.kind == 'parallel'
    assert est_html_info.estimators == tuple(trans[1]
                                             for trans in clf.estimators)
    assert est_html_info.names == ('log_reg', 'mlp')
    assert est_html_info.name_details == (None, None)


def test_get_visual_block_column_transformer():
    ct = ColumnTransformer([
        ('pca', PCA(), ['num1', 'num2']),
        ('svd', TruncatedSVD, [0, 3])
    ])
    est_html_info = _get_visual_block(ct)
    assert est_html_info.kind == 'parallel'
    assert est_html_info.estimators == tuple(
        trans[1] for trans in ct.transformers)
    assert est_html_info.names == ('pca', 'svd')
    assert est_html_info.name_details == (['num1', 'num2'], [0, 3])


def test_display_estimator_pipeline():
    num_trans = Pipeline(steps=[
        ('pass', 'passthrough'),
        ('imputer', SimpleImputer(strategy='median'))
    ])

    cat_trans = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='constant',
                                  missing_values='empty')),
        ('one-hot', OneHotEncoder())
    ])

    preprocess = ColumnTransformer([
        ('num', num_trans, ['a', 'b', 'c', 'd', 'e']),
        ('cat', cat_trans, [0, 1, 2, 3])
    ])

    feat_u = FeatureUnion([
            ('pca', PCA(n_components=1)),
            ('tsvd', Pipeline([('first', TruncatedSVD(n_components=3)),
                               ('select', SelectPercentile())]))
    ])

    clf = VotingClassifier([
        ('lr', LogisticRegression(solver='lbfgs', random_state=1)),
        ('mlp', MLPClassifier(alpha=0.001))
    ])

    pipe = Pipeline([
        ('preprocessor', preprocess), ('feat_u', feat_u), ('classifier', clf)
    ])
    html_output = _estimator_repr_html(pipe)

    expected_strings = [
      'passthrough</label>',
      'div class=\"sk-toggleable__content\"><pre>SimpleImputer'
      '(strategy=\'median\')',
      '<pre>SimpleImputer(missing_values=\'empty\', strategy=\'constant\')'
      '</pre>',
      '(\'one-hot\', OneHotEncoder',
      'preprocessor: ColumnTransformer</label>',
      '<pre>[\'a\', \'b\', \'c\', \'d\', \'e\']</pre>',
      '<pre>LogisticRegression(random_state=1)</pre>',
      '<pre>SelectPercentile()</pre>',
      '>TruncatedSVD</label>',
      '<pre>TruncatedSVD(n_components=3)',
    ]

    for expected_string in expected_strings:
        assert expected_string in html_output


def test_display_estimator_ovo_classifier():
    ovo = OneVsOneClassifier(LinearSVC())
    html_output = _estimator_repr_html(ovo)
    assert "pre>OneVsOneClassifier(estimator=LinearSVC" in html_output
    assert "LinearSVC</label>" in html_output


@pytest.mark.parametrize("final_estimator", [None, LinearSVC()])
def test_stacking_classsifer(final_estimator):
    estimators = [('mlp', MLPClassifier(alpha=0.001)),
                  ('tree', DecisionTreeClassifier())]
    clf = StackingClassifier(
        estimators=estimators, final_estimator=final_estimator)

    html_output = _estimator_repr_html(clf)

    assert "('mlp', MLPClassifier(alpha=0.001)" in html_output
    assert "('tree', DecisionTreeClassifier()" in html_output
    if final_estimator is None:
        assert "LogisticRegression()</pre>" in html_output
    else:
        assert final_estimator.__class__.__name__ in html_output


@pytest.mark.parametrize("final_estimator", [None, LinearSVR()])
def test_stacking_regressor(final_estimator):
    reg = StackingRegressor(
        estimators=[('svr', LinearSVR())], final_estimator=final_estimator)

    html_output = _estimator_repr_html(reg)

    assert "('svr', LinearSVR()" in html_output
    print(html_output)
    if final_estimator is None:
        assert "RidgeCV</label>" in html_output
    else:
        assert final_estimator.__class__.__name__ in html_output
