from contextlib import closing
from io import StringIO

import pytest

from sklearn import config_context
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
from sklearn.multiclass import OneVsOneClassifier
from sklearn._display_estimator import _write_label_html
from sklearn._display_estimator import _estimator_details
from sklearn._display_estimator import _type_of_html_estimator
from sklearn._display_estimator import _estimator_repr_html


@pytest.mark.parametrize('est, expected', [
    ('None', 'None'),
    ('passthrough', 'passthrough'),
    ('hello\nworld', 'hello&#xa;world')
])
def test_estimator_tool_tip(est, expected):
    assert expected == _estimator_details(est)


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
def test_type_of_html_estimator_single_str_none(est):
    est_html_info = _type_of_html_estimator(est)
    assert est_html_info.type == 'single'
    assert est_html_info.estimators[0] == est
    assert est_html_info.names[0] == str(est)
    assert est_html_info.name_details[0] == str(est)


def test_type_of_html_estimator_single_estimator():
    # single estimator prints all the details
    est = LogisticRegression(C=10.0)
    est_html_info = _type_of_html_estimator(est, first_call=True)
    assert est_html_info.type == 'single'
    assert est_html_info.estimators[0] == est
    assert est_html_info.names[0] == est.__class__.__name__
    assert est_html_info.name_details[0] == _estimator_details(est)


def test_type_of_html_estimator_pipeline():
    # multiple estimators in a pipeline prints only the changes
    pipe = Pipeline([
        ('imputer', SimpleImputer()),
        ('classifier', LogisticRegression())
    ])
    est_html_info = _type_of_html_estimator(pipe)
    assert est_html_info.type == 'serial'
    assert est_html_info.estimators == [step[1] for step in pipe.steps]
    assert est_html_info.names == ['imputer', 'classifier']

    with config_context(print_changed_only=True):
        assert est_html_info.name_details == [_estimator_details(step[1])
                                              for step in pipe.steps]


def test_type_of_html_estimator_feature_union():
    f_union = FeatureUnion([
        ('pca', PCA()), ('svd', TruncatedSVD())
    ])
    est_html_info = _type_of_html_estimator(f_union)
    assert est_html_info.type == 'parallel'
    assert est_html_info.names == ['pca', 'svd']
    assert est_html_info.estimators == [trans[1]
                                        for trans in f_union.transformer_list]
    assert est_html_info.name_details == [None, None]


def test_type_of_html_estimator_voting():
    clf = VotingClassifier([
        ('log_reg', LogisticRegression()),
        ('mlp', MLPClassifier())
    ])
    est_html_info = _type_of_html_estimator(clf)
    assert est_html_info.type == 'parallel'
    assert est_html_info.estimators == [trans[1]
                                        for trans in clf.estimators]
    assert est_html_info.names == ['log_reg', 'mlp']
    assert est_html_info.name_details == [None, None]


def test_type_of_html_estimator_column_transformer():
    ct = ColumnTransformer([
        ('pca', PCA(), ['num1', 'num2']),
        ('svd', TruncatedSVD, [0, 3])
    ])
    est_html_info = _type_of_html_estimator(ct)
    assert est_html_info.type == 'parallel'
    assert est_html_info.estimators == [trans[1]
                                        for trans in ct.transformers]
    assert est_html_info.names == ['pca', 'svd']
    assert est_html_info.name_details == [['num1', 'num2'], [0, 3]]


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
      'preprocessor</label>',
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
