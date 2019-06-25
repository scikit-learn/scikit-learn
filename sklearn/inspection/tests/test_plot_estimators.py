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
from sklearn.inspection._plot_estimators import _write_label_html
from sklearn.inspection._plot_estimators import _estimator_tool_tip
from sklearn.inspection._plot_estimators import _type_of_html_estimator
from sklearn.inspection._plot_estimators import export_html
from sklearn.inspection._plot_estimators import _STYLE


@pytest.mark.parametrize('est, expected', [
    ('None', 'None'),
    ('passthrough', 'passthrough'),
    ('hello\nworld', 'hello&#xa;world')
])
def test_estimator_tool_tip(est, expected):
    assert expected == _estimator_tool_tip(est)


def test_write_label_html():
    name = "LogisticRegression"
    tool_tip = "hello-world"

    expected = ('<div class="sk-label" '
                'sk-data-tooltip="hello-world">'
                'LogisticRegression</div>')

    with closing(StringIO()) as out:
        _write_label_html(out, name, tool_tip)
        html_label = out.getvalue()
        assert html_label == expected


def test_type_of_html_estimator_error():
    with pytest.raises(ValueError, match="Invalid estimator"):
        _type_of_html_estimator(100)


@pytest.mark.parametrize('est', ['passthrough', 'drop', None])
def test_type_of_html_estimator_single_str_none(est):
    est_html_info = _type_of_html_estimator(est)
    assert est_html_info.type == 'single'
    assert est_html_info.estimators == est
    assert est_html_info.names == str(est)
    assert est_html_info.name_tips == str(est)


def test_type_of_html_estimator_single_estimator():
    est = LogisticRegression(C=10.0)
    est_html_info = _type_of_html_estimator(est)
    assert est_html_info.type == 'single'
    assert est_html_info.estimators == est
    assert est_html_info.names == est.__class__.__name__
    assert est_html_info.name_tips == _estimator_tool_tip(est)


def test_type_of_html_estimator_pipeline():
    pipe = Pipeline([
        ('imputer', SimpleImputer()),
        ('classifier', LogisticRegression())
    ])
    est_html_info = _type_of_html_estimator(pipe)
    assert est_html_info.type == 'serial'
    assert est_html_info.estimators == [step[1] for step in pipe.steps]
    assert est_html_info.names == ['imputer', 'classifier']
    assert est_html_info.name_tips == [_estimator_tool_tip(step[1])
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
    assert est_html_info.name_tips == [_estimator_tool_tip(trans[1])
                                       for trans in f_union.transformer_list]


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
    assert est_html_info.name_tips == [_estimator_tool_tip(trans[1])
                                       for trans in clf.estimators]


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
    assert est_html_info.name_tips == [['num1', 'num2'], [0, 3]]


expected_export_html = """<html><head><style>
{style}
</style></head>
<body>
<div class="sk-container">
  <div class="sk-serial">
    <div class="sk-serial-item sk-dashed-wrapped">
      <div
        class="sk-label"
        sk-data-tooltip="ColumnTransformer(transformers=[(\'num\',&#xa;
        Pipeline(steps=[(\'pass\', \'passthrough\'),&#xa;
        (\'imputer\',&#xa;
        SimpleImputer(strategy=\'median\'))]),&#xa;
        [\'a\', \'b\', \'c\', \'d\', \'e\']),&#xa;
        (\'cat\',&#xa;
        Pipeline(steps=[(\'imputer\',&#xa;
        SimpleImputer(missing_values=\'empty\',&#xa;
        strategy=\'constant\')),&#xa;
        (\'one-hot\', OneHotEncoder())]),&#xa;[0, 1, 2, 3])])"
      >
        preprocessor
      </div>
      <div class="sk-parallel">
        <div class="sk-parallel-item">
          <div
            class="sk-label"
            sk-data-tooltip="[\'a\', \'b\', \'c\', \'d\', \'e\']"
          >
            num
          </div>
          <div class="sk-serial">
            <div class="sk-serial">
              <div class="sk-serial-item">
                <div class="sk-estimator" sk-data-tooltip="passthrough">
                  passthrough
                </div>
              </div>
              <div class="sk-serial-item">
                <div
                  class="sk-estimator"
                  sk-data-tooltip="SimpleImputer(strategy=\'median\')"
                >
                  SimpleImputer
                </div>
              </div>
            </div>
          </div>
        </div>
        <div class="sk-parallel-item">
          <div class="sk-label" sk-data-tooltip="[0, 1, 2, 3]">cat</div>
          <div class="sk-serial">
            <div class="sk-serial">
              <div class="sk-serial-item">
                <div
                  class="sk-estimator"
                  sk-data-tooltip="SimpleImputer(missing_values=\'empty\',
                  strategy=\'constant\')"
                >
                  SimpleImputer
                </div>
              </div>
              <div class="sk-serial-item">
                <div class="sk-estimator" sk-data-tooltip="OneHotEncoder()">
                  OneHotEncoder
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
    <div class="sk-serial-item sk-dashed-wrapped">
      <div
        class="sk-label"
        sk-data-tooltip="FeatureUnion(transformer_list=[(\'pca\',
        PCA(n_components=1)),&#xa;(\'tsvd\',&#xa;
        Pipeline(steps=[(\'first\',&#xa;
        TruncatedSVD(n_components=3)),&#xa;
        (\'select\',&#xa;
        SelectPercentile())]))])"
      >
        feat_u
      </div>
      <div class="sk-parallel">
        <div class="sk-parallel-item">
          <div class="sk-label" sk-data-tooltip="PCA(n_components=1)">
            pca
          </div>
          <div class="sk-serial">
            <div class="sk-serial-item">
              <div class="sk-estimator" sk-data-tooltip="PCA(n_components=1)">
                PCA
              </div>
            </div>
          </div>
        </div>
        <div class="sk-parallel-item">
          <div
            class="sk-label"
            sk-data-tooltip="Pipeline(steps=[(\'first\',
            TruncatedSVD(n_components=3)),&#xa;
            (\'select\', SelectPercentile())])"
          >
            tsvd
          </div>
          <div class="sk-serial">
            <div class="sk-serial">
              <div class="sk-serial-item">
                <div
                  class="sk-estimator"
                  sk-data-tooltip="TruncatedSVD(n_components=3)"
                >
                  TruncatedSVD
                </div>
              </div>
              <div class="sk-serial-item">
                <div class="sk-estimator"
                     sk-data-tooltip="SelectPercentile()">
                  SelectPercentile
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
    <div class="sk-serial-item sk-dashed-wrapped">
      <div
        class="sk-label"
        sk-data-tooltip="VotingClassifier(estimators=[(\'lr\',
        LogisticRegression(random_state=1)),&#xa;
        (\'mlp\', MLPClassifier(alpha=0.001))])"
      >
        classifier
      </div>
      <div class="sk-parallel">
        <div class="sk-parallel-item">
          <div
            class="sk-label"
            sk-data-tooltip="LogisticRegression(random_state=1)"
          >
            lr
          </div>
          <div class="sk-serial">
            <div class="sk-serial-item">
              <div
                class="sk-estimator"
                sk-data-tooltip="LogisticRegression(random_state=1)"
              >
                LogisticRegression
              </div>
            </div>
          </div>
        </div>
        <div class="sk-parallel-item">
          <div class="sk-label" sk-data-tooltip="MLPClassifier(alpha=0.001)">
            mlp
          </div>
          <div class="sk-serial">
            <div class="sk-serial-item">
              <div
                class="sk-estimator"
                sk-data-tooltip="MLPClassifier(alpha=0.001)"
              >
                MLPClassifier
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</div>
<div class="sk-final-spacer">
  (\'classifier\', VotingClassifier(estimators=[(\'lr\',
  LogisticRegression(random_state=1)),&#xa; (\'mlp\',
  MLPClassifier(alpha=0.001))]))
</div>
</body></html>
""".format(style=_STYLE).replace('\n', '').replace(' ', '')


def test_export_html():
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
    html_output = export_html(pipe)

    # IPython HTML
    if hasattr(html_output, 'data'):
        html_output = html_output.data
    assert expected_export_html == html_output.replace(' ', '')
