import pytest

from sklearn.utils import all_estimators
from sklearn.utils._typing import get_docstring_annotations


TYPING_IGNORED = {
    'ARDRegression', 'AdaBoostClassifier', 'AdaBoostRegressor',
    'AdditiveChi2Sampler', 'AffinityPropagation',
    'AgglomerativeClustering', 'BaggingClassifier', 'BaggingRegressor',
    'BayesianGaussianMixture', 'BayesianRidge', 'BernoulliNB',
    'BernoulliRBM', 'Binarizer', 'Birch', 'CCA', 'CalibratedClassifierCV',
    'CategoricalNB', 'ClassifierChain', 'ColumnTransformer',
    'ComplementNB', 'CountVectorizer', 'DBSCAN', 'DecisionTreeClassifier',
    'DecisionTreeRegressor', 'DictVectorizer', 'DictionaryLearning',
    'DummyClassifier', 'DummyRegressor', 'ElasticNet', 'ElasticNetCV',
    'EllipticEnvelope', 'EmpiricalCovariance', 'ExtraTreeClassifier',
    'ExtraTreeRegressor', 'ExtraTreesClassifier', 'ExtraTreesRegressor',
    'FactorAnalysis', 'FastICA', 'FeatureAgglomeration', 'FeatureHasher',
    'FeatureUnion', 'FunctionTransformer', 'GammaRegressor',
    'GaussianMixture', 'GaussianNB', 'GaussianProcessClassifier',
    'GaussianProcessRegressor', 'GaussianRandomProjection',
    'GenericUnivariateSelect', 'GradientBoostingClassifier',
    'GradientBoostingRegressor', 'GraphicalLasso', 'GraphicalLassoCV',
    'GridSearchCV', 'HashingVectorizer', 'HistGradientBoostingClassifier',
    'HistGradientBoostingRegressor', 'HuberRegressor', 'IncrementalPCA',
    'IsolationForest', 'Isomap', 'IsotonicRegression', 'IterativeImputer',
    'KBinsDiscretizer', 'KMeans', 'KNNImputer', 'KNeighborsClassifier',
    'KNeighborsRegressor', 'KNeighborsTransformer', 'KernelCenterer',
    'KernelDensity', 'KernelPCA', 'KernelRidge', 'LabelBinarizer',
    'LabelEncoder', 'LabelPropagation', 'LabelSpreading', 'Lars', 'LarsCV',
    'Lasso', 'LassoCV', 'LassoLars', 'LassoLarsCV', 'LassoLarsIC',
    'LatentDirichletAllocation', 'LedoitWolf',
    'LinearDiscriminantAnalysis', 'LinearRegression', 'LinearSVC',
    'LinearSVR', 'LocalOutlierFactor', 'LocallyLinearEmbedding',
    'LogisticRegressionCV', 'MDS', 'MLPClassifier',
    'MLPRegressor', 'MaxAbsScaler', 'MeanShift', 'MinCovDet',
    'MinMaxScaler', 'MiniBatchDictionaryLearning', 'MiniBatchKMeans',
    'MiniBatchSparsePCA', 'MissingIndicator', 'MultiLabelBinarizer',
    'MultiOutputClassifier', 'MultiOutputRegressor', 'MultiTaskElasticNet',
    'MultiTaskElasticNetCV', 'MultiTaskLasso', 'MultiTaskLassoCV',
    'MultinomialNB', 'NMF', 'NearestCentroid', 'NearestNeighbors',
    'NeighborhoodComponentsAnalysis', 'Normalizer', 'NuSVC', 'NuSVR',
    'Nystroem', 'OAS', 'OPTICS', 'OneClassSVM', 'OneHotEncoder',
    'OneVsOneClassifier', 'OneVsRestClassifier', 'OrdinalEncoder',
    'OrthogonalMatchingPursuit', 'OrthogonalMatchingPursuitCV',
    'OutputCodeClassifier', 'PCA', 'PLSCanonical', 'PLSRegression',
    'PLSSVD', 'PassiveAggressiveClassifier', 'PassiveAggressiveRegressor',
    'PatchExtractor', 'Perceptron', 'Pipeline', 'PoissonRegressor',
    'PolynomialFeatures', 'PowerTransformer',
    'QuadraticDiscriminantAnalysis', 'QuantileTransformer',
    'RANSACRegressor', 'RBFSampler', 'RFE', 'RFECV',
    'RadiusNeighborsClassifier', 'RadiusNeighborsRegressor',
    'RadiusNeighborsTransformer', 'RandomForestClassifier',
    'RandomForestRegressor', 'RandomTreesEmbedding', 'RandomizedSearchCV',
    'RegressorChain', 'Ridge', 'RidgeCV', 'RidgeClassifier',
    'RidgeClassifierCV', 'RobustScaler', 'SGDClassifier', 'SGDRegressor',
    'SVC', 'SVR', 'SelectFdr', 'SelectFpr', 'SelectFromModel', 'SelectFwe',
    'SelectKBest', 'SelectPercentile', 'ShrunkCovariance', 'SimpleImputer',
    'SkewedChi2Sampler', 'SparseCoder', 'SparsePCA',
    'SparseRandomProjection', 'SpectralBiclustering', 'SpectralClustering',
    'SpectralCoclustering', 'SpectralEmbedding', 'StackingClassifier',
    'StackingRegressor', 'StandardScaler', 'TSNE', 'TfidfTransformer',
    'TfidfVectorizer', 'TheilSenRegressor', 'TransformedTargetRegressor',
    'TruncatedSVD', 'TweedieRegressor', 'VarianceThreshold',
    'VotingClassifier', 'VotingRegressor'
}


@pytest.mark.parametrize(
    'name, Estimator', [
        pytest.param(
            name, Estimator, marks=pytest.mark.skipif(
                name in TYPING_IGNORED,
                reason="Estimator does not have annotations"))
        for name, Estimator in all_estimators()])
def test_estimators_typestring(name, Estimator):
    # Check that docstring's type is formated correctly
    docscrape = pytest.importorskip('numpydoc.docscrape')

    doc = docscrape.ClassDoc(Estimator)
    parameters = doc['Parameters']
    parameter_annnotations = get_docstring_annotations(Estimator.__init__)
    _check_annotations(parameters, parameter_annnotations)

    attributes = doc['Attributes']
    attribute_annotations = get_docstring_annotations(Estimator)
    _check_annotations(attributes, attribute_annotations)


def _check_annotations(docstring_items, expected_annotations):

    assert len(docstring_items) == len(expected_annotations)

    for item in docstring_items:
        name, type_str = item.name, item.type

        # skip annotations with "shape of" for now, this can be added when
        # we support Annotated
        if "of shape" in type_str:
            continue

        # whitespaces are collapsed to one whitespace
        type_str = ' '.join(item.type.split())
        assert type_str.startswith(expected_annotations[name]), (
            f"{name} has incorrectly formated docstring")
