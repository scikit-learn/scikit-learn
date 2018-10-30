from sklearn import datasets # already imported in the test_classification file
from sklearn.model_selection import GridSearchCV
from sklearn.mixture import GaussianMixture
from sklearn.utils.testing import assert_raises
from sklearn.tree import DecisionTreeClassifier

def test_check_estimator_one_class():
    # test with a GMM model
    X,y = datasets.make_classification(n_samples = 10000, n_features=10, n_classes=2)
    gmm_model = GaussianMixture()
    # grid search
    param_grid = {'n_components' : [1,2,3,4],
                  'covariance_type': ['tied','full','spherical']}
    grid_search = GridSearchCV(gmm_model, param_grid, scoring='roc_auc')
    # Fit Grid Search with this dataset
    grid_search.fit(X,y)

    # Test with a DecisionTreeClassifier() with a one class dataset
    X2,y2 = datasets.make_classification(n_samples = 10000, n_features=10, n_classes=1)
    dtc_model = DecisionTreeClassifier()
    parameters = {'max_depth':[5,10,15]}
    grid_search2 = GridSearchCV(dtc_model, parameters, scoring='roc_auc')
    grid_search2.fit(X2,y2)
