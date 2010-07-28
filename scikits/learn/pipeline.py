from .base import BaseEstimator

class Pipeline(BaseEstimator):
    """
    Pipeline of transformers with a final predictor
    Sequentialy apply a list of transformers and a final predictor
    A transformer implements fit & transform methods
    A predictor implements fit & predict methods
    
    Example
    =======
    
    >>> from scikits.learn import svm, datasets
    >>> from scikits.learn.datasets import samples_generator
    >>> from scikits.learn.feature_selection import SelectKBest, f_regression
    >>> from scikits.learn.pipeline import Pipeline

    >>> # generate some data to play with
    >>> X, y = samples_generator.test_dataset_classif(k=5)

    >>> # ANOVA SVM-C
    >>> anova_filter = SelectKBest(f_regression, k=5)
    >>> clf = svm.SVC(kernel='linear')
    
    >>> anova_svm = Pipeline([anova_filter], clf)
    >>> _ = anova_svm.fit(X,y)
    
    >>> prediction = anova_svm.predict(X)
    """
    def __init__(self, transformers=[], estimator=None):
        """
        methods: list of UnivRanking objects,
        ie.: with fit/reduce/getSelectedFeatures methods
        """
        for t in  transformers:
            assert hasattr(t, "fit") and hasattr(t, "transform"), ValueError(
                "All transformers should implement fit and transform",
                "'%s' (type %s) )" % (t, type(t))
        )
        assert hasattr(estimator, "fit") and hasattr(estimator, "predict"), \
            ValueError("Predictor should implement fit and predict",
                "'%s' (type %s) )" % (t, type(t))
        )
        self.transformers = transformers
        self.estimator = estimator

    def fit(self, X, y=None):
        Xt = X
        for transformer in self.transformers:
            Xt = transformer.fit(Xt, y).transform(Xt)
        self.estimator.fit(Xt, y)
        return self
    
    def predict(self, X):
        Xt=X
        for transformer in self.transformers:
            Xt = transformer.transform(Xt)
        return self.estimator.predict(Xt)
        
    def score(self, X, y=None):
        Xt = X
        for transformer in self.transformers:
            Xt = transformer.fit(Xt, y).transform(Xt)
        return self.estimator.score(Xt, y)
 
