from .base import BaseEstimator

class Pipeline(BaseEstimator):
    """
    Univariate ranking Pipeline.
    Sequentialy apply a list of methods
    
    Example:
    import datamind.ml.dimred.univRankingGLM as drglm
    d=data.twoclassses()
    y=d[:,["class"]]
    X=d[:,1:]
    # permute features
    X=X[:,N.random.permutation(X.shape[1])]
    print X
    methods=[drglm.UnivOneWayAnovaFstat(pval=.3),
             drglm.UnivOneWayAnovaFstat(z=1.5),
             drglm.UnivOneWayAnovaFstat(dim=2)]
    pipeline=UnivRankingPipeline(methods)
    pipeline.fit(X,y)
    print pipeline.reduce(X)
    print pipeline.getSelectedFeatures()
    """
    def __init__(self,transformers,estimator):
        """
        methods: list of UnivRanking objects,
        ie.: with fit/reduce/getSelectedFeatures methods
        """
        self.transformers=transformers
        self.estimator=estimator

    def fit(self, X,y):
        Xt=X
        for transformer in self.transformers:
            Xt=transformer.fit(Xt,y).transform(Xt)
        estimator.fit(Xt,y)
        return self
    
    def predict(X)
        Xt=X
        for transformer in self.transformers:
            Xt=transformer.transform(Xt)
        return estimator.predict(Xt,y)
        
