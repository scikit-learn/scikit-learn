import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.base import RegressorMixin
from sklearn.base import ClassifierMixin

   
class BaseSupervisedPCA(object):
    """
    Supervised PCA algorithm proposed by Bair et al. (2006).
    
    
    Parameters
    ----------
    
    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).
        
    model : The supervised learning model that will be used to conduct supervised PCA.
    
    Attributes
    ----------
        
    
    References
    ----------
    Bair, Eric, et al. "Prediction by supervised principal components." Journal of the American Statistical Association 101.473 (2006).

    
    """
    
    def __init__(self, fit_intercept=True, model=None):
        self.fit_intercept = fit_intercept
        self._model=model
        
        self._pca=None
        self._leavouts=None
        
    def fit(self,X,y,threshold=0,n_components=2):
        """
        Fit the supervised PCA model
        .
        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples, n_targets]
            Target values
        threshold : the threshold for the coefficient below which it is discarded.
        n_components : the number of components to keep, after running PCA
        
        Returns
        -------
        self : returns an instance of self.
        """
        
        #these are the columns that will be removed
        self._leavouts=[]        
                
        dummy_X=X[:,np.newaxis]
        
        #test the absolute value of the coefficient for each variable. If it
        #is below a the threshold, then add it to the leaveouts        
        for i in range(0,dummy_X.shape[2]):
            current_X=dummy_X[:,:,i]
            self._model.fit(current_X, y)
            #the all([]) syntax is there in order to support both linear and logistic
            #regression. Logistic regression coefficients for multiclass problems
            #come in multi-dimensional arrays.
            if(all([abs(self._model.coef_[0])<threshold])):
                self._leavouts.append(i)
                
        #delete the variables that were below the threshold
        dummy_X=np.delete(dummy_X,self._leavouts,2)
        
        #conduct PCA for the designated number of components
        self._pca = PCA(n_components=n_components)
        dummy_X=self._pca.fit_transform(dummy_X[:,0,:])
        
        self._model=self._model.fit(dummy_X,y)
        
        return self
        
    def predict(self,X):
        """Predict using the supervised PCA model
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.
        Returns
        -------
        C : array, shape = (n_samples,)
            Returns predicted values.        
        """
        #remove the leavouts, transform the data and fit the regression model
        transformed_X=self.get_transformed_data(X)
        return self._model.predict(transformed_X)
    
    def get_transformed_data(self,X):
        """Calculates the components on a new matrix.
        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            
        Returns
        -------
        transformed_X: Returns a transformed numpy array or sparse matrix. The
        leavouts have been removed and the remaining variables are transformed into
        components using the weights of the PCA object.
        
        Notes
        -------
        The algorithm should have first been executed on a dataset.
        
        """
        transformed_X=np.delete(X,self._leavouts,1)
        transformed_X=self._pca.transform(transformed_X)
        return transformed_X
    
    
    #I am implementing a function here to get the components in order to avoid
    #the user having to access the pca object. Another option would be to 
    #copy the components from the pca to a variable located at 'self'. However,
    #this might be too redundant.
    def get_components(self):
        """Returns the components formerly calculated on a training dataset.
            
        Returns
        -------
        components: A numpy matrix with the loadings of the PCA components.
        
        Notes
        -------
        The algorithm should have first been executed on a dataset.
        
        """
        return self._pca.components_
    
    #same principle as in the get_components function
    def get_coefs(self):
        return self._model.coefs_
    
        
        
class SupervisedPCARegressor(BaseSupervisedPCA,RegressorMixin):
    """
    Implementation of supervisedPCA for regression. The underlying model
    is a linear regression model.
    
    Parameters
    ----------
    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.
    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.
    n_jobs : int, optional, default 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. This will only provide speedup for
        n_targets > 1 and sufficient large problems.
    Attributes
    ----------
    coef_ : array, shape (n_features, ) or (n_targets, n_features)
        Estimated coefficients for the linear regression problem.
        If multiple targets are passed during the fit (y 2D), this
        is a 2D array of shape (n_targets, n_features), while if only
        one target is passed, this is a 1D array of length n_features.
    intercept_ : array
        Independent term in the linear model.
    
    """
    def __init__(self,fit_intercept=True, normalize=False, copy_X=True,n_jobs=1):
        model=LinearRegression(copy_X=copy_X,normalize=normalize,n_jobs=n_jobs)  
        super(SupervisedPCARegressor,self).__init__(fit_intercept=fit_intercept,model=model)


class SupervisedPCAClassifier(BaseSupervisedPCA,ClassifierMixin):
    """Implementation of supervisedPCA for classification. The underlying model
    is a logistic regression model.

    Parameters
    ----------
    penalty : str, 'l1' or 'l2'
        Used to specify the norm used in the penalization. The newton-cg and
        lbfgs solvers support only l2 penalties.
    dual : bool
        Dual or primal formulation. Dual formulation is only implemented for
        l2 penalty with liblinear solver. Prefer dual=False when
        n_samples > n_features.
    C : float, optional (default=1.0)
        Inverse of regularization strength; must be a positive float.
        Like in support vector machines, smaller values specify stronger
        regularization.
    fit_intercept : bool, default: True
        Specifies if a constant (a.k.a. bias or intercept) should be
        added the decision function.
    intercept_scaling : float, default: 1
        Useful only if solver is liblinear.
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased.
    class_weight : dict or 'balanced', optional
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.
        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``
    max_iter : int
        Useful only for the newton-cg and lbfgs solvers. Maximum number of
        iterations taken for the solvers to converge.
    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.
    solver : {'newton-cg', 'lbfgs', 'liblinear'}
        Algorithm to use in the optimization problem.
    tol : float, optional
        Tolerance for stopping criteria.
    multi_class : str, {'ovr', 'multinomial'}
        Multiclass option can be either 'ovr' or 'multinomial'. If the option
        chosen is 'ovr', then a binary problem is fit for each label. Else
        the loss minimised is the multinomial loss fit across
        the entire probability distribution. Works only for the 'lbfgs'
        solver.
    verbose : int
        For the liblinear and lbfgs solvers set verbose to any positive
        number for verbosity.

    """
    def __init__(self,fit_intercept=True, normalize=False, copy_X=True,penalty='l2', dual=False, tol=1e-4, C=1.0,
                 intercept_scaling=1, class_weight=None,
                 random_state=None, solver='liblinear', max_iter=100,
                 multi_class='ovr', verbose=0):
        model=LogisticRegression()  
        super(SupervisedPCAClassifier,self).__init__(fit_intercept=fit_intercept,model=model)