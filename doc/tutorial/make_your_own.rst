
Rule your own scikit-learn estimator
====================================

To implement new algorithms to preprocess, infer, or classify data, one
should implement their own scikit-learn estimator. These estimators need
to follow the scikit-learn API to benefit of all utilities (pipeline,
grid-search, etc.). This tutorial shows how to rule your own
scikit-learn estimator.

Is it necessary to make your own estimator
------------------------------------------

Before to implement your own estimator, be aware that scikit-learn
provides a ``FunctionTransformer`` which will convert a python function
to a scikit-learn transformer. We willustrate the use of this
transformer with a quick example::

    from sklearn.preprocessing import FunctionTransformer
    help(FunctionTransformer)

.. parsed-literal::

    Help on class FunctionTransformer in module sklearn.preprocessing._function_transformer:
    
    class FunctionTransformer(sklearn.base.BaseEstimator, sklearn.base.TransformerMixin)
     |  Constructs a transformer from an arbitrary callable.
     |  
     |  A FunctionTransformer forwards its X (and optionally y) arguments to a
     |  user-defined function or function object and returns the result of this
     |  function. This is useful for stateless transformations such as taking the
     |  log of frequencies, doing custom scaling, etc.
     |  
     |  Note: If a lambda is used as the function, then the resulting
     |  transformer will not be pickleable.
     |  
     |  .. versionadded:: 0.17
     |  
     |  Read more in the :ref:`User Guide <function_transformer>`.
     |  
     |  Parameters
     |  ----------
     |  func : callable, optional default=None
     |      The callable to use for the transformation. This will be passed
     |      the same arguments as transform, with args and kwargs forwarded.
     |      If func is None, then func will be the identity function.
     |  
     |  inverse_func : callable, optional default=None
     |      The callable to use for the inverse transformation. This will be
     |      passed the same arguments as inverse transform, with args and
     |      kwargs forwarded. If inverse_func is None, then inverse_func
     |      will be the identity function.
     |  
     |  validate : bool, optional default=True
     |      Indicate that the input X array should be checked before calling
     |      func. If validate is false, there will be no input validation.
     |      If it is true, then X will be converted to a 2-dimensional NumPy
     |      array or sparse matrix. If this conversion is not possible or X
     |      contains NaN or infinity, an exception is raised.
     |  
     |  accept_sparse : boolean, optional
     |      Indicate that func accepts a sparse matrix as input. If validate is
     |      False, this has no effect. Otherwise, if accept_sparse is false,
     |      sparse matrix inputs will cause an exception to be raised.
     |  
     |  pass_y : bool, optional default=False
     |      Indicate that transform should forward the y argument to the
     |      inner callable.
     |  
     |      .. deprecated::0.19
     |  
     |  check_inverse : bool, default=True
     |     Whether to check that or ``func`` followed by ``inverse_func`` leads to
     |     the original inputs. It can be used for a sanity check, raising a
     |     warning when the condition is not fulfilled.
     |  
     |     .. versionadded:: 0.20
     |  
     |  kw_args : dict, optional
     |      Dictionary of additional keyword arguments to pass to func.
     |  
     |  inv_kw_args : dict, optional
     |      Dictionary of additional keyword arguments to pass to inverse_func.
     |  
     |  Method resolution order:
     |      FunctionTransformer
     |      sklearn.base.BaseEstimator
     |      sklearn.base.TransformerMixin
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self, func=None, inverse_func=None, validate=True, accept_sparse=False, pass_y='deprecated', check_inverse=True, kw_args=None, inv_kw_args=None)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  fit(self, X, y=None)
     |      Fit transformer by checking X.
     |      
     |      If ``validate`` is ``True``, ``X`` will be checked.
     |      
     |      Parameters
     |      ----------
     |      X : array-like, shape (n_samples, n_features)
     |          Input array.
     |      
     |      Returns
     |      -------
     |      self
     |  
     |  inverse_transform(self, X, y='deprecated')
     |      Transform X using the inverse function.
     |      
     |      Parameters
     |      ----------
     |      X : array-like, shape (n_samples, n_features)
     |          Input array.
     |      
     |      y : (ignored)
     |          .. deprecated::0.19
     |      
     |      Returns
     |      -------
     |      X_out : array-like, shape (n_samples, n_features)
     |          Transformed input.
     |  
     |  transform(self, X, y='deprecated')
     |      Transform X using the forward function.
     |      
     |      Parameters
     |      ----------
     |      X : array-like, shape (n_samples, n_features)
     |          Input array.
     |      
     |      y : (ignored)
     |          .. deprecated::0.19
     |      
     |      Returns
     |      -------
     |      X_out : array-like, shape (n_samples, n_features)
     |          Transformed input.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from sklearn.base.BaseEstimator:
     |  
     |  __getstate__(self)
     |  
     |  __repr__(self)
     |      Return repr(self).
     |  
     |  __setstate__(self, state)
     |  
     |  get_params(self, deep=True)
     |      Get parameters for this estimator.
     |      
     |      Parameters
     |      ----------
     |      deep : boolean, optional
     |          If True, will return the parameters for this estimator and
     |          contained subobjects that are estimators.
     |      
     |      Returns
     |      -------
     |      params : mapping of string to any
     |          Parameter names mapped to their values.
     |  
     |  set_params(self, **params)
     |      Set the parameters of this estimator.
     |      
     |      The method works on simple estimators as well as on nested objects
     |      (such as pipelines). The latter have parameters of the form
     |      ``<component>__<parameter>`` so that it's possible to update each
     |      component of a nested object.
     |      
     |      Returns
     |      -------
     |      self
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from sklearn.base.BaseEstimator:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from sklearn.base.TransformerMixin:
     |  
     |  fit_transform(self, X, y=None, **fit_params)
     |      Fit to data, then transform it.
     |      
     |      Fits transformer to X and y with optional parameters fit_params
     |      and returns a transformed version of X.
     |      
     |      Parameters
     |      ----------
     |      X : numpy array of shape [n_samples, n_features]
     |          Training set.
     |      
     |      y : numpy array of shape [n_samples]
     |          Target values.
     |      
     |      Returns
     |      -------
     |      X_new : numpy array of shape [n_samples, n_features_new]
     |          Transformed array.
    


Define a function which will square the input data::

     def square_X(X):
         return X ** 2

Create a transformer using the ``FunctionTransformer``::

    transformer = FunctionTransformer(func=square_X, validate=False)

As any other transformer in scikit-learn, ``transformer`` implements the
``fit_transform`` method::

    import numpy as np
    
    X = np.random.randn(3, 2)
    X
    array([[ 0.37164319,  0.49252007],
           [ 0.38457574, -0.35885232],
           [ 0.66372047,  0.43601824]])

Call ``fit_transform``::

    transformer.fit_transform(X)
    array([[0.13811866, 0.24257602],
           [0.1478985 , 0.12877499],
           [0.44052487, 0.19011191]])



One of the limitations of this transformer is that the ``fit`` method is
actually stateless and one might want to embed some supervision to
transform the data. In this latter case, you need to implement your own
transformer.

Then, make your own estimator
-----------------------------

We will show how to create your own transformer, regressor, and
classifier as well as a quick example to illustrate their integrations
with the scikit-learn utilities.

Scikit-learn base estimator
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The central piece of transformer, regressor, and classifier is
``BaseEstimator``. All estimators in scikit-learn are derived from this
class. In more details, this base class enables to set and get
parameters of the estimator::

    from sklearn.base import BaseEstimator
    help(BaseEstimator)


.. parsed-literal::

    Help on class BaseEstimator in module sklearn.base:
    
    class BaseEstimator(builtins.object)
     |  Base class for all estimators in scikit-learn
     |  
     |  Notes
     |  -----
     |  All estimators should specify all the parameters that can be set
     |  at the class level in their ``__init__`` as explicit keyword
     |  arguments (no ``*args`` or ``**kwargs``).
     |  
     |  Methods defined here:
     |  
     |  __getstate__(self)
     |  
     |  __repr__(self)
     |      Return repr(self).
     |  
     |  __setstate__(self, state)
     |  
     |  get_params(self, deep=True)
     |      Get parameters for this estimator.
     |      
     |      Parameters
     |      ----------
     |      deep : boolean, optional
     |          If True, will return the parameters for this estimator and
     |          contained subobjects that are estimators.
     |      
     |      Returns
     |      -------
     |      params : mapping of string to any
     |          Parameter names mapped to their values.
     |  
     |  set_params(self, **params)
     |      Set the parameters of this estimator.
     |      
     |      The method works on simple estimators as well as on nested objects
     |      (such as pipelines). The latter have parameters of the form
     |      ``<component>__<parameter>`` so that it's possible to update each
     |      component of a nested object.
     |      
     |      Returns
     |      -------
     |      self
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    


Build a scikit-learn transformer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transformers are scikit-lean estimators which implement a ``transform``
method. The use case is the following:

-  at ``fit``, some parameters can be learned from ``X`` and ``y``.
-  at ``transform``, ``X`` will be transformed, using the parameters
   learned during ``fit``.

In addition, scikit-learn provides a
`mixin <https://en.wikipedia.org/wiki/Mixin>`__, i.e.
``TransformerMixin``, which implement the combination of ``fit`` and
``transform`` called ``fit_transform``::

    from sklearn.base import TransformerMixin
    help(TransformerMixin)


.. parsed-literal::

    Help on class TransformerMixin in module sklearn.base:
    
    class TransformerMixin(builtins.object)
     |  Mixin class for all transformers in scikit-learn.
     |  
     |  Methods defined here:
     |  
     |  fit_transform(self, X, y=None, **fit_params)
     |      Fit to data, then transform it.
     |      
     |      Fits transformer to X and y with optional parameters fit_params
     |      and returns a transformed version of X.
     |      
     |      Parameters
     |      ----------
     |      X : numpy array of shape [n_samples, n_features]
     |          Training set.
     |      
     |      y : numpy array of shape [n_samples]
     |          Target values.
     |      
     |      Returns
     |      -------
     |      X_new : numpy array of shape [n_samples, n_features_new]
     |          Transformed array.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    


Therefore, when creating a transformer, you need to create a class which
inherates from both ``BaseEstimator`` and ``TransformerMixin``. The
scikit-learn API imposed ``fit`` to return ``self``. The reason is that
it allows to pipeline ``fit`` and ``transform`` imposed by the
``TransformerMixin``. The ``fit`` method is expected to have ``X`` and
``y`` as inputs. Note that ``transform`` take only ``X`` as input and is
expected to return the transformed version of ``X``::

    class MyOwnTransformer(BaseEstimator, TransformerMixin):
        
        def fit(self, X, y=None):
            return self
        
        def transform(self, X):
            return X

We build a basic example to show that our ``MyOwnTransformer`` is
working within a scikit-learn ``pipeline``::

    from sklearn.datasets import load_iris
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression    

    X, y = load_iris(return_X_y=True)
    pipe = make_pipeline(MyOwnTransformer(), LogisticRegression())
    pipe.fit(X, y)
    Pipeline(memory=None,
         steps=[('myowntransformer', MyOwnTransformer()), ('logisticregression', LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
              intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
              penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
              verbose=0, warm_start=False))])

We can call the ``predict`` method of the pipeline which is equivalent to call
``transform`` of the transformer and ``predict`` of the classifier::

    pipe.predict(X)
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2,
           2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])



Build a scikit-learn regressor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly, regressors are scikit-lean estimators which implement a
``predict`` method. The use case is the following:

-  at ``fit``, some parameters can be learned from ``X`` and ``y``.
-  at ``predict``, predictions will be computed using ``X`` using the
   parameters learned during ``fit``.

In addition, scikit-learn provides a
`mixin <https://en.wikipedia.org/wiki/Mixin>`__, i.e.
``RegressorMixin``, which implement the ``score`` method which compute
the :math:`R^2` score of the predictions::

    from sklearn.base import RegressorMixin
    help(RegressorMixin)


.. parsed-literal::

    Help on class RegressorMixin in module sklearn.base:
    
    class RegressorMixin(builtins.object)
     |  Mixin class for all regression estimators in scikit-learn.
     |  
     |  Methods defined here:
     |  
     |  score(self, X, y, sample_weight=None)
     |      Returns the coefficient of determination R^2 of the prediction.
     |      
     |      The coefficient R^2 is defined as (1 - u/v), where u is the residual
     |      sum of squares ((y_true - y_pred) ** 2).sum() and v is the total
     |      sum of squares ((y_true - y_true.mean()) ** 2).sum().
     |      The best possible score is 1.0 and it can be negative (because the
     |      model can be arbitrarily worse). A constant model that always
     |      predicts the expected value of y, disregarding the input features,
     |      would get a R^2 score of 0.0.
     |      
     |      Parameters
     |      ----------
     |      X : array-like, shape = (n_samples, n_features)
     |          Test samples.
     |      
     |      y : array-like, shape = (n_samples) or (n_samples, n_outputs)
     |          True values for X.
     |      
     |      sample_weight : array-like, shape = [n_samples], optional
     |          Sample weights.
     |      
     |      Returns
     |      -------
     |      score : float
     |          R^2 of self.predict(X) wrt. y.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    


Therefore, we create a regressor, ``MyOwnRegressor`` which inherates
from both ``BaseEstimator`` and ``RegressorMixin``. The method ``fit``
gets ``X`` and ``y`` as input and should return ``self``. It should
implement the ``predict`` function which should output the predictions
of your regressor::

    import numpy as np
    
    
    class MyOwnRegressor(BaseEstimator, RegressorMixin):
        
        def fit(self, X, y):
            return self
        
        def predict(self, X):
            return np.mean(X, axis=1)

We illustrate that this regressor is working within a scikit-learn
pipeline::

    from sklearn.datasets import load_diabetes
    from sklearn.pipeline import make_pipeline
    
    X, y = load_diabetes(return_X_y=True)
    pipe = make_pipeline(MyOwnTransformer(), MyOwnRegressor())
    pipe.fit(X, y)
    Pipeline(memory=None,
         steps=[('myowntransformer', MyOwnTransformer()), ('myownregressor', MyOwnRegressor())])



As we defined the ``predict`` method, we can call it::

    pipe.predict(X)
    array([ 4.95495135e-03, -2.77553225e-02,  3.69509479e-03, -1.33173475e-02,
           -1.07322419e-02, -5.18397864e-02, -2.62834231e-02,  3.86272696e-02,
            7.13601945e-03, -1.30130115e-02, -5.98097614e-02, -4.87315957e-03,
           -1.48189099e-02,  6.16239115e-03, -8.51760901e-03,  3.19136466e-02,
            1.57383528e-02,  2.94368433e-02, -2.53873858e-02, -2.10159471e-02,
           -3.18151715e-02, -2.33380378e-02, -2.29718695e-02,  4.55023173e-02,
           -2.34701457e-02, -1.18886337e-02, -5.39000598e-02, -1.38638642e-02,
           -2.67205855e-02,  1.24185650e-02, -1.50151977e-02, -4.03884229e-02,
            2.06405750e-02, -1.73470330e-02, -4.14075501e-02,  1.85332429e-02,
            5.01600480e-03, -1.65430288e-02,  5.26183537e-02,  3.02657474e-04,
            3.13694636e-02, -5.51624756e-02, -1.94368651e-02, -3.57954358e-03,
            1.69678050e-02, -1.35414972e-03, -3.71529852e-02, -5.78078677e-02,
            1.65562706e-03,  3.13810990e-03, -1.40628247e-02,  2.96302984e-03,
           -3.57248948e-03,  6.30322734e-03, -4.35062345e-03, -2.15455741e-02,
           -1.50952217e-02, -5.27260758e-02, -2.41080182e-03,  2.98740205e-02,
           -4.14279677e-02, -6.47157617e-03, -1.37173389e-02, -2.32307637e-02,
           -1.21493671e-02,  1.35303307e-02, -1.62876677e-05,  3.84719150e-03,
           -6.13526030e-03, -3.02634193e-02, -1.29545502e-02,  2.38035044e-02,
            3.69630716e-02,  1.45607805e-02,  2.45509199e-02, -1.77687143e-03,
           -2.33218087e-02, -4.17374997e-02, -3.54088988e-02, -2.74514575e-02,
            2.75519997e-02,  5.62741521e-03, -1.61994174e-02, -3.40471746e-02,
           -4.45693164e-02, -7.41763370e-03, -3.74660291e-02,  1.46907433e-02,
           -1.46332177e-02, -1.51625140e-02, -2.32399153e-02,  2.07593319e-02,
            6.98580978e-03, -3.50625762e-02, -2.68214953e-02, -3.79436091e-02,
            2.33455251e-02,  1.38298052e-02, -9.92051132e-03, -5.73252056e-03,
            1.34808666e-02,  1.24903916e-02, -5.28170044e-03,  2.32403031e-02,
           -1.92504771e-02, -5.92669442e-03, -5.15939626e-02, -4.84849670e-03,
            6.81343495e-03,  1.96774107e-02, -1.72309569e-02, -5.57438087e-03,
           -1.36757890e-02,  3.07224255e-02,  2.51591276e-02,  2.77188658e-02,
            3.25280933e-02,  3.96861632e-02,  1.44441741e-02, -9.04655596e-03,
           -2.15330990e-02,  1.17377309e-02,  4.08419487e-02,  6.53532198e-02,
           -1.37693091e-02, -1.31079027e-03, -5.14965431e-02,  6.95289037e-03,
           -3.08233883e-02,  1.69130824e-02,  2.67024490e-02, -4.59218716e-02,
            1.27661783e-02, -3.85384349e-02, -1.46468415e-02, -4.89811364e-03,
           -5.21083802e-02,  7.08801004e-03,  2.85229547e-02,  2.81918843e-03,
            5.36616060e-03,  2.46769309e-02,  2.89234323e-02, -1.01910118e-02,
            2.09712031e-02, -5.50697264e-03,  1.37188028e-02,  7.40440158e-03,
            5.64728084e-03, -6.68674148e-03, -3.80270495e-03, -9.45410287e-03,
            3.34837933e-02,  8.38324638e-03,  1.85784144e-02,  4.22735904e-02,
           -2.26737716e-02,  1.15758713e-02, -2.66270154e-02,  1.09370526e-02,
           -3.44857200e-02,  6.33903536e-02, -2.13710392e-02,  2.41752691e-02,
            8.21804910e-03, -3.41164445e-02, -4.97742212e-02,  3.43952096e-02,
            4.80956593e-02,  1.96507708e-02, -1.05580750e-02, -5.00306692e-02,
            4.02125016e-02, -3.59116993e-02, -2.12730935e-02, -1.55994786e-02,
            2.23024617e-02,  1.70510584e-02, -4.40236510e-03, -3.94423448e-03,
            4.98302602e-03, -5.49750114e-03,  3.49330653e-02,  1.77573410e-02,
           -4.39429304e-03, -5.95807614e-03,  1.36771710e-02, -5.70135040e-02,
            1.83508644e-02, -1.48734291e-02, -3.73114087e-03, -2.67290594e-02,
            3.01524510e-03,  2.34874743e-02, -2.63345831e-02,  2.18524182e-02,
           -6.54209011e-03,  2.07556248e-02, -3.45265135e-02,  1.34974387e-02,
            9.76431748e-03, -1.00333887e-02,  4.94707322e-02,  3.13225121e-02,
            2.28773469e-02, -4.20163753e-03,  1.57974348e-02,  1.06402322e-02,
            2.81970885e-02,  1.93096741e-04,  3.46228546e-03, -1.79320588e-03,
            1.92487755e-02, -1.67469212e-02, -2.65273685e-02,  3.00180006e-02,
            4.46264655e-02,  3.36870518e-02, -1.98988927e-02, -3.95716583e-02,
           -2.49333863e-02, -2.41593649e-02, -9.76169063e-03, -4.46992275e-02,
           -4.25931859e-02,  1.11034563e-02, -2.79082135e-02,  1.22030734e-02,
           -2.31601157e-02, -2.95718574e-03,  5.28697349e-02, -5.32651479e-04,
            2.21608010e-02, -1.04281011e-02,  1.85073543e-02,  6.22970486e-03,
            2.77065916e-02, -1.69235022e-02,  7.57776639e-03,  2.83968723e-02,
            2.73032929e-02,  1.47624457e-02, -4.60585464e-02, -1.16170731e-02,
           -1.39082793e-02, -4.40869080e-02,  6.57132367e-03, -3.36041702e-02,
            4.64088989e-02,  2.27492625e-02,  1.76150421e-02,  4.67532432e-02,
            1.47721218e-03,  3.08392277e-02,  4.97098049e-02, -1.34910358e-02,
           -2.16294293e-04, -8.08198743e-03,  1.52753602e-02, -8.53221586e-03,
           -1.14473072e-02,  1.80587401e-02,  3.44046102e-02, -3.52590921e-02,
           -1.09791707e-02,  1.94556764e-03, -1.35800485e-02,  1.74579993e-02,
            4.54087164e-02, -2.95235822e-03,  8.57994170e-03,  2.30560223e-03,
           -3.93032590e-02,  4.35492052e-02, -1.29938524e-03,  6.39434215e-03,
            3.70590644e-02, -2.84159699e-02,  2.22353573e-03, -4.97332287e-03,
            1.60899857e-02, -2.27060906e-02,  2.71266512e-02, -3.46137835e-02,
           -1.29082495e-02,  3.48811534e-03, -2.95363876e-02,  3.09420252e-02,
            1.73627562e-02,  1.52752948e-02,  2.84967963e-02,  2.20616864e-02,
           -3.02468975e-02, -1.40778188e-04, -1.04608252e-02, -2.99920888e-03,
           -2.50406687e-02, -3.33299703e-02, -2.92862193e-02, -4.25134204e-03,
            5.84139106e-03,  1.76901785e-02,  1.57343070e-02,  3.08863866e-02,
            2.95471318e-02, -7.42793313e-03, -8.32377954e-03,  2.62464335e-02,
           -4.42167237e-02,  3.30140231e-02,  6.45401580e-03,  1.91455127e-02,
           -2.75800826e-02,  3.47150182e-02, -1.07313113e-02, -5.31260665e-03,
            1.79540602e-02,  2.57891323e-02,  1.80486489e-02,  6.67745784e-03,
            2.20510771e-02,  4.98713314e-02,  4.65704618e-02,  4.59080535e-02,
            3.09903202e-02,  2.45069709e-02,  1.01006837e-02,  2.33357163e-02,
           -1.23566716e-02, -6.98817291e-03,  1.64175755e-02, -1.35402050e-02,
            2.14887825e-02,  1.19620118e-02, -3.72709334e-02, -1.92104470e-02,
            3.82481084e-02,  1.26971983e-02, -2.58571830e-02,  1.47183136e-02,
           -4.25995302e-03,  9.89301433e-03,  1.61247687e-02,  4.77712592e-04,
           -9.24574723e-03,  2.18848908e-02,  4.88895679e-02, -6.53461476e-03,
           -4.69193475e-03, -3.89582405e-02,  3.86853858e-02, -4.05990949e-02,
           -2.20298556e-02,  7.41182142e-03,  4.63494463e-02, -1.15757874e-02,
           -7.07425355e-03, -1.74009096e-03, -3.85478948e-02,  2.35961578e-02,
            3.15672974e-03, -1.45018831e-02,  1.99298074e-02, -1.04644536e-02,
            1.29652238e-02,  7.37550960e-03,  3.46775206e-02,  3.69625382e-02,
            2.40025073e-02, -8.94537344e-03, -2.56033245e-02,  2.96875003e-02,
            1.52983078e-04, -3.44701298e-02, -4.99483427e-02,  3.30667234e-02,
            3.72462611e-02,  1.14273518e-02,  3.67889210e-03, -3.64328371e-02,
           -1.73823225e-02, -3.55313864e-02,  4.01151519e-03, -8.48642248e-03,
            2.57246654e-02, -9.45999589e-03, -1.22328353e-02, -4.27570610e-02,
            8.73683684e-03, -1.86306736e-02,  3.47139519e-02, -4.72378555e-02,
            6.83410322e-03, -3.46709748e-02,  1.61652765e-02, -1.75695952e-02,
           -3.57307845e-02, -1.12243514e-02, -4.35853283e-04,  2.53651044e-02,
            1.76337340e-04, -2.23408141e-02,  1.90796343e-02,  7.54756037e-03,
           -6.58189685e-03,  3.57774299e-03, -5.40304523e-02,  3.41587890e-03,
            2.31578394e-02,  2.28609434e-02,  1.36224432e-02, -3.04393462e-02,
            3.11416917e-02, -1.13605347e-02,  4.82939948e-02,  2.71697797e-03,
            1.14343316e-02, -1.77919528e-02,  1.60539870e-02, -4.03849177e-02,
           -2.81879503e-04,  3.04866725e-02,  3.01188342e-02,  1.64881087e-02,
            5.85579826e-03, -5.82946559e-02,  3.63936593e-02, -1.79164511e-02,
            5.39543907e-02, -8.73795876e-03, -1.15842817e-02,  6.33430025e-04,
            2.98884968e-02, -7.32320512e-03, -2.53349135e-03, -2.30077446e-02,
           -3.88933705e-02,  1.70663693e-02,  1.22125367e-02, -2.48803176e-03,
           -1.75204817e-04,  1.57815550e-05])



Since we inherite from the ``RegressorMixin``, we can call the ``score``
method which will return the :math:`R^2` score::

    pipe.score(X, y)
    -3.90271854560383



Build a scikit-learn classifier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly to regressors, classifiers implement ``predict``. In addition,
they output the probabilities of the prediction using the
``predict_proba``:

-  at ``fit``, some parameters can be learned from ``X`` and ``y``.
-  at ``predict``, predictions will be computed using ``X`` using the
   parameters learned during ``fit``. It corresponds to the class for
   each sample.
-  ``predict_proba`` will give a 2D matrix where each column corresponds
   to the class and each entry will be the probability to be the
   associated class.

In addition, scikit-learn provides a
`mixin <https://en.wikipedia.org/wiki/Mixin>`__, i.e.
``ClassifierMixin``, which implement the ``score`` method which compute
the accuracy score of the predictions::

    from sklearn.base import ClassifierMixin
    help(ClassifierMixin)


.. parsed-literal::

    Help on class ClassifierMixin in module sklearn.base:
    
    class ClassifierMixin(builtins.object)
     |  Mixin class for all classifiers in scikit-learn.
     |  
     |  Methods defined here:
     |  
     |  score(self, X, y, sample_weight=None)
     |      Returns the mean accuracy on the given test data and labels.
     |      
     |      In multi-label classification, this is the subset accuracy
     |      which is a harsh metric since you require for each sample that
     |      each label set be correctly predicted.
     |      
     |      Parameters
     |      ----------
     |      X : array-like, shape = (n_samples, n_features)
     |          Test samples.
     |      
     |      y : array-like, shape = (n_samples) or (n_samples, n_outputs)
     |          True labels for X.
     |      
     |      sample_weight : array-like, shape = [n_samples], optional
     |          Sample weights.
     |      
     |      Returns
     |      -------
     |      score : float
     |          Mean accuracy of self.predict(X) wrt. y.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    


Therefore, we create a classifier, ``MyOwnClassifier`` which inherates
from both ``BaseEstimator`` and ``ClassifierMixin``. The method ``fit``
gets ``X`` and ``y`` as input and should return ``self``. It should
implement the ``predict`` function which should output the class infered
by the classifier. ``predict_proba`` will output some probabilities
instead::

    import numpy as np
    
    
    class MyOwnClassifier(BaseEstimator, ClassifierMixin):
        
        def fit(self, X, y):
            self.classes_ = np.unique(y)
            return self
        
        def predict(self, X):
            return np.random.randint(0, self.classes_.size, size=X.shape[0])
        
        def predict_proba(self, X):
            pred = np.random.rand(X.shape[0], self.classes_.size)
            return pred / np.sum(pred, axis=1)[:, np.newaxis]

We illustrate that this regressor is working within a scikit-learn
pipeline::

    from sklearn.datasets import load_iris
    from sklearn.pipeline import make_pipeline
    
    X, y = load_iris(return_X_y=True)
    pipe = make_pipeline(MyOwnTransformer(), MyOwnClassifier())
    pipe.fit(X, y)
    Pipeline(memory=None,
         steps=[('myowntransformer', MyOwnTransformer()), ('myownclassifier', MyOwnClassifier())])



Then, you can call ``predict`` and ``predict_proba``::

    pipe.predict(X)
    array([1, 0, 0, 2, 1, 2, 1, 0, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 0,
           1, 1, 2, 2, 0, 0, 0, 1, 1, 2, 1, 1, 2, 2, 0, 2, 2, 2, 1, 0, 2, 2,
           0, 0, 0, 2, 2, 1, 2, 2, 0, 0, 1, 0, 1, 2, 2, 1, 1, 1, 2, 2, 0, 1,
           1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 2, 1, 0, 0, 2, 1, 1, 1, 2, 1,
           2, 2, 0, 1, 2, 2, 1, 0, 2, 2, 1, 0, 2, 1, 2, 0, 2, 1, 1, 1, 1, 2,
           1, 2, 0, 0, 0, 1, 0, 2, 1, 0, 0, 2, 1, 1, 1, 0, 1, 2, 1, 2, 1, 0,
           2, 2, 2, 2, 0, 2, 2, 1, 0, 1, 0, 2, 2, 0, 2, 0, 1, 0])

    pipe.predict_proba(X)
    array([[1.48898954e-01, 2.74767300e-01, 5.76333746e-01],
           [1.73879547e-01, 5.42471890e-01, 2.83648563e-01],
           [3.40340578e-01, 2.04503856e-01, 4.55155566e-01],
           [2.73222522e-01, 3.40804800e-01, 3.85972678e-01],
           [4.55841326e-01, 1.26782683e-01, 4.17375991e-01],
           [4.07918869e-01, 7.10331214e-02, 5.21048010e-01],
           [3.89824021e-01, 9.53583469e-02, 5.14817632e-01],
           [3.34717625e-02, 6.06390796e-01, 3.60137441e-01],
           [5.38870842e-01, 3.30227721e-01, 1.30901437e-01],
           [4.51187958e-01, 2.41165811e-01, 3.07646231e-01],
           [5.02403349e-01, 2.17128334e-01, 2.80468316e-01],
           [4.17673966e-01, 1.23342082e-01, 4.58983953e-01],
           [1.65928476e-01, 4.00305006e-01, 4.33766518e-01],
           [1.41536182e-01, 3.71696454e-01, 4.86767364e-01],
           [1.62912805e-01, 4.42621802e-01, 3.94465393e-01],
           [7.80468447e-01, 6.64741320e-02, 1.53057421e-01],
           [1.48624816e-01, 1.24632447e-01, 7.26742738e-01],
           [5.23250715e-01, 8.01713046e-02, 3.96577981e-01],
           [2.53780185e-01, 5.65176803e-01, 1.81043012e-01],
           [2.61667911e-01, 4.27884633e-01, 3.10447456e-01],
           [1.21001070e-01, 7.83007673e-01, 9.59912567e-02],
           [1.70184523e-01, 4.09286445e-01, 4.20529032e-01],
           [1.13115488e-01, 7.26555438e-01, 1.60329074e-01],
           [2.16316478e-01, 3.79334940e-01, 4.04348582e-01],
           [3.76914968e-01, 5.16871717e-01, 1.06213315e-01],
           [4.53046131e-01, 3.20084822e-01, 2.26869047e-01],
           [1.43569949e-01, 8.38202897e-02, 7.72609761e-01],
           [3.70413698e-01, 3.35522492e-01, 2.94063809e-01],
           [3.38715613e-01, 1.70426456e-01, 4.90857931e-01],
           [4.52648140e-01, 4.88974767e-01, 5.83770923e-02],
           [1.08176110e-01, 6.99976862e-01, 1.91847028e-01],
           [5.57591627e-01, 2.70846629e-01, 1.71561744e-01],
           [2.55013573e-01, 3.12993395e-01, 4.31993032e-01],
           [4.03957154e-01, 3.75145549e-01, 2.20897297e-01],
           [2.35594332e-01, 2.92020985e-01, 4.72384683e-01],
           [2.84544466e-01, 2.97413490e-01, 4.18042044e-01],
           [2.79385976e-01, 3.30411606e-01, 3.90202419e-01],
           [1.32840434e-02, 7.80741838e-01, 2.05974118e-01],
           [1.45895347e-01, 4.99359142e-01, 3.54745511e-01],
           [2.89221905e-01, 5.23204803e-01, 1.87573292e-01],
           [3.56456301e-01, 3.95208886e-01, 2.48334813e-01],
           [3.72995472e-01, 5.17896994e-01, 1.09107534e-01],
           [1.46861418e-01, 2.24022597e-01, 6.29115985e-01],
           [6.65255942e-01, 3.07417284e-02, 3.04002329e-01],
           [3.61284287e-01, 1.56733057e-01, 4.81982656e-01],
           [6.48294156e-02, 8.13166650e-01, 1.22003934e-01],
           [3.40012707e-01, 9.84115412e-02, 5.61575752e-01],
           [6.39786237e-02, 5.19155271e-01, 4.16866106e-01],
           [4.32514738e-01, 1.53632946e-02, 5.52121968e-01],
           [7.35941156e-01, 9.51042521e-03, 2.54548418e-01],
           [5.37303626e-01, 4.09293716e-01, 5.34026581e-02],
           [4.03674350e-01, 5.52535758e-01, 4.37898923e-02],
           [3.51789299e-01, 5.90980201e-01, 5.72304999e-02],
           [4.06431970e-01, 1.36587365e-01, 4.56980665e-01],
           [8.26829364e-01, 1.43241972e-01, 2.99286641e-02],
           [3.34666284e-01, 3.31664400e-01, 3.33669316e-01],
           [3.19039435e-01, 3.44069000e-01, 3.36891566e-01],
           [6.44112115e-01, 7.38850136e-02, 2.82002871e-01],
           [2.06994523e-01, 3.72068916e-01, 4.20936561e-01],
           [6.69881801e-01, 7.43099910e-02, 2.55808208e-01],
           [4.62031286e-01, 3.21817873e-02, 5.05786926e-01],
           [5.65133819e-01, 1.64425383e-01, 2.70440798e-01],
           [3.52255585e-01, 1.74876113e-01, 4.72868302e-01],
           [4.92328463e-01, 4.96106788e-01, 1.15647493e-02],
           [1.37283500e-01, 6.59843624e-01, 2.02872876e-01],
           [3.50740744e-01, 1.17685058e-02, 6.37490751e-01],
           [4.90337630e-01, 3.23346873e-01, 1.86315498e-01],
           [7.98536545e-02, 3.83129645e-01, 5.37016700e-01],
           [1.93241203e-02, 2.72971335e-01, 7.07704545e-01],
           [2.00924763e-01, 7.94905024e-02, 7.19584735e-01],
           [2.96847602e-01, 3.72535888e-01, 3.30616510e-01],
           [3.25660366e-01, 4.58555752e-01, 2.15783882e-01],
           [6.06678796e-01, 2.63787113e-01, 1.29534091e-01],
           [1.21525442e-01, 4.36753098e-01, 4.41721460e-01],
           [4.12912148e-01, 2.71237916e-01, 3.15849935e-01],
           [3.72959038e-01, 1.45348986e-01, 4.81691976e-01],
           [7.86603572e-02, 5.58061881e-01, 3.63277761e-01],
           [4.05114827e-01, 2.80981381e-01, 3.13903792e-01],
           [1.19049841e-01, 4.29268271e-01, 4.51681888e-01],
           [2.56598240e-01, 3.03001728e-01, 4.40400032e-01],
           [3.51620949e-01, 3.70244784e-01, 2.78134267e-01],
           [3.33911157e-01, 5.37566937e-01, 1.28521905e-01],
           [7.35379599e-04, 6.00140416e-01, 3.99124205e-01],
           [3.07387456e-01, 1.79414474e-01, 5.13198070e-01],
           [5.73680705e-02, 5.27766539e-01, 4.14865390e-01],
           [4.42004157e-02, 4.00301891e-01, 5.55497693e-01],
           [2.38346308e-01, 5.96633930e-02, 7.01990299e-01],
           [8.41616001e-02, 6.58366806e-01, 2.57471594e-01],
           [4.87858190e-01, 4.67632187e-02, 4.65378592e-01],
           [2.88355346e-01, 2.10520131e-01, 5.01124523e-01],
           [6.68212253e-02, 7.21724450e-02, 8.61006330e-01],
           [3.73506884e-01, 1.40626429e-01, 4.85866687e-01],
           [7.84455432e-01, 2.08725085e-01, 6.81948302e-03],
           [3.51774080e-01, 3.51818476e-01, 2.96407444e-01],
           [1.89650591e-01, 5.25535923e-01, 2.84813486e-01],
           [4.45909016e-01, 3.11768689e-01, 2.42322295e-01],
           [5.57534911e-01, 7.29017551e-02, 3.69563334e-01],
           [4.46838200e-01, 3.73624622e-01, 1.79537179e-01],
           [4.77553384e-01, 1.90452065e-01, 3.31994551e-01],
           [4.04189524e-01, 1.53355875e-01, 4.42454601e-01],
           [2.88342819e-02, 7.98497501e-01, 1.72668217e-01],
           [3.80362363e-01, 4.71491910e-01, 1.48145728e-01],
           [3.14584055e-01, 2.50005203e-01, 4.35410742e-01],
           [4.72201899e-02, 4.15702028e-01, 5.37077782e-01],
           [3.18343962e-01, 1.98397909e-01, 4.83258129e-01],
           [6.75514300e-03, 5.21752707e-01, 4.71492150e-01],
           [5.97604127e-01, 3.10738334e-01, 9.16575395e-02],
           [3.43189922e-01, 3.61407776e-01, 2.95402302e-01],
           [7.51405160e-01, 6.05593770e-02, 1.88035463e-01],
           [3.69234805e-01, 5.15296565e-02, 5.79235538e-01],
           [3.57791649e-01, 5.10333011e-01, 1.31875340e-01],
           [1.53787744e-01, 4.02174912e-01, 4.44037344e-01],
           [9.53375282e-02, 5.39248412e-01, 3.65414060e-01],
           [5.27205914e-01, 1.34321961e-01, 3.38472125e-01],
           [4.94381695e-01, 4.86386452e-02, 4.56979660e-01],
           [3.47911885e-01, 5.97250963e-01, 5.48371518e-02],
           [5.16615519e-01, 2.70684769e-01, 2.12699712e-01],
           [3.21296395e-01, 4.95534985e-01, 1.83168620e-01],
           [1.60978331e-01, 4.09468666e-01, 4.29553003e-01],
           [2.74577501e-01, 2.36423295e-01, 4.88999204e-01],
           [6.20578839e-01, 3.52039413e-01, 2.73817480e-02],
           [2.74417432e-01, 3.62480342e-01, 3.63102226e-01],
           [3.97506959e-01, 4.72721462e-01, 1.29771580e-01],
           [3.30688511e-01, 4.08164424e-01, 2.61147065e-01],
           [2.72902049e-01, 4.08329729e-01, 3.18768222e-01],
           [5.98235030e-01, 7.29950347e-02, 3.28769935e-01],
           [1.15081073e-01, 9.64678612e-02, 7.88451065e-01],
           [2.25129254e-01, 3.64459123e-01, 4.10411622e-01],
           [5.26159390e-01, 4.51765893e-01, 2.20747172e-02],
           [3.28611157e-01, 5.44294205e-01, 1.27094638e-01],
           [4.10476284e-01, 4.07370671e-01, 1.82153045e-01],
           [3.35105803e-01, 2.09223129e-01, 4.55671068e-01],
           [4.45056257e-01, 3.45589636e-01, 2.09354107e-01],
           [1.21535967e-01, 1.62457311e-01, 7.16006722e-01],
           [4.40941371e-01, 3.08953704e-01, 2.50104925e-01],
           [3.61750556e-01, 2.37802560e-01, 4.00446884e-01],
           [2.71954263e-01, 3.22637237e-01, 4.05408500e-01],
           [8.31254904e-02, 3.46861591e-01, 5.70012919e-01],
           [3.33654044e-01, 2.56878434e-01, 4.09467522e-01],
           [1.86404801e-01, 7.00382576e-01, 1.13212622e-01],
           [3.32215493e-01, 5.82748755e-01, 8.50357520e-02],
           [3.42720814e-01, 4.99159225e-01, 1.58119961e-01],
           [2.18128822e-01, 4.88987489e-01, 2.92883689e-01],
           [5.33163504e-01, 2.23235840e-01, 2.43600656e-01],
           [4.11810931e-01, 8.07591603e-02, 5.07429909e-01],
           [3.28487999e-01, 1.09963062e-02, 6.60515695e-01],
           [3.78658861e-01, 1.63968828e-01, 4.57372310e-01],
           [4.59974973e-01, 1.52603969e-01, 3.87421057e-01],
           [7.88297624e-02, 2.85798761e-01, 6.35371476e-01],
           [3.92767253e-01, 1.88362006e-01, 4.18870742e-01]])



Since our classifier inherites from ``ClassifierMixin``, we can compute
the accuracy by calling the ``score`` method::

    pipe.score(X, y)
    0.36

