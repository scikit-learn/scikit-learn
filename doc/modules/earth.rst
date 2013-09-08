.. _earth:

========================================
Multivariate Adaptive Regression Splines
========================================

.. currentmodule:: sklearn.earth

Multivariate adaptive regression splines, implemented by the :class:`Earth` class, is a method for supervised 
learning that is most commonly used for feature extraction and selection.  ``Earth`` models can be thought of as linear models in a higher dimensional 
basis space.  ``Earth`` automatically searches for interactions and non-linear relationships.  Each term in an ``Earth`` model is a 
product of so called "hinge functions".  A hinge function is a function that's equal to its argument where that argument 
is greater than zero and is zero everywhere else.

.. math::
	\text{h}\left(x-t\right)=\left[x-t\right]_{+}=\begin{cases}
	x-t, & x>t\\
	0, & x\leq t
	\end{cases} 

.. image:: ../images/hinge.png

An ``Earth`` model is a linear combination of basis functions, each of which is a product of one 
or more of the following:

	1. A constant
	2. Linear functions of input variables
	3. Hinge functions of input variables  

For example, a simple piecewise linear function in one variable can be expressed 
as a linear combination of two hinge functions and a constant (see below).  During fitting, the ``Earth`` class 
automatically determines which variables and basis functions to use.  
The algorithm has two stages.  First, the 
forward pass searches for terms that locally minimize squared error loss on the training set.  Next, a pruning pass selects a subset of those 
terms that produces a locally minimal generalized cross-validation (GCV) score.  The GCV 
score is not actually based on cross-validation, but rather is meant to approximate a true
cross-validation score by penalizing model complexity.  The final result is a set of basis functions
that is nonlinear in the original feature space, may include interactions, and is likely to 
generalize well.


.. math::
	y=1-2\text{h}\left(1-x\right)+\frac{1}{2}\text{h}\left(x-1\right)


.. image:: ../images/piecewise_linear.png


A Simple Earth Example
----------------------


::

	import numpy
	from pyearth import Earth
	from matplotlib import pyplot
    
	#Create some fake data
	numpy.random.seed(0)
	m = 1000
	n = 10
	X = 80*numpy.random.uniform(size=(m,n)) - 40
	y = numpy.abs(X[:,6] - 4.0) + 1*numpy.random.normal(size=m)
    
	#Fit an Earth model
	model = Earth()
	model.fit(X,y)
    
	#Print the model
	print model.trace()
	print model.summary()
    
	#Plot the model
	y_hat = model.predict(X)
	pyplot.figure()
	pyplot.plot(X[:,6],y,'r.')
	pyplot.plot(X[:,6],y_hat,'b.')
	pyplot.xlabel('x_6')
	pyplot.ylabel('y')
	pyplot.title('Simple Earth Example')
	pyplot.show()

.. image:: ../images/simple_earth_example.png


.. topic:: Bibliography:

	1. Friedman, J. (1991). Multivariate adaptive regression splines. The annals of statistics, 
	   19(1), 1–67. http://www.jstor.org/stable/10.2307/2241837
	2. Stephen Milborrow. Derived from mda:mars by Trevor Hastie and Rob Tibshirani.
	   (2012). earth: Multivariate Adaptive Regression Spline Models. R package
	   version 3.2-3.
	3. Friedman, J. (1993). Fast MARS. Stanford University Department of Statistics, Technical Report No 110. 
	   http://statistics.stanford.edu/~ckirby/techreports/LCS/LCS%20110.pdf
	4. Friedman, J. (1991). Estimating functions of mixed ordinal and categorical variables using adaptive splines.
	   Stanford University Department of Statistics, Technical Report No 108. 
	   http://statistics.stanford.edu/~ckirby/techreports/LCS/LCS%20108.pdf
	5. Stewart, G.W. Matrix Algorithms, Volume 1: Basic Decompositions. (1998). Society for Industrial and Applied 
	   Mathematics.
	6. Bjorck, A. Numerical Methods for Least Squares Problems. (1996). Society for Industrial and Applied 
	   Mathematics.
	7. Hastie, T., Tibshirani, R., & Friedman, J. The Elements of Statistical Learning (2nd Edition). (2009).  
	   Springer Series in Statistics
	8. Golub, G., & Van Loan, C. Matrix Computations (3rd Edition). (1996). Johns Hopkins University Press.
	   

	References 7, 2, 1, 3, and 4 contain discussions likely to be useful to users.  References 1, 2, 6, 5, 
	8, 3, and 4 are useful in understanding the implementation.
