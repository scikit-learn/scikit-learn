"""
A Support Vector Machine, this module defines the following classes:

- `LibSvmCClassificationModel`, a model for C-SV classification
- `LibSvmNuClassificationModel`, a model for nu-SV classification
- `LibSvmEpsilonRegressionModel`, a model for epsilon-SV regression
- `LibSvmNuRegressionModel`, a model for nu-SV regression
- `LibSvmOneClassModel`, a model for distribution estimation
  (one-class SVM)

Kernel classes:

- `LinearKernel`, a linear kernel
- `PolynomialKernel`, a polynomial kernel
- `RBFKernel`, a radial basis function kernel
- `SigmoidKernel`, a sigmoid kernel
- `CustomKernel`, a kernel that wraps any callable

Dataset classes:

- `LibSvmClassificationDataSet`, a dataset for training classification
  models
- `LibSvmRegressionDataSet`, a dataset for training regression models
- `LibSvmOneClassDataSet`, a dataset for training distribution
  estimation (one-class SVM) models
- `LibSvmTestDataSet`, a dataset for testing with any model

Data type classes:

- `svm_node_dtype`, the libsvm data type for its arrays

How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import svm`` or ``from svm import ...``.

2. Create a training dataset for your problem::

       traindata = LibSvmClassificationDataSet(labels, x)
       traindata = LibSvmRegressionDataSet(y, x)
       traindata = LibSvmOneClassDataSet(x)

   where x is sequence of NumPy arrays containing scalars or
   svm_node_dtype entries.

3. Create a test dataset::

       testdata = LibSvmTestDataSet(u)

4. Create a model and fit it to the training data::

       model = LibSvmCClassificationModel(kernel)
       results = model.fit(traindata)

5. Use the results to make predictions with the test data::

       p = results.predict(testdata)
       v = results.predict_values(testdata)
"""

from classification import *
from regression import *
from oneclass import *
from dataset import *
from kernel import *
from predict import *
