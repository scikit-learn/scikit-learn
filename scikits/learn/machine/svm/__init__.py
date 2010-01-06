"""
A Support Vector Machine, this module defines the following classes:

- `CClassificationModel`, a model for C-SV classification
- `NuClassificationModel`, a model for nu-SV classification
- `EpsilonRegressionModel`, a model for epsilon-SV regression
- `NuRegressionModel`, a model for nu-SV regression
- `OneClassModel`, a model for distribution estimation
  (one-class SVM)

Kernel classes:

- `Linear`, a linear kernel
- `Polynomial`, a polynomial kernel
- `RBF`, a radial basis function kernel
- `Sigmoid`, a sigmoid kernel
- `Custom`, a kernel that wraps any callable

Dataset classes:

- `ClassificationDataSet`, a dataset for training classification
  models
- `RegressionDataSet`, a dataset for training regression models
- `OneClassDataSet`, a dataset for training distribution
  estimation (one-class SVM) models
- `TestDataSet`, a dataset for testing with any model

Data type classes:

- `svm_node_dtype`, the libsvm data type for its arrays

How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import svm`` or ``from svm import ...``.

2. Create a training dataset for your problem::

       traindata = ClassificationDataSet(labels, x)
       traindata = RegressionDataSet(y, x)
       traindata = OneClassDataSet(x)

   where x is sequence of NumPy arrays containing scalars or
   svm_node_dtype entries.

3. Create a test dataset::

       testdata = TestDataSet(u)

4. Create a model and fit it to the training data::

       model = CClassificationModel(kernel)
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
