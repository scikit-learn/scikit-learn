"""
A Support Vector Machine, this module defines the following classes:

- `CSVCModel`, a model for C-SV classification
- `NuSVCModel`, a model for nu-SV classification
- `EpsilonSVRModel`, a model for epsilon-SV regression
- `NuSVRModel`, a model for nu-SV regression
- `OneClassModel`, a model for distribution estimation (one-class SVM)


How To Use This Module
======================
(See the individual classes, methods, and attributes for details.)

1. Import it: ``import svm`` or ``from svm import ...``.
"""

from classification import *
from regression import *
from oneclass import *
from dataset import *
from kernel import *
from predict import *
