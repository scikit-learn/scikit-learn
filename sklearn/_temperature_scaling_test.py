'''
This file is created to test if the custom 'TemperatureScaling' class runs properly,
and serves as proof of work for the changes made to the scikit-learn repository.

The file also includes examples related to developing a temperature scaling method
for probability calibration in multi-class classification.


References:
-----------
  .. [1] https://github.com/scikit-learn/scikit-learn/issues/28574. Original issue
         on Github.
'''

from sklearn.calibration_temperature import CalibratedClassifierCV_test
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier

# Load the Iris dataset
X, y = datasets.load_iris(return_X_y=True)
X_train, X_calib, y_train, y_calib = train_test_split(X, y)

# Load the following classifiers for testing
SV_classifier: SVC = SVC(probability=False)
Logistic_classifier: LogisticRegression = LogisticRegression()
Tree_classifier: DecisionTreeClassifier = DecisionTreeClassifier()


# Initiate the calibrators for the classifiers
SVC_scaled: CalibratedClassifierCV_test = CalibratedClassifierCV_test(SV_classifier,
                                                                      cv=3,
                                                                      method='temperature'
                                                                      )
Logistic_scaled: CalibratedClassifierCV_test = CalibratedClassifierCV_test(Logistic_classifier,
                                                                           cv=7,
                                                                           method='temperature'
                                                                           )
Tree_scaled: CalibratedClassifierCV_test = CalibratedClassifierCV_test(Tree_classifier,
                                                                       cv=3,
                                                                       method='temperature'
                                                                       )


# Fit all classifier-calibrator pairs
SVC_scaled.fit(X_train,y_train)
Logistic_scaled.fit(X_train,y_train)
Tree_scaled.fit(X_train,y_train)
# compare_classifier.fit(X_train, y_train)

print("Optimal Temperatures For Each Classifiers")
print(f"{SVC_scaled.calibrated_classifiers_[0].calibrators[0].T_=}")
print(f"{Logistic_scaled.calibrated_classifiers_[0].calibrators[0].T_=}")
print(f"{Tree_scaled.calibrated_classifiers_[0].calibrators[0].T_=}")

print("Printing calibrated probabilities...")
print(f"{SVC_scaled.predict_proba(X_calib) = }")
print(f"{Logistic_scaled.predict_proba(X_calib) = }")
print(f"{Tree_scaled.predict_proba(X_calib)=}")
# print(f"{compare_classifier.predict_proba(X_calib)=}")
print(f"{y_calib=}")
