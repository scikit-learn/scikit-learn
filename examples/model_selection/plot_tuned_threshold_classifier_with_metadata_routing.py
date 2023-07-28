# %%
from sklearn.datasets import fetch_openml

credit_card = fetch_openml(data_id=1597, as_frame=True, parser="pandas")
credit_card.frame

# %%
data = credit_card.frame.drop(columns=["Class", "Amount"])
target = credit_card.frame["Class"].astype(int).to_numpy()
amount = credit_card.frame["Amount"].to_numpy()

# %%
from collections import Counter

Counter(target)

# %%
import matplotlib.pyplot as plt

fraud = target == 1
amount_fraud = amount[fraud]
_, ax = plt.subplots()
ax.hist(amount_fraud, bins=100)
ax.set_title("Amount of fraud transaction")
_ = ax.set_xlabel("Amount ($)")

# %%


def business_metric(y_true, y_pred, amount):
    mask_true_positive = (y_true == 1) & (y_pred == 1)
    mask_true_negative = (y_true == 0) & (y_pred == 0)
    mask_false_positive = (y_true == 0) & (y_pred == 1)
    mask_false_negative = (y_true == 1) & (y_pred == 0)
    fraudulent_refuse = (mask_true_positive.sum() * 50) + amount[
        mask_true_positive
    ].sum()
    fraudulent_accept = -amount[mask_false_negative].sum()
    legitimate_refuse = mask_false_positive.sum() * -5
    legitimate_accept = (amount[mask_true_negative] * 0.02).sum()
    return fraudulent_refuse + fraudulent_accept + legitimate_refuse + legitimate_accept


# %%
import numpy as np

benefit_cost = business_metric(target, np.zeros_like(target), amount)
print(f"Benefit/cost of not detecting fraud: ${benefit_cost:,.2f}")

# %%
benefit_cost = business_metric(target, np.ones_like(target), amount)
print(f"Benefit/cost of tagging everything as fraud: ${benefit_cost:,.2f}")

# %%
from sklearn.model_selection import train_test_split

data_train, data_test, target_train, target_test, amount_train, amount_test = (
    train_test_split(
        data, target, amount, stratify=target, test_size=0.5, random_state=42
    )
)

# %%
import sklearn
from sklearn.metrics import make_scorer

sklearn.set_config(enable_metadata_routing=True)

business_scorer = make_scorer(business_metric).set_score_request(amount=True)

# %%
from sklearn.dummy import DummyClassifier

easy_going_classifier = DummyClassifier(strategy="constant", constant=0)
easy_going_classifier.fit(data_train, target_train)

# %%
benefit_cost = business_scorer(
    easy_going_classifier, data_test, target_test, amount=amount_test
)
print(f"Benefit/cost of our easy-going classifier: ${benefit_cost:,.2f}")

# %%
intolerant_classifier = DummyClassifier(strategy="constant", constant=1)
intolerant_classifier.fit(data_train, target_train)

# %%
benefit_cost = business_scorer(
    intolerant_classifier, data_test, target_test, amount=amount_test
)
print(f"Benefit/cost of our intolerant classifier: ${benefit_cost:,.2f}")

# %%
from sklearn.linear_model import LogisticRegression

logistic_regression = LogisticRegression(max_iter=1_000).fit(data_train, target_train)
benefit_cost = business_scorer(
    logistic_regression, data_test, target_test, amount=amount_test
)
print(f"Benefit/cost of our logistic regression: ${benefit_cost:,.2f}")

# %%
from sklearn.metrics import get_scorer

balanced_accuracy_scorer = get_scorer("balanced_accuracy")
balanced_accuracy = balanced_accuracy_scorer(
    logistic_regression, data_test, target_test
)
print(f"Balanced accuracy of our logistic regression: {balanced_accuracy:.2f}")

# %%
from sklearn.model_selection import TunedThresholdClassifier

tuned_model = TunedThresholdClassifier(
    estimator=logistic_regression,
    objective_metric=business_scorer,
    n_thresholds=1_000,
    n_jobs=-1,
)
tuned_model.fit(data_train, target_train, amount=amount_train)

# %%
benefit_cost = business_scorer(tuned_model, data_test, target_test, amount=amount_test)
print(f"Benefit/cost of our tuned model: ${benefit_cost:,.2f}")

# %%
balanced_accuracy = balanced_accuracy_scorer(tuned_model, data_test, target_test)
print(f"Balanced accuracy of our tuned model: {balanced_accuracy:.2f}")

# %%
tuned_model.set_params(objective_metric="balanced_accuracy").fit(
    data_train, target_train
)
balanced_accuracy = balanced_accuracy_scorer(tuned_model, data_test, target_test)
print(f"Balanced accuracy of our tuned model: {balanced_accuracy:.2f}")

# %%
logistic_regression.set_params(class_weight="balanced").fit(data_train, target_train)
balanced_accuracy = balanced_accuracy_scorer(
    logistic_regression, data_test, target_test
)
print(f"Balanced accuracy of our logistic regression: {balanced_accuracy:.2f}")

# %%
from sklearn.calibration import CalibrationDisplay

CalibrationDisplay.from_estimator(
    logistic_regression, data_test, target_test, n_bins=10
)

# %%
CalibrationDisplay.from_estimator(tuned_model, data_test, target_test, n_bins=10)

# %%
from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(n_jobs=-1).fit(data_train, target_train)
# %%
balanced_accuracy = balanced_accuracy_scorer(rf, data_test, target_test)
print(f"Balanced accuracy of our random forest: {balanced_accuracy:.2f}")
benefit_cost = business_scorer(rf, data_test, target_test, amount=amount_test)
print(f"Benefit/cost of our random forest: ${benefit_cost:,.2f}")

# %%
from imblearn.ensemble import BalancedRandomForestClassifier

brf = BalancedRandomForestClassifier(
    sampling_strategy="all", replacement=True, bootstrap=False, n_jobs=-1
)
brf.fit(data_train, target_train)

# %%
balanced_accuracy = balanced_accuracy_scorer(brf, data_test, target_test)
print(f"Balanced accuracy of our balanced random forest: {balanced_accuracy:.2f}")
benefit_cost = business_scorer(brf, data_test, target_test, amount=amount_test)
print(f"Benefit/cost of our balanced random forest: ${benefit_cost:,.2f}")

# %%
tuned_model = TunedThresholdClassifier(
    estimator=brf,
    objective_metric=business_scorer,
    n_thresholds=1_000,
    n_jobs=-1,
)
tuned_model.fit(data_train, target_train, amount=amount_train)

# %%
balanced_accuracy = balanced_accuracy_scorer(tuned_model, data_test, target_test)
print(f"Balanced accuracy of our balanced random forest: {balanced_accuracy:.2f}")
benefit_cost = business_scorer(tuned_model, data_test, target_test, amount=amount_test)
print(f"Benefit/cost of our balanced random forest: ${benefit_cost:,.2f}")

# %%
tuned_model = TunedThresholdClassifier(
    estimator=rf,
    objective_metric=business_scorer,
    n_thresholds=1_000,
    n_jobs=-1,
)
tuned_model.fit(data_train, target_train, amount=amount_train)

# %%
balanced_accuracy = balanced_accuracy_scorer(tuned_model, data_test, target_test)
print(f"Balanced accuracy of our balanced random forest: {balanced_accuracy:.2f}")
benefit_cost = business_scorer(tuned_model, data_test, target_test, amount=amount_test)
print(f"Benefit/cost of our balanced random forest: ${benefit_cost:,.2f}")

# %%
