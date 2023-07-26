# %%
import pandas as pd

filename_training = "~/Downloads/kdd98/epsilon_mirror/cup98lrn.zip"
index = ["CONTROLN"]
df_train = pd.read_csv(
    filename_training, compression="zip", encoding="latin-1", low_memory=False
).set_index(index)
train_indices = df_train.index

# %%
filename_data_test = "~/Downloads/kdd98/epsilon_mirror/cup98val.zip"
data_test = pd.read_csv(
    filename_data_test, compression="zip", encoding="latin-1", low_memory=False
).set_index(index)

# %%
filename_target_test = "~/Downloads/kdd98/epsilon_mirror/valtargt.txt"
target_test = pd.read_csv(filename_target_test).set_index(index)

# %%
df_test = pd.concat([data_test, target_test], axis=1)
test_indices = df_test.index
# Do not convert string to avoid the pd.NA bug:
# xref: https://github.com/scikit-learn/scikit-learn/issues/26890
df = pd.concat([df_train, df_test], axis=0)
# df = df.convert_dtypes(convert_string=False)

# %%
# convert to categorical in case we don't want to use the TableVectorizer
categorical_columns = {}
for col_idx, col in enumerate(df.columns):
    series = df[col]
    dtype_series = series.dtype
    # if isinstance(dtype_series, pd.StringDtype):
    if dtype_series.kind == "O":
        categorical_columns[col] = "category"
        continue  # skip the rest of the loop
    unique_values = series.value_counts()
    if len(unique_values) < 60:
        # low-cardinality features are considered as categorical
        categorical_columns[col] = "category"

df = df.astype(categorical_columns)

# %%
target_continuous_name = "TARGET_D"
target_binary_name = "TARGET_B"
neg_label, pos_label = 0, 1

# %%
from sklearn.model_selection import train_test_split

# df_train, df_test = train_test_split(df, test_size=0.5, random_state=42)
df_train = df.loc[train_indices]
df_test = df.loc[test_indices]

# %%
data_train = df_train.drop(columns=[target_continuous_name, target_binary_name])

# %%
data_train

# %%
target_continuous_train = df_train[target_continuous_name]
target_binary_train = df_train[target_binary_name]

# %%


def cost_metric(y_true, y_pred, neg_label, pos_label, donation_amount):
    """Compute the business cost related to the prediction.

    The real cost is to send a mail to a person and is evaluated to $0.68.
    The gain is the donation amount if the person donate.
    """
    mask_true_positive = (y_true == pos_label) & (y_pred == pos_label)
    mask_false_positive = (y_true == neg_label) & (y_pred == pos_label)
    # mask_false_negative = (y_true == pos_label) & (y_pred == neg_label)
    # mask_true_negative = (y_true == neg_label) & (y_pred == neg_label)
    cost_sending_mail = -0.68
    cost_false_positive = cost_sending_mail * mask_false_positive.sum()
    gain_true_positive = (cost_sending_mail + donation_amount[mask_true_positive]).sum()
    return gain_true_positive + cost_false_positive
    # loss_false_negative = donation_amount[mask_false_negative].sum()
    # fp = - mask_false_positive.sum() * cost_sending_mail
    # fn = - (donation_amount[mask_false_negative] - cost_sending_mail).sum()
    # return fp + fn


# %%
import numpy as np

gain = cost_metric(
    y_true=target_binary_train,
    y_pred=np.zeros_like(target_binary_train),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_train,
)
print(f"Gain if we don't send any mail: ${gain:,.2f}")

# %%
gain = cost_metric(
    y_true=target_binary_train,
    y_pred=np.ones_like(target_binary_train),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_train,
)
print(f"Gain if we send mails to everyone: ${gain:,.2f}")

# %%
gain = cost_metric(
    y_true=target_binary_train,
    y_pred=target_binary_train,
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_train,
)
print(f"Maximum gain on the training set: ${gain:,.2f}")

# %%
import sklearn
from sklearn.metrics import make_scorer

sklearn.set_config(enable_metadata_routing=True)

cost_scorer = make_scorer(
    cost_metric, neg_label=neg_label, pos_label=pos_label
).set_score_request(donation_amount=True)

# %%
from sklearn.dummy import DummyClassifier

stingy_classifier = DummyClassifier(strategy="constant", constant=0)
stingy_classifier.fit(data_train, target_binary_train)

# %%
bling_bling_classifier = DummyClassifier(strategy="constant", constant=1)
bling_bling_classifier.fit(data_train, target_binary_train)

# %%
gain = cost_scorer(
    stingy_classifier,
    data_train,
    target_binary_train,
    donation_amount=target_continuous_train,
)
print(f"Gain of the stingy classifier: ${gain:,.2f}")

# %%
gain = cost_scorer(
    bling_bling_classifier,
    data_train,
    target_binary_train,
    donation_amount=target_continuous_train,
)
print(f"Gain of the bling-bling classifier: ${gain:,.2f}")

# %%
data_test = df_test.drop(columns=[target_continuous_name, target_binary_name])

# %%
target_continuous_test = df_test[target_continuous_name]
target_binary_test = df_test[target_binary_name]

# %%
from sklearn.compose import ColumnTransformer
from sklearn.compose import make_column_selector as selector
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OrdinalEncoder

sklearn.set_config(transform_output="pandas")

categorical_columns = selector(dtype_include="category")
preprocessing = ColumnTransformer(
    transformers=[
        (
            "categorical",
            OrdinalEncoder(
                handle_unknown="use_encoded_value",
                unknown_value=-1,
                min_frequency=0.05,
                max_categories=255,
            ),
            categorical_columns,
        ),
    ],
    remainder="passthrough",
    n_jobs=-1,
    verbose_feature_names_out=False,
)
model = Pipeline(
    steps=[
        ("preprocessing", preprocessing),
        (
            "classifier",
            HistGradientBoostingClassifier(
                max_iter=1_000,
                early_stopping=True,
                categorical_features=categorical_columns(data_train),
            ),
            # DecisionTreeClassifier(max_leaf_nodes=200, random_state=42),
        ),
    ]
)
model.fit(data_train, target_binary_train)

# %%
model.score(data_test, target_binary_test)

# %%
from sklearn.metrics import balanced_accuracy_score

balanced_accuracy_score(target_binary_test, model.predict(data_test))

# %%
gain = cost_metric(
    y_true=target_binary_test,
    y_pred=model.predict(data_test),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_test,
)
print(f"Gain of the model: ${gain:,.2f}")

# %%
from sklearn.calibration import CalibrationDisplay

CalibrationDisplay.from_estimator(model, data_test, target_binary_test, n_bins=10)

# %%
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, TunedThresholdClassifier

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
tuned_model = TunedThresholdClassifier(
    # estimator=model,
    estimator=CalibratedClassifierCV(model, cv=5, method="isotonic"),
    n_thresholds=np.linspace(0.02, 0.05, 1_000),
    # n_thresholds=1_000,
    # strategy="constant",
    # constant_threshold=0.0302,
    pos_label=pos_label,
    objective_metric=cost_scorer,
    cv=cv,
    n_jobs=-1,
)
tuned_model.fit(
    data_train, target_binary_train, donation_amount=target_continuous_train
)

# %%
CalibrationDisplay.from_estimator(tuned_model, data_test, target_binary_test, n_bins=10)

# %%
tuned_model.score(data_test, target_binary_test)

# %%
balanced_accuracy_score(target_binary_test, tuned_model.predict(data_test))

# %%
gain = cost_metric(
    y_true=target_binary_test,
    y_pred=tuned_model.predict(data_test),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_test,
)
print(f"Gain of the model: ${gain:,.2f}")

# %%
import matplotlib.pyplot as plt

for thresholds, scores in zip(tuned_model.cv_thresholds_, tuned_model.cv_scores_):
    plt.semilogx(thresholds, scores, alpha=0.5)

# %%
import plotly.graph_objects as go

fig = go.Figure()
for thresholds, scores in zip(tuned_model.cv_thresholds_, tuned_model.cv_scores_):
    fig.add_trace(go.Scatter(x=thresholds, y=scores, mode="lines", name="lines"))
fig.update_xaxes(type="log")
fig.show()

# %%

# %%
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer

categorical_columns = selector(dtype_include="category")
preprocessing = ColumnTransformer(
    transformers=[
        (
            "categorical",
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1),
            categorical_columns,
        ),
    ],
    remainder="passthrough",
    n_jobs=-1,
    verbose_feature_names_out=False,
)
model = Pipeline(
    steps=[
        ("preprocessing", preprocessing),
        ("imputer", SimpleImputer(strategy="mean")),
        ("classifier", RandomForestClassifier(n_jobs=-1)),
    ]
)
tuned_model = TunedThresholdClassifier(
    estimator=CalibratedClassifierCV(model, cv=5, method="isotonic"),
    pos_label=pos_label,
    objective_metric=cost_scorer,
    cv=0.8,
)
tuned_model.fit(
    data_train, target_binary_train, donation_amount=target_continuous_train
)

# %%
balanced_accuracy_score(target_binary_test, tuned_model.predict(data_test))

# %%
gain = cost_metric(
    y_true=target_binary_test,
    y_pred=tuned_model.predict(data_test),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_test,
)
print(f"Gain of the model: ${gain:,.2f}")

# %%
from imblearn.ensemble import BalancedRandomForestClassifier

sklearn.set_config(enable_metadata_routing=False)

categorical_columns = selector(dtype_include="category")
preprocessing = ColumnTransformer(
    transformers=[
        (
            "categorical",
            OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1),
            categorical_columns,
        ),
    ],
    remainder="passthrough",
    n_jobs=-1,
    verbose_feature_names_out=False,
)
model = Pipeline(
    steps=[
        ("preprocessing", preprocessing),
        ("imputer", SimpleImputer(strategy="mean")),
        (
            "classifier",
            BalancedRandomForestClassifier(
                sampling_strategy="all", replacement=True, bootstrap=False, n_jobs=-1
            ),
        ),
    ]
)

# %%
balanced_accuracy_score(target_binary_test, model.predict(data_test))

# %%
gain = cost_metric(
    y_true=target_binary_test,
    y_pred=model.predict(data_test),
    neg_label=neg_label,
    pos_label=pos_label,
    donation_amount=target_continuous_test,
)
print(f"Gain of the model: ${gain:,.2f}")


# %%
# %%
# %%
from sklearn.datasets import fetch_openml

credit_card = fetch_openml(data_id=1597, as_frame=True, parser="pandas")
credit_card.frame

# %%
data = credit_card.frame.drop(columns=["Class", "Amount"])
target = credit_card.frame["Class"].astype(int).to_numpy()
amount = credit_card.frame["Amount"].to_numpy()

# %%
target.value_counts()

# %%
target.value_counts(normalize=True)

# %%
fraud = target == 1
amount_fraud = amount[fraud]
ax = amount_fraud.plot.hist(bins=100)
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
