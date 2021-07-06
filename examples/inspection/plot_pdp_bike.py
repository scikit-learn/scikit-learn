# %%
import pandas as pd

bikes = pd.read_csv(
    "https://raw.githubusercontent.com/christophM/interpretable-ml-book/master"
    "/data/bike.csv"
)
bikes.head()

# %%
target_name = "cnt"
X = bikes.drop(columns=[target_name, "days_since_2011"])
y = bikes[target_name]

# %%
from sklearn.model_selection import train_test_split

y -= y.mean()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)

# %%
# Create a column transformer that will preprocess the numerical data with
# standard scaler and the categorical data with one hot encoder.
from sklearn.compose import ColumnTransformer
from sklearn.compose import make_column_selector as selector
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import OneHotEncoder

numerical_features = selector(dtype_include=["float", "int"])
categorical_features = selector(dtype_include=["object"])

preprocessor = ColumnTransformer(
    transformers=[
        ("num", StandardScaler(), numerical_features),
        ("cat", OneHotEncoder(handle_unknown="ignore"), categorical_features),
    ]
)

# %%
# Create a pipeline that will preprocess the data and then use a MLP regressor
# to predict the target.
from time import time
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline

print("Training MLPRegressor...")
tic = time()
mlp_model = make_pipeline(
    preprocessor,
    MLPRegressor(
        hidden_layer_sizes=(100, 50), learning_rate_init=0.01, early_stopping=True
    ),
)
mlp_model.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {mlp_model.score(X_test, y_test):.2f}")

# %%
from sklearn.inspection import plot_partial_dependence

print("Computing partial dependence plots...")
features = numerical_features(X_train)
tic = time()
display = plot_partial_dependence(
    mlp_model,
    X_train,
    features,
    kind="both",
    subsample=50,
    n_jobs=-1,
    grid_resolution=20,
    random_state=0,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "XXX",
)
display.figure_.subplots_adjust(wspace=0.4, hspace=0.3)

# %%
# Create a pipeline that will preprocess the data and then use a Histogram
# graph regressor to predict the target.
# The preprocessor will only be made of an ordinal encoder for the categorical
# data.
from sklearn.preprocessing import OrdinalEncoder

preprocessor = ColumnTransformer(
    transformers=[
        ("cat", OrdinalEncoder(), categorical_features),
    ],
    remainder="passthrough",
    sparse_threshold=1,
)

# %%
from sklearn.ensemble import HistGradientBoostingRegressor

print("Training HistGradientBoostingRegressor...")
tic = time()
hgbdt_model = make_pipeline(preprocessor, HistGradientBoostingRegressor())
hgbdt_model.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {hgbdt_model.score(X_test, y_test):.2f}")

# %%
print("Computing partial dependence plots...")
features = numerical_features(X_train)
tic = time()
display = plot_partial_dependence(
    hgbdt_model,
    X_train,
    features,
    kind="both",
    subsample=50,
    n_jobs=-1,
    grid_resolution=20,
    random_state=0,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "XXX",
)
display.figure_.subplots_adjust(wspace=0.4, hspace=0.3)

# %%
