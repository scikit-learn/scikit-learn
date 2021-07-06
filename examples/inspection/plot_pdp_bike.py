# %%
from sklearn.datasets import fetch_openml

X, y = fetch_openml("Bike_Sharing_Demand", version=2, as_frame=True, return_X_y=True)

# %%
X.head()

# %%
y.head()

# %%
from sklearn.model_selection import train_test_split

y -= y.mean()
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.1, shuffle=False, random_state=0
)

# %%
numerical_features = [
    "temp",
    "feel_temp",
    "humidity",
    "windspeed",
]
categorical_features = X_train.columns.drop(numerical_features)

# %%
# Create a column transformer that will preprocess the numerical data with
# standard scaler and the categorical data with one hot encoder.
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import OneHotEncoder

mlp_preprocessor = ColumnTransformer(
    transformers=[
        ("num", QuantileTransformer(n_quantiles=100), numerical_features),
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
    mlp_preprocessor,
    MLPRegressor(
        hidden_layer_sizes=(50, 50), learning_rate_init=0.01, early_stopping=True
    ),
)
mlp_model.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {mlp_model.score(X_test, y_test):.2f}")

# %%
from sklearn.inspection import plot_partial_dependence

print("Computing partial dependence plots...")
features = ["temp", "humidity", "windspeed"]
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
display.figure_.subplots_adjust(wspace=0.4, hspace=0.5)

# %%
# Create a pipeline that will preprocess the data and then use a Histogram
# graph regressor to predict the target.
# The preprocessor will only be made of an ordinal encoder for the categorical
# data.
from sklearn.preprocessing import OrdinalEncoder

hgbdt_preprocessor = ColumnTransformer(
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
hgbdt_model = make_pipeline(hgbdt_preprocessor, HistGradientBoostingRegressor())
hgbdt_model.fit(X_train, y_train)
print(f"done in {time() - tic:.3f}s")
print(f"Test R2 score: {hgbdt_model.score(X_test, y_test):.2f}")

# %%
print("Computing partial dependence plots...")
features = ["temp", "humidity", "windspeed"]
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
display.figure_.subplots_adjust(hspace=0.5)

# %%
print("Computing partial dependence plots...")
features = ["season", "weather"]
tic = time()
display = plot_partial_dependence(
    hgbdt_model,
    X_train,
    features,
    kind="average",
    subsample=50,
    n_jobs=-1,
    grid_resolution=20,
    random_state=0,
    categorical_features=categorical_features,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "XXX",
)
display.figure_.subplots_adjust(hspace=1.0)

# %%
print("Computing partial dependence plots...")
features = ["temp", "humidity", ("temp", "humidity")]
tic = time()
display = plot_partial_dependence(
    hgbdt_model,
    X_train,
    features,
    kind="average",
    subsample=50,
    n_jobs=-1,
    grid_resolution=20,
    random_state=0,
)
print(f"done in {time() - tic:.3f}s")
display.figure_.suptitle(
    "XXX",
)
display.figure_.subplots_adjust(wspace=0.5, hspace=0.5)

# %%
