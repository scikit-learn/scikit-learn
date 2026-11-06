"""
==================================
Visualising metadata routing trees
==================================

.. currentmodule:: sklearn

When metadata routing is enabled and a composite estimator routes metadata such
as ``sample_weight`` or ``groups`` through a hierarchy of pipelines, column
transformers, scorers and cross-validation splitters, it can be hard to tell
which parameter is consumed where, which is renamed (aliased), and which raises
or warns.

:func:`~sklearn.utils.metadata_routing.get_routing_for_object` returns the
raw routing object, but its default ``repr`` is not designed to be read at a
glance. The helpers :func:`~sklearn.utils._metadata_routing_visualise.format_routing`
and :func:`~sklearn.utils._metadata_routing_visualise.visualise_routing` render
it as a human-readable tree with a per-parameter summary.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Build a representative estimator
# --------------------------------
#
# A :class:`~sklearn.model_selection.RandomizedSearchCV` wraps a
# :class:`~sklearn.pipeline.Pipeline` whose preprocessing step is a
# :class:`~sklearn.compose.ColumnTransformer`. The CV uses
# :class:`~sklearn.model_selection.GroupKFold` (so ``groups`` must be routed)
# and a scorer with an aliased ``sample_weight`` (so ``my_w`` must be routed
# to the scorer's ``sample_weight``).

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import get_scorer
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_routing_visualise import visualise_routing
from sklearn.utils.metadata_routing import get_routing_for_object

set_config(enable_metadata_routing=True)

numeric = Pipeline(
    steps=[
        ("imputer", SimpleImputer(strategy="median")),
        (
            "scaler",
            StandardScaler()
            .set_fit_request(sample_weight=True)
            .set_transform_request(copy=True),
        ),
    ]
)
categorical = Pipeline(
    steps=[
        ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ("selector", SelectPercentile(chi2, percentile=50)),
    ]
)
preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric, ["age", "fare"]),
        ("cat", categorical, ["embarked", "sex", "pclass"]),
    ]
)
clf = Pipeline(
    steps=[
        ("preprocessor", preprocessor),
        ("classifier", LogisticRegression().set_fit_request(sample_weight=False)),
    ]
)
scorer = get_scorer("accuracy").set_score_request(sample_weight="my_w")
search = RandomizedSearchCV(clf, {}, cv=GroupKFold(), scoring=scorer, random_state=0)

routing = get_routing_for_object(search)

# %%
# Default view: legend, tree, and per-method summary
# --------------------------------------------------

visualise_routing(routing)

# %%
# Filter by parameter
# -------------------
#
# Restrict the output to a single user-facing parameter. The parameter name
# used for filtering is the *alias* (what the user would pass as a keyword
# when calling ``fit``/``score``), so filtering by ``"my_w"`` picks up the
# scorer even though the component name is ``sample_weight``.

visualise_routing(routing, param="my_w")

# %%
# Filter by root method
# ---------------------
#
# Restrict the summary to a single root method — for example, only what
# happens when ``search.fit(...)`` is called.

visualise_routing(routing, method="fit")
