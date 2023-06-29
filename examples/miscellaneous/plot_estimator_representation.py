"""
===========================================
Displaying estimators and complex pipelines
===========================================

This example illustrates different ways estimators and pipelines can be
displayed.
"""

from sklearn.compose import make_column_transformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

# %%
# Compact text representation
# ---------------------------
#
# Estimators will only show the parameters that have been set to non-default
# values when displayed as a string. This reduces the visual noise and makes it
# easier to spot what the differences are when comparing instances.

lr = LogisticRegression(penalty="l1")
print(lr)

# %%
# Rich HTML representation
# ------------------------
# In notebooks estimators and pipelines will use a rich HTML representation.
# This is particularly useful to summarise the
# structure of pipelines and other composite estimators, with interactivity to
# provide detail.  Click on the example image below to expand Pipeline
# elements.  See :ref:`visualizing_composite_estimators` for how you can use
# this feature.

num_proc = make_pipeline(SimpleImputer(strategy="median"), StandardScaler())

cat_proc = make_pipeline(
    SimpleImputer(strategy="constant", fill_value="missing"),
    OneHotEncoder(handle_unknown="ignore"),
)

preprocessor = make_column_transformer(
    (num_proc, ("feat1", "feat3")), (cat_proc, ("feat0", "feat2"))
)

clf = make_pipeline(preprocessor, LogisticRegression())
clf
