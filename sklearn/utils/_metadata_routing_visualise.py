# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from io import StringIO
from pprint import pprint

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_requests import SIMPLE_METHODS
from sklearn.utils.metadata_routing import (
    MetadataRequest,
    MetadataRouter,
    get_routing_for_object,
)

# Enable metadata routing
set_config(enable_metadata_routing=True)


def _param_names(router):
    res = set()
    for method in SIMPLE_METHODS:
        res = res.union(
            router._get_param_names(
                method=method, return_alias=True, ignore_self_request=False
            )
        )
    return res


def visualise_routing(routing_info):
    params = _param_names(routing_info)
    for param in params:
        output = StringIO()
        print(f"Visualising routing for {param}", file=output)
        routed = _visualise_param(param, routing_info, indent=0, output=output)
        if routed:
            print(output.getvalue())


def _visualise_param(param, routing_info, indent, output):
    if isinstance(routing_info, MetadataRouter):
        return _visualise_metadata_router(param, routing_info, indent, output)
    elif isinstance(routing_info, MetadataRequest):
        return _visualise_param_in_metadata_request(
            param, routing_info, indent=indent + 1, output=output
        )
    else:
        raise ValueError("Unknown routing info type")


def _visualise_param_in_metadata_request(param, request, indent, output):
    print(" ." * indent, request.owner, file=output)
    routed = False
    for method in SIMPLE_METHODS:
        if request.consumes(method, [param]):
            print(" ." * (indent + 1), f"{method} consumes {param}", file=output)
            routed = True
    return routed


def _visualise_metadata_router(param, router, indent, output):
    routed = False
    print(" ." * indent, router.owner, file=output)
    indent += 1
    if router._self_request:
        sub_output = StringIO()
        sub_routed = _visualise_param_in_metadata_request(
            param, router._self_request, indent=indent + 1, output=sub_output
        )
        if sub_routed:
            print(" ." * indent, "Self request", file=output)
            print(sub_output.getvalue(), end="", file=output)
            routed = True

    for name, router_mapping_pair in router._route_mappings.items():
        inner_router = router_mapping_pair.router
        for mappings in router_mapping_pair.mapping:
            caller, callee = mappings.caller, mappings.callee
            if inner_router.consumes(method=callee, params=[param]):
                sub_output = StringIO()
                sub_routed = _visualise_param(
                    param, inner_router, indent=indent + 1, output=sub_output
                )
                if sub_routed:
                    if caller == "fit" and callee == "fit" and name == "classifier":
                        temp = 1
                    print(" ." * indent, f"{caller} -> {name}.{callee}", file=output)
                    print(sub_output.getvalue(), end="", file=output)
                    routed = True

    return routed


# test1 = get_routing_for_object(StandardScaler())
# print("empty StandardScaler")
# pprint(test1)
# visualise_routing(test1)

# test2 = get_routing_for_object(
#     StandardScaler()
#     .set_fit_request(sample_weight=True)
#     .set_inverse_transform_request(copy=False)
#     .set_transform_request(copy=True)
# )
# print("StandardScaler with requests")
# pprint(test2)
# visualise_routing(test2)

# test3 = get_routing_for_object(
#     Pipeline(steps=[("scaler", StandardScaler().set_fit_request(sample_weight=True))])
# )
# print("Pipeline with StandardScaler")
# pprint(test3._serialize())
# visualise_routing(test3)

# test4 = get_routing_for_object(
#     Pipeline(
#         steps=[
#             (
#                 "scaler",
#                 StandardScaler()
#                 .set_fit_request(sample_weight=True)
#                 .set_inverse_transform_request(copy=False)
#                 .set_transform_request(copy=True),
#             )
#         ]
#     )
# )
# print("Pipeline with StandardScaler with multiple requests")
# pprint(test4._serialize())
# visualise_routing(test4)

numeric_features = ["age", "fare"]
numeric_transformer = Pipeline(
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

categorical_features = ["embarked", "sex", "pclass"]
categorical_transformer = Pipeline(
    steps=[
        ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ("selector", SelectPercentile(chi2, percentile=50)),
    ]
)
preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, numeric_features),
        ("cat", categorical_transformer, categorical_features),
    ]
)

# %%
# Append classifier to preprocessing pipeline.
# Now we have a full prediction pipeline.
clf = Pipeline(
    steps=[
        ("preprocessor", preprocessor),
        ("classifier", LogisticRegression().set_fit_request(sample_weight=True)),
    ]
)

param_grid = {
    "preprocessor__num__imputer__strategy": ["mean", "median"],
    "preprocessor__cat__selector__percentile": [10, 30, 50, 70],
    "classifier__C": [0.1, 1.0, 10, 100],
}

search_cv = RandomizedSearchCV(clf, param_grid, cv=GroupKFold(), random_state=0)


# Get the routing information
test5 = get_routing_for_object(search_cv)
print("RandomizedSearchCV")
pprint(test5._serialize())
visualise_routing(test5)
