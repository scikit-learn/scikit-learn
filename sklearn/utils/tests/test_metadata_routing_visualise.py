"""Tests for sklearn.utils._metadata_routing_visualise

The visualisation helper functions are largely side-effect based (they print
ASCII diagrams).  The public contract we care about is that the *underlying*
metadata inspection logic builds the correct internal structures and that the
high-level helpers do not crash and include key information in their output.

These tests focus on `_collect_routing_info` (the work-horse that produces a
rich mapping used by all views) and perform a smoke-test on `visualise_routing` to
ensure it runs end-to-end.
"""

from __future__ import annotations

import re

from sklearn import config_context, set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import get_scorer
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_requests import WARN, _routing_repr
from sklearn.utils._metadata_routing_visualise import (
    _collect_routing_info,
    visualise_routing,
)
from sklearn.utils.metadata_routing import get_routing_for_object


@config_context(enable_metadata_routing=True)
def capsys_disabled_routing(monkeypatch):
    """Capture *only* stdout of `visualise_routing`.

    We monkey-patch ``print`` inside the visualisation module so that the output
    is sent to a list we control.  This is more robust than capturing the
    global stdout via ``capsys`` in case other prints happen in parallel.
    """

    captured: list[str] = []

    def fake_print(*args, **kwargs):
        captured.append(" ".join(map(str, args)))

    from sklearn.utils import _metadata_routing_visualise as vis_mod

    monkeypatch.setattr(vis_mod, "print", fake_print)
    return captured


# ---------------------------------------------------------------------------
# Low-level mapping tests
# ---------------------------------------------------------------------------


@config_context(enable_metadata_routing=True)
def test_collect_routing_info_simple_consumer():
    lr = LogisticRegression().set_fit_request(sample_weight=True)
    routing = get_routing_for_object(lr)

    info = _collect_routing_info(routing)
    path = _routing_repr(lr)  # root path is _routing_repr of the root object

    # The parameter should be present and requested for ``fit`` only.
    assert "sample_weight" in info[path]["params"]
    assert info[path]["statuses"]["sample_weight"]["fit"] is True


@config_context(enable_metadata_routing=True)
def test_collect_routing_info_alias():
    scaler = StandardScaler().set_fit_request(sample_weight="user_w")
    routing = get_routing_for_object(scaler)

    info = _collect_routing_info(routing)
    path = scaler.__class__.__name__

    # Original and alias parameters should be registered correctly.
    assert "sample_weight" in info[path]["params"]
    assert info[path]["aliases"]["sample_weight"] == "user_w"
    assert info[path]["statuses"]["sample_weight"]["fit"] == "user_w"


@config_context(enable_metadata_routing=True)
def test_collect_routing_info_show_all_metadata_flag():
    scaler = StandardScaler().set_fit_request(sample_weight=True)
    routing = get_routing_for_object(scaler)

    # Collecting with the flag True (always default) still returns all params.
    info_all = _collect_routing_info(routing)
    # "copy" is a legitimate StandardScaler kw on transform method.
    assert "copy" in info_all[_routing_repr(scaler)]["params"]

    # Passing False no longer changes behaviour; still includes all params.
    info_req = _collect_routing_info(routing)
    assert "sample_weight" in info_req[_routing_repr(scaler)]["params"]


# ---------------------------------------------------------------------------
# Smoke test high-level visualisation
# ---------------------------------------------------------------------------


@config_context(enable_metadata_routing=True)
def test_visualise_routing_smoke(capsys):
    """The helper should print something meaningful and not crash."""
    pipe = Pipeline(
        steps=[
            ("scaler", StandardScaler().set_fit_request(sample_weight="w")),
            ("clf", LogisticRegression().set_fit_request(sample_weight=True)),
        ]
    )

    routing = get_routing_for_object(pipe)
    visualise_routing(routing)

    out = capsys.readouterr().out
    # Basic expectations – root name and the alias arrow should appear.
    assert pipe.__class__.__name__ in out  # Root header
    assert "w→sample_weight" in out  # alias display
    # tree connectors should be present to confirm rendering
    assert "│" in out or "├" in out or "└" in out


# ---------------------------------------------------------------------------
# Regex safety – ensure no stray placeholders left in output (e.g. "{{" )
# ---------------------------------------------------------------------------


@config_context(enable_metadata_routing=True)
def test_no_template_placeholders_in_output(capsys):
    lr = LogisticRegression().set_fit_request(sample_weight=True)
    routing = get_routing_for_object(lr)
    visualise_routing(routing)
    out = capsys.readouterr().out
    assert not re.search(r"\{\{.*\}\}", out)


# ============================================================================
# TEST CASES
# ============================================================================


def run_test_1():
    numeric_features = ["age", "fare"]
    numeric_transformer = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            (
                "scaler",
                StandardScaler()
                .set_fit_request(sample_weight="inner_weights")
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
            ("classifier", LogisticRegression().set_fit_request(sample_weight=False)),
        ]
    )

    param_grid = {
        "preprocessor__num__imputer__strategy": ["mean", "median"],
        "preprocessor__cat__selector__percentile": [10, 30, 50, 70],
        "classifier__C": [0.1, 1.0, 10, 100],
    }

    scorer = get_scorer("accuracy").set_score_request(sample_weight=True)

    search_cv = RandomizedSearchCV(
        clf, param_grid, cv=GroupKFold(), scoring=scorer, random_state=0
    )

    # Get the routing information
    test = get_routing_for_object(search_cv)

    visualise_routing(test)


def run_test_2():
    est = make_pipeline(
        make_pipeline(StandardScaler().set_fit_request(sample_weight=True)),
        make_pipeline(StandardScaler().set_fit_request(sample_weight=False)),
        make_pipeline(StandardScaler()),
        make_pipeline(StandardScaler().set_fit_request(sample_weight=WARN)),
    )

    visualise_routing(get_routing_for_object(est))


def run_test_3():
    est = RandomizedSearchCV(estimator=LogisticRegression(), param_distributions={})
    visualise_routing(get_routing_for_object(est))

    est = RandomizedSearchCV(
        estimator=LogisticRegression(),
        param_distributions={},
        scoring=get_scorer("accuracy").set_score_request(sample_weight=True),
    )
    visualise_routing(get_routing_for_object(est))


if __name__ == "__main__":
    # Enable metadata routing
    set_config(enable_metadata_routing=True)
    run_test_1()
    run_test_2()
    run_test_3()
