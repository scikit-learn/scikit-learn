# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for ``sklearn.utils._metadata_routing_visualise``."""

import pytest

from sklearn import config_context
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import get_scorer
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_requests import UNUSED, WARN, MetadataRouter
from sklearn.utils._metadata_routing_visualise import (
    _build_tree,
    _status_category,
    format_routing,
    visualise_routing,
)
from sklearn.utils.metadata_routing import get_routing_for_object


@pytest.fixture(autouse=True)
def _enable_metadata_routing():
    with config_context(enable_metadata_routing=True):
        yield


def _routing_of(est):
    return get_routing_for_object(est)


# --- Classification --------------------------------------------------------


def test_status_category_maps_known_values():
    assert _status_category(True) == "requested"
    assert _status_category(False) == "ignored"
    assert _status_category(None) == "errors"
    assert _status_category(WARN) == "warns"
    assert _status_category(UNUSED) == "unused"
    # Aliases are also "requested"
    assert _status_category("my_alias") == "requested"


# --- Traversal: node tree structure ----------------------------------------


def test_build_tree_single_consumer_sets_status():
    lr = LogisticRegression().set_fit_request(sample_weight=True)
    root = _build_tree(_routing_of(lr))
    assert root.owner_repr == "LogisticRegression"
    assert root.step_name is None
    assert root.children == []
    assert root.statuses["sample_weight"]["fit"] is True


def test_build_tree_alias_recorded_as_string_status():
    scaler = StandardScaler().set_fit_request(sample_weight="user_w")
    root = _build_tree(_routing_of(scaler))
    assert root.statuses["sample_weight"]["fit"] == "user_w"


def test_build_tree_nested_pipeline_depth():
    inner = Pipeline([("sc", StandardScaler().set_fit_request(sample_weight=True))])
    outer = Pipeline([("inner", inner)])
    root = _build_tree(_routing_of(outer))
    assert len(root.children) == 1
    assert root.children[0].step_name == "inner"
    assert root.children[0].children[0].step_name == "sc"


def test_build_tree_column_transformer_branches():
    ct = ColumnTransformer(
        [
            ("num", StandardScaler().set_fit_request(sample_weight=True), ["a"]),
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["b"]),
        ]
    )
    root = _build_tree(_routing_of(ct))
    names = [c.step_name for c in root.children]
    assert "num" in names and "cat" in names


def test_build_tree_origins_include_indirect_roots():
    # A scorer inside a RandomizedSearchCV is reachable from *both* ``fit``
    # (the CV loop calls the scorer while fitting) and ``score``.
    est = RandomizedSearchCV(
        LogisticRegression(),
        {},
        scoring=get_scorer("accuracy").set_score_request(sample_weight=True),
    )
    root = _build_tree(_routing_of(est))
    scorer = next(c for c in root.children if c.step_name == "scorer")
    roots = scorer.origins[("sample_weight", "score")]
    assert "fit" in roots
    assert "score" in roots


# --- End-to-end output ------------------------------------------------------


def test_format_routing_returns_string():
    out = format_routing(_routing_of(LogisticRegression()))
    assert isinstance(out, str)
    assert out  # non-empty


def test_visualise_routing_prints(capsys):
    visualise_routing(_routing_of(LogisticRegression()))
    captured = capsys.readouterr()
    assert "LogisticRegression" in captured.out
    assert captured.err == ""


def test_legend_present_by_default_and_togglable():
    routing = _routing_of(LogisticRegression())
    with_legend = format_routing(routing)
    without_legend = format_routing(routing, legend=False)
    assert "Legend:" in with_legend
    assert "Legend:" not in without_legend


def test_alias_arrow_in_tree():
    scaler = StandardScaler().set_fit_request(sample_weight="w")
    out = format_routing(_routing_of(scaler))
    assert "w → sample_weight" in out


def test_alias_uses_requested_glyph_not_a_dedicated_one():
    # An aliased request is always "requested"; no dedicated alias glyph is
    # needed because the ``alias → component`` prefix already carries the
    # renaming information.
    scaler = StandardScaler().set_fit_request(sample_weight="w")
    out = format_routing(_routing_of(scaler))
    assert "fit✓" in out
    assert "↗" not in out


def test_colour_off_by_default_in_format_routing():
    out = format_routing(_routing_of(LogisticRegression()))
    # No ANSI escape sequences.
    assert "\x1b[" not in out


def test_colour_flag_emits_ansi_escapes():
    out = format_routing(
        _routing_of(LogisticRegression().set_fit_request(sample_weight=True)),
        colour=True,
    )
    # Green escape wraps the requested glyph.
    assert "\x1b[32m✓\x1b[0m" in out


def test_colour_categories_match_severity():
    # Error, warn and requested should each get a distinct ANSI colour.
    routing = _routing_of(LogisticRegression().set_fit_request(sample_weight=WARN))
    out = format_routing(routing, colour=True)
    assert "\x1b[33m⚠\x1b[0m" in out  # yellow warn
    assert "\x1b[31m⛔\x1b[0m" in out  # red error (score is None by default)


def test_no_color_env_disables_auto_colour(monkeypatch, capsys):
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setattr("sys.stdout.isatty", lambda: True)
    visualise_routing(
        _routing_of(LogisticRegression().set_fit_request(sample_weight=True))
    )
    out = capsys.readouterr().out
    assert "\x1b[" not in out


def test_warn_glyph_rendered():
    est = make_pipeline(StandardScaler().set_fit_request(sample_weight=WARN))
    out = format_routing(_routing_of(est))
    assert "⚠" in out
    # And the summary should list it under "warns".
    assert "warns" in out


def test_unused_sentinel_never_surfaces_in_output():
    # ``UNUSED`` is only valid in class-level ``__metadata_request__*`` attrs,
    # where it means "remove this param" — entries tagged UNUSED are filtered
    # out before they reach the MethodMetadataRequest storage. The viz should
    # therefore never show the literal sentinel.
    out = format_routing(_routing_of(LogisticRegression()))
    assert UNUSED not in out


def test_self_request_rendered_at_router_path():
    # A router with a self-request should show its own consumed params on its
    # own node, not on a child.
    router = MetadataRouter(owner="my_router").add_self_request(
        StandardScaler().set_fit_request(sample_weight=True)
    )
    out = format_routing(router)
    assert "my_router" in out
    assert "sample_weight" in out


def test_summary_present_by_default_and_togglable():
    routing = _routing_of(
        LogisticRegression().set_fit_request(sample_weight=True),
    )
    with_summary = format_routing(routing)
    without_summary = format_routing(routing, summary=False)
    assert "Parameter summary:" in with_summary
    assert "Parameter summary:" not in without_summary


def test_tree_togglable():
    routing = _routing_of(
        LogisticRegression().set_fit_request(sample_weight=True),
    )
    # With tree off, there should be no root class line followed by indented
    # params — only the summary + legend.
    without_tree = format_routing(routing, tree=False)
    assert "LogisticRegression\n" not in without_tree
    assert "Parameter summary:" in without_tree


# --- Filtering --------------------------------------------------------------


def test_method_filter_restricts_summary_sections():
    est = RandomizedSearchCV(
        LogisticRegression().set_fit_request(sample_weight=True),
        {},
        scoring=get_scorer("accuracy").set_score_request(sample_weight=True),
    )
    routing = _routing_of(est)
    full = format_routing(routing)
    fit_only = format_routing(routing, method="fit")
    # ``score`` section appears by default but not under method="fit".
    assert "\nscore\n" in full
    assert "\nscore\n" not in fit_only
    # The fit section is still present.
    assert "\nfit\n" in fit_only


def test_method_filter_accepts_composite_names():
    # Passing ``"fit_transform"`` should be treated as its simple parts.
    est = make_pipeline(StandardScaler().set_fit_request(sample_weight=True))
    routing = _routing_of(est)
    out = format_routing(routing, method="fit_transform")
    assert "fit\n" in out  # fit summary section present


def test_param_filter_restricts_to_named_parameter():
    est = RandomizedSearchCV(
        LogisticRegression().set_fit_request(sample_weight=True),
        {},
        cv=GroupKFold(),
    )
    routing = _routing_of(est)
    only_groups = format_routing(routing, param="groups")
    assert "groups" in only_groups
    assert "sample_weight" not in only_groups


def test_param_filter_uses_alias_name():
    est = RandomizedSearchCV(
        LogisticRegression().set_fit_request(sample_weight=True),
        {},
        scoring=get_scorer("accuracy").set_score_request(sample_weight="my_w"),
    )
    routing = _routing_of(est)
    out = format_routing(routing, param="my_w")
    # Aliased consumption appears; non-aliased consumption of the same
    # component param on the estimator is filtered out.
    assert "my_w" in out
    assert "├─ ✓ sample_weight" not in out


# --- Large integration smoke -----------------------------------------------


def test_searchcv_pipeline_column_transformer_integration():
    """Big-case smoke: no crash and all expected landmarks show up."""
    num = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            (
                "scaler",
                StandardScaler()
                .set_fit_request(sample_weight="inner_weights")
                .set_transform_request(copy=True),
            ),
        ]
    )
    cat = Pipeline(
        [
            ("encoder", OneHotEncoder(handle_unknown="ignore")),
            ("selector", SelectPercentile(chi2, percentile=50)),
        ]
    )
    pre = ColumnTransformer([("num", num, ["a"]), ("cat", cat, ["b"])])
    clf = Pipeline(
        [
            ("preprocessor", pre),
            ("classifier", LogisticRegression().set_fit_request(sample_weight=False)),
        ]
    )
    scorer = get_scorer("accuracy").set_score_request(sample_weight=True)
    est = RandomizedSearchCV(clf, {}, cv=GroupKFold(), scoring=scorer, random_state=0)

    out = format_routing(_routing_of(est))

    # Landmarks from each branch of the tree.
    for token in [
        "RandomizedSearchCV",
        "preprocessor (ColumnTransformer)",
        "num (Pipeline)",
        "scaler (StandardScaler)",
        "classifier (LogisticRegression)",
        "splitter (GroupKFold)",
        "inner_weights → sample_weight",
        "groups[split✓]",
    ]:
        assert token in out, f"missing token: {token!r}"


@pytest.mark.parametrize("tree", [True, False])
@pytest.mark.parametrize("summary", [True, False])
@pytest.mark.parametrize("legend", [True, False])
def test_output_sections_are_independently_togglable(tree, summary, legend):
    routing = _routing_of(
        LogisticRegression().set_fit_request(sample_weight=True),
    )
    out = format_routing(routing, tree=tree, summary=summary, legend=legend)
    assert ("Legend:" in out) == legend
    assert ("Parameter summary:" in out) == summary
    assert ("LogisticRegression\n" in out) == tree
