# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn import config_context
from sklearn.callback._testing.callbacks import (
    RecordingCallback,
    SampleWeightCallback,
)
from sklearn.callback._testing.estimators import MaxIterEstimator, MetaEstimator


@config_context(enable_metadata_routing=True)
def test_metadata_routing_callback_consumer():
    """Test the metadata routing to consumer callbacks."""

    cb = (
        RecordingCallback()
        .set_on_fit_task_begin_request(requested_arg_begin="foo")
        .set_on_fit_task_end_request(requested_arg_end=True)
    )
    cb2 = (
        RecordingCallback()
        .set_on_fit_task_begin_request(requested_arg_begin="bar")
        .set_on_fit_task_end_request(requested_arg_end="sample_weight")
    )

    MaxIterEstimator().set_callbacks(cb, cb2).fit(
        foo="val_1", requested_arg_end="val_2", bar="bar", sample_weight="sample_weight"
    )

    for record, val_begin, val_end in [
        (cb.record, "val_1", "val_2"),
        (cb2.record, "bar", "sample_weight"),
    ]:
        task_begin_metadatas = [
            rec["kwargs"]["requested_arg_begin"]
            for rec in record
            if rec["name"] == "on_fit_task_begin"
        ]
        task_end_metadatas = [
            rec["kwargs"]["requested_arg_end"]
            for rec in record
            if rec["name"] == "on_fit_task_end"
        ]
        assert task_begin_metadatas
        assert task_end_metadatas
        assert all([m == val_begin for m in task_begin_metadatas])
        assert all([m == val_end for m in task_end_metadatas])


@config_context(enable_metadata_routing=True)
@pytest.mark.parametrize("n_jobs", [1, 2])
def test_metadata_routing_callback_consumer_in_metaestimator(n_jobs):
    """Test a callback that requests metadata in a meta-estimator."""

    cb = (
        RecordingCallback()
        .set_on_fit_task_begin_request(requested_arg_begin="foo")
        .set_on_fit_task_end_request(requested_arg_end=True)
    )

    estimator = MaxIterEstimator().set_callbacks(cb)
    MetaEstimator(estimator=estimator, n_jobs=n_jobs).fit(
        foo="val_1", requested_arg_end="val_2"
    )
    task_begin_metadatas = [
        rec["kwargs"]["requested_arg_begin"]
        for rec in cb.record
        if rec["name"] == "on_fit_task_begin"
    ]
    task_end_metadatas = [
        rec["kwargs"]["requested_arg_end"]
        for rec in cb.record
        if rec["name"] == "on_fit_task_end"
    ]
    assert task_begin_metadatas
    assert task_end_metadatas
    assert all([m == "val_1" for m in task_begin_metadatas])
    assert all([m == "val_2" for m in task_end_metadatas])


@config_context(enable_metadata_routing=True)
@pytest.mark.parametrize("with_callback", [True, False])
def test_callback_support_keeps_estimator_own_requests(with_callback):
    """Check that an estimator with callback support keeps its own requests.

    `MaxIterEstimator` consumes `sample_weight` and relies on the routing of
    `CallbackSupportMixin`: the mixin must not drop the requests of the estimator
    itself, otherwise a meta-estimator would not route `sample_weight` to it.
    """
    estimator = MaxIterEstimator().set_fit_request(sample_weight=True)
    if with_callback:
        estimator.set_callbacks(RecordingCallback())

    routing = estimator.get_metadata_routing()
    assert routing.consumes("fit", ["sample_weight"]) == {"sample_weight"}


def _assert_recorded_routed_arg(callback, expected):
    for hook_name, key in (
        ("on_fit_task_begin", "requested_arg_begin"),
        ("on_fit_task_end", "requested_arg_end"),
    ):
        values = [
            rec["kwargs"][key] for rec in callback.record if rec["name"] == hook_name
        ]
        assert values
        assert all(value == expected for value in values)


def test_sample_weight_without_routing():
    """Forward sample_weight without metadata routing only to callbacks that accept it.

    Fit goes through both manual routing steps. `_get_manual_callback_params` builds
    one mapping shared by every callback, then `_get_manual_routing_params` delivers
    `sample_weight` only to callbacks that accept it. Other callbacks are still called.
    """
    weight = "sample_weight"
    accepting = SampleWeightCallback()
    ignoring = RecordingCallback()
    MaxIterEstimator(max_iter=1, computation_intensity=0).set_callbacks(
        ignoring, accepting
    ).fit(sample_weight=weight)
    _assert_recorded_routed_arg(accepting, weight)
    _assert_recorded_routed_arg(ignoring, None)

    # Nothing to forward: no sample_weight, or no callback accepts it.
    for callbacks, fit_kwargs in (
        ([SampleWeightCallback()], {}),
        ([RecordingCallback()], {"sample_weight": weight}),
    ):
        (callback,) = callbacks
        MaxIterEstimator(max_iter=1, computation_intensity=0).set_callbacks(
            *callbacks
        ).fit(**fit_kwargs)
        _assert_recorded_routed_arg(callback, None)
