import re
from collections import namedtuple

import numpy as np
import pytest

from sklearn._config import config_context
from sklearn._engine import (
    convert_attributes,
    get_engine_classes,
    list_engine_provider_names,
)
from sklearn._engine.base import EngineSpec, _get_engine_classes, _parse_entry_point
from sklearn.base import EngineAwareMixin


class FakeDefaultEngine:
    pass


class FakeEngine:
    pass


class FakeEngineHolder:
    class NestedFakeEngine:
        pass


# Dummy classes used to test engine resolution
class DefaultEngine:
    engine_name = "test-engine"

    def __init__(self, estimator):
        self.estimator = estimator

    def accepts(self, X, y=None, sample_weight=None):
        return True


class NeverAcceptsEngine:
    engine_name = "test-engine"

    def __init__(self, estimator):
        self.estimator = estimator

    def accepts(self, X, y=None, sample_weight=None):
        return False


class AlwaysAcceptsEngine:
    engine_name = "test-engine"

    def __init__(self, estimator):
        self.estimator = estimator

    def accepts(self, X, y=None, sample_weight=None):
        return True


class AlsoAlwaysAcceptsEngine(AlwaysAcceptsEngine):
    pass


class FakeEstimator(EngineAwareMixin):
    _engine_name = "test-engine"
    _default_engine = DefaultEngine


FakeEntryPoint = namedtuple("FakeEntryPoint", ["name", "value"])


def test_parse_entry_point():
    fake_entry_point = FakeEntryPoint(
        name="fake_engine",
        value="sklearn._engine.tests.test_engines:FakeEngine",
    )
    spec = _parse_entry_point(fake_entry_point)
    assert spec.name == "fake_engine"
    assert spec.provider_name == "sklearn"  # or should it be scikit-learn?
    assert spec.get_engine_class() is FakeEngine


def test_parse_entry_point_for_nested_engine_class():
    fake_entry_point = FakeEntryPoint(
        name="nested_fake_engine",
        value="sklearn._engine.tests.test_engines:FakeEngineHolder.NestedFakeEngine",
    )
    spec = _parse_entry_point(fake_entry_point)
    assert spec.name == "nested_fake_engine"
    assert spec.provider_name == "sklearn"  # or should it be scikit-learn?
    assert spec.get_engine_class() is FakeEngineHolder.NestedFakeEngine


def test_list_engine_provider_names():
    provider_names = list_engine_provider_names()
    for provider_name in provider_names:
        assert isinstance(provider_name, str)


def test_get_engine_class_with_default():
    # Use config_context with an empty provider tuple to make sure that not provider
    # are available for test_missing_engine_name
    with config_context(engine_provider=()):
        engine_classes = list(
            get_engine_classes("test_missing_engine_name", default=FakeEngine)
        )
    assert engine_classes == [("default", FakeEngine)]


def test_get_engine_class():
    engine_specs = (
        EngineSpec(
            "kmeans", "provider3", "sklearn._engine.tests.test_engines", "FakeEngine"
        ),
        EngineSpec(
            "kmeans",
            "provider4",
            "sklearn._engine.tests.test_engines",
            "FakeEngineHolder.NestedFakeEngine",
        ),
    )

    engine_class = list(
        _get_engine_classes(
            engine_name="missing",
            provider_names=("provider1", "provider3"),
            engine_specs=engine_specs,
            default=FakeDefaultEngine,
        )
    )
    assert engine_class == [("default", FakeDefaultEngine)]

    engine_class = list(
        _get_engine_classes(
            engine_name="kmeans",
            provider_names=("provider3", "provider4"),
            engine_specs=engine_specs,
            default=FakeDefaultEngine,
        )
    )
    assert engine_class == [
        ("provider3", FakeEngine),
        ("provider4", FakeEngineHolder.NestedFakeEngine),
        ("default", FakeDefaultEngine),
    ]

    engine_class = list(
        _get_engine_classes(
            engine_name="kmeans",
            provider_names=("provider4", "provider3"),
            engine_specs=engine_specs,
            default=FakeDefaultEngine,
        )
    )
    assert engine_class == [
        ("provider4", FakeEngineHolder.NestedFakeEngine),
        ("provider3", FakeEngine),
        ("default", FakeDefaultEngine),
    ]

    engine_specs = engine_specs + (
        EngineSpec(
            "kmeans",
            "provider1",
            "sklearn.provider1.somewhere",
            "OtherEngine",
        ),
    )

    # Invalid imports are delayed until they are actually needed.
    engine_classes = _get_engine_classes(
        engine_name="kmeans",
        provider_names=("provider4", "provider3", "provider1"),
        engine_specs=engine_specs,
        default=FakeDefaultEngine,
    )

    next(engine_classes)
    next(engine_classes)
    with pytest.raises(ImportError, match=re.escape("sklearn.provider1")):
        next(engine_classes)


@pytest.mark.parametrize(
    "attribute_types,converted", [("sklearn_types", True), ("engine_types", False)]
)
def test_attribute_conversion(attribute_types, converted):
    """Test attribute conversion logic

    The estimator uses Numpy Array API arrays as its native type.
    """
    np_array_api = pytest.importorskip("numpy.array_api")

    class Engine:
        @staticmethod
        def convert_to_sklearn_types(name, value):
            return np.asarray(value)

    class Estimator:
        _default_engine = DefaultEngine
        # Setup attribute as if `Engine` had previously been selected,
        # we want to test attribute conversion, not engine resolution.
        _engine_class = Engine

        @convert_attributes
        def fit(self, X):
            self.X_ = np_array_api.asarray(X)

    X = np.array([1, 2, 3])
    est = Estimator()
    with config_context(engine_attributes=attribute_types):
        est.fit(X)

    assert isinstance(est.X_, np.ndarray) == converted
    if converted:
        assert est._engine_class == est._default_engine == DefaultEngine
    else:
        assert est._engine_class == Engine


def test_engine_selection():
    """Check that the correct engine is selected."""
    # Values aren't important, just need something to pass as argument
    # to _get_engine
    X = [[1, 2], [3, 4]]

    # If no engine accepts, default engine should be selected
    with config_context(engine_provider=NeverAcceptsEngine):
        est = FakeEstimator()
        engine = est._get_engine(X)
        assert isinstance(engine, DefaultEngine)

    with config_context(engine_provider=AlwaysAcceptsEngine):
        est = FakeEstimator()
        engine = est._get_engine(X)
        assert isinstance(engine, AlwaysAcceptsEngine)

    # Engine with second priority (AlwaysAccepts) is selected
    with config_context(engine_provider=(NeverAcceptsEngine, AlwaysAcceptsEngine)):
        est = FakeEstimator()
        engine = est._get_engine(X)
        assert isinstance(engine, AlwaysAcceptsEngine)


def test_engine_selection_is_fozen():
    """Check that a previously selected engine keeps being used.

    Engine selection is only performed once, after that the same engine
    is used. Re-reselecting the engine is possible when explicitly requested.
    """
    # Values aren't important, just need something to pass as argument
    # to _get_engine
    X = [[1, 2], [3, 4]]

    est = FakeEstimator()

    with config_context(engine_provider=(NeverAcceptsEngine, AlwaysAcceptsEngine)):
        engine = est._get_engine(X)
        assert isinstance(engine, AlwaysAcceptsEngine)

    # Even though `AlsoAlwaysAcceptsEngine` is listed first, it should not
    # be selected
    with config_context(engine_provider=(AlsoAlwaysAcceptsEngine, AlwaysAcceptsEngine)):
        engine = est._get_engine(X)
        assert isinstance(engine, AlwaysAcceptsEngine)

    # Explicitly ask for engine re-selection
    with config_context(engine_provider=(AlsoAlwaysAcceptsEngine, AlwaysAcceptsEngine)):
        engine = est._get_engine(X, reset=True)
        assert isinstance(engine, AlsoAlwaysAcceptsEngine)


def test_missing_engine_raises():
    """Check an exception is raised when a previously configured engine is
    no longer available.
    """
    # Values aren't important, just need something to pass as argument
    # to _get_engine
    X = [[1, 2], [3, 4]]

    est = FakeEstimator()

    with config_context(engine_provider=(NeverAcceptsEngine, AlwaysAcceptsEngine)):
        engine = est._get_engine(X)
        assert isinstance(engine, AlwaysAcceptsEngine)

    # Raise an exception because the previously selected engine isn't available
    with config_context(engine_provider=(AlsoAlwaysAcceptsEngine,)):
        with pytest.raises(RuntimeError, match="Previously selected engine.*"):
            est._get_engine(X)

    # Doesn't raise because `reset=True`
    with config_context(engine_provider=(NeverAcceptsEngine, AlsoAlwaysAcceptsEngine)):
        engine = est._get_engine(X, reset=True)
        assert isinstance(engine, AlsoAlwaysAcceptsEngine)


def test_default_engine_always_works():
    """Check that an estimator that uses the default engine works, even when
    no engines are explicitly configured.
    """
    # Values aren't important, just need something to pass as argument
    # to _get_engine
    X = [[1, 2], [3, 4]]

    est = FakeEstimator()

    with config_context(engine_provider=NeverAcceptsEngine):
        engine = est._get_engine(X)
        assert isinstance(engine, DefaultEngine)

    assert est._engine_class == DefaultEngine

    # With no explicit config, the default engine should still be selected
    engine = est._get_engine(X)
    assert isinstance(engine, DefaultEngine)
