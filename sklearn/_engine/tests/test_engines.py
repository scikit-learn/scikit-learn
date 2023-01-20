import re
from collections import namedtuple
import pytest

import numpy as np

from sklearn._engine import convert_attributes
from sklearn._engine import list_engine_provider_names
from sklearn._engine import get_engine_classes
from sklearn._engine.base import _parse_entry_point
from sklearn._engine.base import _get_engine_classes
from sklearn._engine.base import EngineSpec
from sklearn._config import config_context
from sklearn.base import EngineAwareMixin


class FakeDefaultEngine:
    pass


class FakeEngine:
    pass


class FakeEngineHolder:
    class NestedFakeEngine:
        pass


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
        _engine_class = Engine

        @convert_attributes
        def fit(self, X):
            self.X_ = np_array_api.asarray(X)

    X = np.array([1, 2, 3])
    est = Estimator()
    with config_context(engine_attributes=attribute_types):
        est.fit(X)

    assert isinstance(est.X_, np.ndarray) == converted


def test_engine_selection():
    """Check that the correct engine is selected"""

    class DefaultEngine:
        def __init__(self, estimator):
            self.estimator = estimator

        def accepts(self, X, y=None, sample_weight=None):
            return True

    class NeverAcceptsEngine:
        def __init__(self, estimator):
            self.estimator = estimator

        def accepts(self, X, y=None, sample_weight=None):
            return False

    class AlwaysAcceptsEngine:
        def __init__(self, estimator):
            self.estimator = estimator

        def accepts(self, X, y=None, sample_weight=None):
            return True

    class Estimator(EngineAwareMixin):
        _engine_name = "test-engine"
        _default_engine = DefaultEngine

        def resolve_engine(self, X, y):
            # Test that `_get_engine` works properly
            engine = self._get_engine(X, y, reset=True)
            return engine

    with config_context(engine_provider=(NeverAcceptsEngine)):
        est = Estimator()
        engine = est.resolve_engine([[1, 2], [3, 4]], [0, 1])
        assert isinstance(engine, DefaultEngine)

    with config_context(engine_provider=(AlwaysAcceptsEngine)):
        est = Estimator()
        engine = est.resolve_engine([[1, 2], [3, 4]], [0, 1])
        assert isinstance(engine, AlwaysAcceptsEngine)

    with config_context(engine_provider=(NeverAcceptsEngine, AlwaysAcceptsEngine)):
        est = Estimator()
        engine = est.resolve_engine([[1, 2], [3, 4]], [0, 1])
        assert isinstance(engine, AlwaysAcceptsEngine)
