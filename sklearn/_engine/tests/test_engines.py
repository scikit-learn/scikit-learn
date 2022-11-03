import re
from collections import namedtuple
import pytest

from sklearn._engine import list_engine_provider_names
from sklearn._engine import get_engine_classes
from sklearn._engine.base import _parse_entry_point
from sklearn._engine.base import _get_engine_classes
from sklearn._engine.base import EngineSpec
from sklearn._config import config_context


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
    assert engine_classes == [FakeEngine]


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
    assert engine_class == [FakeDefaultEngine]

    engine_class = list(
        _get_engine_classes(
            engine_name="kmeans",
            provider_names=("provider3", "provider4"),
            engine_specs=engine_specs,
            default=FakeDefaultEngine,
        )
    )
    assert engine_class == [
        FakeEngine,
        FakeEngineHolder.NestedFakeEngine,
        FakeDefaultEngine,
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
        FakeEngineHolder.NestedFakeEngine,
        FakeEngine,
        FakeDefaultEngine,
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
