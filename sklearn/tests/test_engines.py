import re
from collections import namedtuple
import pytest

from sklearn._engine import list_engine_provider_names
from sklearn._engine import _parse_entry_point
from sklearn._engine import get_engine_class
from sklearn._engine import _get_engine_class
from sklearn._engine import EngineSpec
from sklearn._config import config_context


class FakeEngine:
    pass


class FakeEngineHolder:
    class NestedFakeEngine:
        pass


FakeEntryPoint = namedtuple("FakeEntryPoint", ["name", "value"])


def test_get_engine_class():
    fake_entry_point = FakeEntryPoint(
        name="fake_engine",
        value="sklearn.tests.test_engines:FakeEngine",
    )
    spec = _parse_entry_point(fake_entry_point)
    assert spec.name == "fake_engine"
    assert spec.provider_name == "sklearn"  # or should it be scikit-learn?
    assert spec.get_engine_class() is FakeEngine


def test_get_nested_engine_class():
    fake_entry_point = FakeEntryPoint(
        name="nested_fake_engine",
        value="sklearn.tests.test_engines:FakeEngineHolder.NestedFakeEngine",
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
        engine_class = get_engine_class("test_missing_engine_name", default=FakeEngine)
    assert engine_class is FakeEngine


def test_get_engine_class_for_invalid_provider():
    expected_message = re.escape(
        "Could not find any provider for the sklearn_engines entry point with"
        " name(s): 'invalid_provider_name'"
    )
    with pytest.raises(RuntimeError, match=expected_message):
        with config_context(engine_provider="invalid_provider_name"):
            get_engine_class("kmeans")

    expected_message = re.escape(
        "Could not find any provider for the sklearn_engines entry point with"
        " name(s): 'invalid_provider_name_1', 'invalid_provider_name_2'"
    )
    with pytest.raises(RuntimeError, match=expected_message):
        with config_context(
            engine_provider=("invalid_provider_name_1", "invalid_provider_name_2")
        ):
            get_engine_class("kmeans")


def test_get_engine_class():
    engine_specs = (
        EngineSpec("kmeans", "provider1", "sklearn.provider1.module", "KMeansEngine"),
        EngineSpec("other", "provider1", "sklearn.provider1.module", "OtherEngine"),
        EngineSpec("kmeans", "provider2", "sklearn.provider2.module", "KMeansEngine"),
        EngineSpec("kmeans", "provider3", "sklearn.tests.test_engines", "FakeEngine"),
        EngineSpec(
            "kmeans",
            "provider4",
            "sklearn.tests.test_engines",
            "FakeEngineHolder.NestedFakeEngine",
        ),
    )

    engine_class = _get_engine_class(
        engine_name="missing",
        provider_names=("provider1", "provider3"),
        engine_specs=engine_specs,
    )
    assert engine_class is None

    engine_class = _get_engine_class(
        engine_name="kmeans",
        provider_names=("provider3", "provider4", "provider1", "provider2"),
        engine_specs=engine_specs,
    )
    assert engine_class == FakeEngine

    engine_class = _get_engine_class(
        engine_name="kmeans",
        provider_names=("provider4", "provider3", "provider1", "provider2"),
        engine_specs=engine_specs,
    )
    assert engine_class == FakeEngineHolder.NestedFakeEngine

    with pytest.raises(ImportError, match=re.escape("sklearn.provider1")):
        # Invalid imports are delayed until they are actually needed.
        _get_engine_class(
            engine_name="kmeans",
            provider_names=("provider1", "provider3"),
            engine_specs=engine_specs,
        )
