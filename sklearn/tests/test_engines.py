from sklearn._engine import list_engine_provider_names
from sklearn._engine import _parse_entry_point


class FakeEngine:
    pass


class FakeEngineHolder:
    class NestedFakeEngine:
        pass


def test_get_engine_class():
    fake_entry_point = {
        "name": "fake_engine",
        "value": "sklearn.tests.test_engines:FakeEngine"
    }
    spec = _parse_entry_point(fake_entry_point)
    assert spec.name == "fake_engine"
    assert spec.provider_name == "sklearn"  # or should it be scikit-learn?
    assert spec.get_engine_class() is FakeEngine


def test_get_nested_engine_class():
    fake_entry_point = {
        "name": "nested_fake_engine",
        "value": "sklearn.tests.test_engines:FakeEngineHolder.NestedFakeEngine"
    }
    spec = _parse_entry_point(fake_entry_point)
    assert spec.name == "nested_fake_engine"
    assert spec.provider_name == "sklearn"  # or should it be scikit-learn?
    assert spec.get_engine_class() is FakeEngineHolder.NestedFakeEngine


def test_list_engine_provider_names():
    provider_names = list_engine_provider_names()
    for provider_name in provider_names:
        assert isinstance(provider_name, str)


