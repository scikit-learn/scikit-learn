from sklearn._engines import list_engine_provider_names
from sklearn._engines import load_engines
import pytest


def test_list_engine_provider_names():
    provider_names = list_engine_provider_names()
    for provider_name in provider_names:
        assert isinstance(provider_name, str)


def test_load_engines():
    all_engines = load_engines()
    # TODO: write me

def test_load_engines_invalid_provider():
    provider_name = "some_invalid_test_provider_name"
    assert provider_name not in list_engine_provider_names()
    assert load_engines(provider_name=provider_name) == []
