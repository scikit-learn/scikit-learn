from __future__ import absolute_import, division, print_function

import pytest


def pytest_addoption(parser):
    group = parser.getgroup("debugconfig")
    group.addoption('--setupplan', '--setup-plan', action="store_true",
                    help="show what fixtures and tests would be executed but "
                    "don't execute anything.")


@pytest.hookimpl(tryfirst=True)
def pytest_fixture_setup(fixturedef, request):
    # Will return a dummy fixture if the setuponly option is provided.
    if request.config.option.setupplan:
        fixturedef.cached_result = (None, None, None)
        return fixturedef.cached_result


@pytest.hookimpl(tryfirst=True)
def pytest_cmdline_main(config):
    if config.option.setupplan:
        config.option.setuponly = True
        config.option.setupshow = True
