from pytest import xfail

from sklearn import config_context

from sklearn.exceptions import FeatureNotCoveredByPluginError


# TODO: document this pytest plugin + write a tutorial on how to develop a new plugin
# and explain good practices regarding testing against sklearn test modules.
def pytest_addoption(parser):
    group = parser.getgroup("Sklearn plugin testing")
    group.addoption(
        "--sklearn-engine-provider",
        action="store",
        nargs=1,
        type=str,
        help="Name of the an engine provider for sklearn to activate for all tests.",
    )


def pytest_runtest_call(item):
    engine_provider = item.config.getoption("sklearn_engine_provider")
    if engine_provider is None:
        return item.runtest()

    with config_context(engine_provider=engine_provider):
        try:
            item.runtest()
        except FeatureNotCoveredByPluginError:
            xfail(
                reason=f"This test cover features that are not supported by the "
                f"engine provided by {engine_provider}."
            )
