from pytest import hookimpl

from sklearn import config_context


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


@hookimpl(hookwrapper=True)
def pytest_pyfunc_call(pyfuncitem):
    engine_provider = pyfuncitem.config.getoption("sklearn_engine_provider")
    if engine_provider is None:
        yield
        return

    with config_context(engine_provider=engine_provider):
        outcome = yield
        outcome.get_result()
