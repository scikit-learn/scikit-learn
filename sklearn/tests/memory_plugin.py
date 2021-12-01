import os
import pytest

try:
    import psutil
except ImportError:
    psutil = None


def pytest_addoption(parser):
    parser.addoption(
        "--min-ram",
        type=int,
        help="Minimal amount of memory (MB) incremented to display in summary",
    )


def pytest_configure(config):
    min_ram = config.getoption("--min-ram")
    if min_ram is not None:
        if psutil is None:
            raise ValueError("psutil needs to be installed to use --min-ram")
        memory_usage_plugin = MemoryUsage(config, min_ram)
        config.pluginmanager.register(memory_usage_plugin)


def is_main_process(config):
    """True if the code running the given pytest.config object is running in a xdist main
    node or not running xdist at all.
    """
    return not hasattr(config, "workerinput")


def get_usage_memory_mb():
    """
    Measures memory usage per Python process

    Returns
    -------
    memory_usage : float
    """
    process = psutil.Process(os.getpid())
    memory_use = process.memory_info()
    return memory_use.rss / 1e6  # to MB


class MemoryUsage:
    """Measure memory usage during test execution.

    This plugin can also collect memory usage data from pytest-xdist workers.

    This code is directly adapted from:
    https://gist.github.com/DKorytkin/8a186693af9a015abe89f6b874ca0795 by
    Dmytro Korytkin.
    """

    SHARED_MEMORY_USAGE_INFO = "memory_usage"

    def __init__(self, config, min_ram):
        """Setup the plugin.

        Parameters
        ----------
        config : _pytest.config.Config
        min_ram : int
        """
        self.config = config
        self.min_ram = min_ram
        self.is_main_process = is_main_process(config)
        self.stats = {}

    def add(self, name):
        self.stats[name] = self.stats.get(name) or {}
        return self.stats[name]

    def pytest_runtest_setup(self, item):
        """Record maxrss for pre-setup."""
        self.add(item.nodeid)["setup"] = get_usage_memory_mb()

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_call(self, item):
        """Track test memory usage

        This collects memory usage data before and after each test run.

        However, we do not measure the peak memory usage during the test run
        itself. This would require some kind of sampling mechanism (e.g. using
        `memory_profiler`) but this does not see to be easily achievable using
        the pytest plugin API.

        Parameters
        ----------
        item : _pytest.main.Item
        """
        start = get_usage_memory_mb()
        yield
        end = get_usage_memory_mb()
        node_stats = self.add(item.nodeid)
        if "setup" in node_stats:
            reference = node_stats["setup"]
        else:
            reference = start

        node_stats["diff"] = end - reference
        node_stats["end"] = end
        node_stats["start"] = start

    def pytest_terminal_summary(self, terminalreporter):
        tr = terminalreporter
        if self.stats:
            stats_filtered = filter(
                lambda x: x[-1]["diff"] > self.min_ram, self.stats.items()
            )
            stats = sorted(stats_filtered, key=lambda x: x[-1]["diff"], reverse=True)
            tr._tw.sep(
                "=",
                f"{len(stats)} tests that incremented RAM more than {self.min_ram} MB",
                yellow=True,
            )
            for test_name, info in stats:
                line = "before setup: {:.3f} MB, increment: {:.3f} MB - {}".format(
                    info["setup"], info["diff"], test_name
                )
                tr._tw.line(line)

    def pytest_testnodedown(self, node, error):
        """Collect memory usage stats from pytest-xdist nodes.

        All stats are merged into the main process stats.
        """
        node_stats = node.workeroutput[self.SHARED_MEMORY_USAGE_INFO]
        self.stats.update(node_stats)

    @pytest.hookimpl(hookwrapper=True, trylast=True)
    def pytest_sessionfinish(self, session, exitstatus):
        """Dump memory usage statistics to `workeroutput`

        Executed once per node if with xdist and will gen from mater node

        Parameters
        ----------
        session : _pytest.Session
        exitstatus : int
        """
        yield
        if not self.is_main_process:
            self.config.workeroutput[self.SHARED_MEMORY_USAGE_INFO] = self.stats
