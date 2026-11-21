# sklearn/mixture/tests/conftest.py

# default values
COVARIANCE_TYPE = ["full", "tied", "diag", "tied-diag", "spherical", "tied-spherical"]
PRIOR_TYPE = ["dirichlet_process", "dirichlet_distribution"]


def pytest_addoption(parser):
    parser.addoption(
        "--covariance-types",
        action="store",
        default=",".join(COVARIANCE_TYPE),
        help="Comma-separated covariance types (e.g., full,diag,spherical)",
    )


def pytest_configure(config):
    # overwrite the global variable dynamically
    global COVARIANCE_TYPE
    COVARIANCE_TYPE = [
        t.strip() for t in config.getoption("--covariance-types").split(",")
    ]
