"""Fixture module to skip the datasets loading when offline

The RCV1 data is rather large and some CI workers such as travis are
stateless hence will not cache the dataset as regular scikit-learn users would do.

The following will skip the execution of the rcv1.rst doctests
if the proper environment variable is configured (see the source code of
check_skip_network for more details).

"""
from sklearn.utils.testing import check_skip_network, SkipTest
import os
from sklearn.datasets import get_data_home


def setup_module():
    check_skip_network()

    # skip the test in rcv1.rst if the dataset is not already loaded
    rcv1_dir = os.path.join(get_data_home(), "RCV1")
    if not os.path.exists(rcv1_dir):
        raise SkipTest("Download RCV1 dataset to run this test.")

