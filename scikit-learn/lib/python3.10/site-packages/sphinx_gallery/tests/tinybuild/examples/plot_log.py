"""
Logging test
============
"""

# %%

import logging
import sys


def _create_logger(name=None):
    logger = logging.getLogger(name=name)
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)
    return logger


# %%
logger = _create_logger("first_logger")
logger.info("is in the same cell")

# %%
logger.info("is not in the same cell")
