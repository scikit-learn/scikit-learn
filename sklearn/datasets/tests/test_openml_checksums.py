"""Test the openml loader with MD5 checksums.

Clearly this is only a proof of concept, and must be replaced
with the, somewhat convoluted, existing "monkey-mocking-patching"
testsuite...
"""

from sklearn.datasets import fetch_openml


def test():
    fetch_openml('iris', version=1, verify_checksum=True)
