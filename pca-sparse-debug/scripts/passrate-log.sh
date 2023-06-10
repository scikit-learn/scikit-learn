#!/bin/bash

# writes pca sparse pass/fail results to stdout

SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all" OMP_NUM_THREADS=1 pytest --color=no -v --tb=no -n "$(nproc --all)" sklearn/decomposition/tests/test_pca.py::test_pca_sparse
