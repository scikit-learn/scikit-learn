#!/bin/bash

set -e
set -x

PYTEST="pytest -n $CI_CPU_COUNT" make test-doc
