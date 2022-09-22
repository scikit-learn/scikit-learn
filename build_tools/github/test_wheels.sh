#!/bin/bash

set -e
set -x

if [[ "$OSTYPE" != "linux-gnu" ]]; then
    # The Linux test environment is run in a Docker container and
    # it is not possible to copy the test configuration file (yet)
    cp $CONFTEST_PATH $CONFTEST_NAME
fi

# Test that there are no links to system libraries in the
# threadpoolctl output section of the show_versions output:
python -c "import sklearn; sklearn.show_versions()"
pytest --pyargs sklearn
