#!/bin/bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"

# Skip tests that require large downloads over the network to save bandwith
# usage as travis workers are stateless and therefore traditional local
# disk caching does not work.
export SKLEARN_SKIP_NETWORK_TESTS=1

# Do not use "make test" or "make test-coverage" as they enable verbose mode
# which renders travis output too slow to display in a browser.
if [[ "$COVERAGE" == "true" ]]; then
    nosetests -s --with-coverage sklearn
else
    nosetests -s sklearn
fi

make test-doc test-sphinxext
