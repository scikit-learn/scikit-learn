#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]] || [[ "$DISTRIB" == "scipy-dev" ]]; then
    export PATH=$HOME/miniconda3/bin:$PATH
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
fi

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
python -c "\
try:
    import pandas
    print('pandas %s' % pandas.__version__)
except ImportError:
    pass
"

mkdir -p $TEST_DIR
cp setup.cfg $TEST_DIR
cd $TEST_DIR

export SKLEARN_SKIP_NETWORK_TESTS=1

pytest --showlocals --durations=20 --pyargs sklearn/preprocessing/tests/
