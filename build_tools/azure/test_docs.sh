#!/bin/bash

set -e

if [[ "$DISTRIB" == "conda" ]]; then
    source activate $VIRTUALENV
elif [[ "$DISTRIB" == "ubuntu" ]]; then
    source $VIRTUALENV/bin/activate
elif [[ "$DISTRIB" == "32bit" ]]; then
    # make test-doc will currently not do anything, but for consistency, let us
    # activate the proper environment and call it in the event we may want to
    # run it in the future.
    source $VIRTUALENV/bin/activate
fi

make test-doc
