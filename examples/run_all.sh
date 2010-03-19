#!/bin/bash
#
# small scripts that runs all examples for testing purposes

for file in `ls *.py`; do
    python $file
done

