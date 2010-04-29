#!/bin/bash
#
# small scripts that runs all examples for testing purposes

for file in `ls *.py`; do
    python $file
done

for file in `ls svm/*.py`; do
    python $file
done

for file in `ls em/*.py`; do
    python $file
done
