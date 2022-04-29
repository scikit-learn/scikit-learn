#!/bin/bash

set -e
set -x

echo $EVENT_OBJ

echo $GITHUB_REF | sed -e "s/^refs\/\(heads\/\)\{0,1\}//" -e "s/\/merge$/\/head/" > branch.txt
