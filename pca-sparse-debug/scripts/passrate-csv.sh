#!/bin/bash

LOGFILE=$1
echo 'seed,rtol,solver,format,k,density,outcome'
grep -P 'PASSED|FAILED' $LOGFILE | sed -E -e 's/^.*(FAILED|PASSED).*\[(.*)\]/\2 \1/' -e 's/1e-/1e@/g' -e 's/-/ /g' -e 's/@/-/g' -e 's/ $//' -e 's/ /,/g'
