#!/bin/bash

LOGFILE=$1
paste -d ',' \
      <(echo 'seed,rtol,solver,layout,k,density'; grep -Po '(?<=test_pca_sparse\[).+?(?=\])' $LOGFILE | sed -e 's/1e-/1e@/g' -e 's/-/,/g' -e 's/@/-/g') \
      <(echo 'bad,total'; grep -Po '(?<=Mismatched elements: )\d+ / \d+' $LOGFILE | sed -e 's/ //g' -e 's/\//,/g')
