LOGFILE=$1
paste -d ',' \
      <(echo 'seed,solver,format,k,density'; grep -Po '(?<=_ test_pca_sparse\[).+?(?=\])' $LOGFILE | sed 's/-/,/g') \
      <(echo 'bad,total'; grep -Po '(?<=Mismatched elements: )\d+ / \d+' $LOGFILE | sed -e 's/ //g' -e 's/\//,/g')
