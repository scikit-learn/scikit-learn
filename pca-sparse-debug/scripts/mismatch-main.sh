set -euxo pipefail

SCRIPTDIR=$(dirname "$0")
DATADIR=$SCRIPTDIR/../data
PLOTDIR=$SCRIPTDIR/../plots
TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
GITHASH=$(git rev-parse --short HEAD)

BASENAME=pca-sparse-mismatch-$GITHASH-$TIMESTAMP
LOGFILE=$DATADIR/$BASENAME.log
CSVFILE=$DATADIR/$BASENAME.csv
PLOTFILE=$PLOTDIR/$BASENAME.png

mkdir -p $DATADIR
mkdir -p $PLOTDIR

bash $SCRIPTDIR/mismatch-log.sh > $LOGFILE || true
bash $SCRIPTDIR/mismatch-csv.sh $LOGFILE > $CSVFILE
python $SCRIPTDIR/mismatch-plot.py $CSVFILE $PLOTFILE
