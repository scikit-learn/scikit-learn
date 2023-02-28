set -euxo pipefail

SCRIPTDIR=$(dirname "$0")
DATADIR=$SCRIPTDIR/../data
PLOTDIR=$SCRIPTDIR/../plots
TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
GITHASH=$(git rev-parse --short HEAD)

BASENAME=pca-sparse-passrate-$GITHASH-$TIMESTAMP
LOGFILE=$DATADIR/$BASENAME.log
CSVFILE=$DATADIR/$BASENAME.csv
PLOTFILE=$PLOTDIR/$BASENAME.png

mkdir -p $DATADIR
mkdir -p $PLOTDIR

bash $SCRIPTDIR/passrate-log.sh > $LOGFILE || true
bash $SCRIPTDIR/passrate-csv.sh $LOGFILE > $CSVFILE
python $SCRIPTDIR/passrate-plot.py $CSVFILE $PLOTFILE
