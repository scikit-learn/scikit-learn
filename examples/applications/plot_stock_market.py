"""
=======================================
Visualizing the stock market structure
=======================================

This example employs several unsupervised learning techniques to extract
the stock market structure from variations in historical quotes.
"""
print __doc__

# Author: Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD

import datetime

import numpy as np
import pylab as pl
from matplotlib import finance
from matplotlib.collections import LineCollection

from sklearn import cluster, covariance, manifold

###############################################################################
# Retrieve the data from Internet

# Choose a time period reasonnably calm (not too long ago so that we get
# high-tech firms, and before the 2008 crash)
d1 = datetime.datetime(2003, 01, 01)
d2 = datetime.datetime(2008, 01, 01)

symbol_dict = {
        'TOT'  : 'Total',
        'XOM'  : 'Exxon',
        'CVX'  : 'Chevron',
        'COP'  : 'ConocoPhillips',
        'VLO'  : 'Valero Energy',
        'MSFT' : 'Microsoft',
        'IBM'  : 'IBM',
        'TWX'  : 'Time Warner',
        'CMCSA': 'Comcast',
        'CVC'  : 'Cablevision',
        'YHOO' : 'Yahoo',
        'DELL' : 'Dell',
        'HPQ'  : 'HP',
        'AMZN' : 'Amazon',
        'TM'   : 'Toyota',
        'CAJ'  : 'Canon',
        'MTU'  : 'Mitsubishi',
        'SNE'  : 'Sony',
        'F'    : 'Ford',
        'HMC'  : 'Honda',
        'NAV'  : 'Navistar',
        'NOC'  : 'Northrop Grumman',
        'BA'   : 'Boeing',
        'KO'   : 'Coca Cola',
        'MMM'  : '3M',
        'MCD'  : 'Mc Donalds',
        'PEP'  : 'Pepsi',
        'KFT'  : 'Kraft Foods',
        'K'    : 'Kellogg',
        'UN'   : 'Unilever',
        'MAR'  : 'Marriott',
        'PG'   : 'Procter Gamble',
        'CL'   : 'Colgate-Palmolive',
        'NWS'  : 'News Corp',
        'GE'   : 'General Electrics',
        'WFC'  : 'Wells Fargo',
        'JPM'  : 'JPMorgan Chase',
        'AIG'  : 'AIG',
        'AXP'  : 'American express',
        'BAC'  : 'Bank of America',
        'GS'   : 'Goldman Sachs',
        'AAPL' : 'Apple',
        'SAP'  : 'SAP',
        'CSCO' : 'Cisco',
        'TXN'  : 'Texas instruments',
        'XRX'  : 'Xerox',
        'LMT'  : 'Lookheed Martin',
        'WMT'  : 'Wal-Mart',
        'WAG'  : 'Walgreen',
        'HD'   : 'Home Depot',
        'GSK'  : 'GlaxoSmithKline',
        'PFE'  : 'Pfizer',
        'SNY'  : 'Sanofi-Aventis',
        'NVS'  : 'Novartis',
        'KMB'  : 'Kimberly-Clark',
        'R'    : 'Ryder',
        'GD'   : 'General Dynamics',
        'RTN'  : 'Raytheon',
        'CVS'  : 'CVS',
        'CAT'  : 'Caterpillar',
        'DD'   : 'DuPont de Nemours',
    }

symbols, names = np.array(symbol_dict.items()).T

quotes = [finance.quotes_historical_yahoo(symbol, d1, d2, asobject=True)
                for symbol in symbols]

open    = np.array([q.open   for q in quotes]).astype(np.float)
close   = np.array([q.close  for q in quotes]).astype(np.float)
# The daily variations of the quotes are what carry most information
variation = close - open

###############################################################################
# Cluster using affinity propagation

correlations = np.corrcoef(variation)
_, labels = cluster.affinity_propagation(correlations)
n_labels = labels.max()

for i in range(n_labels+1):
    print 'Cluster %i: %s' % ((i+1),
                              ', '.join(names[labels==i]))

###############################################################################
# Learn a graphical structure from the correlations
model = covariance.GLassoCV()
# standardize the time series: using correlations rather than covariance
# is more efficient for structure recovery
X = variation.copy().T
X /= X.std(axis=0)
model.fit(X)

###############################################################################
# Find a low-dimension embedding for visualization

# We use a dense eigen_solver to achieve reproducibility (arpack is
# initiated with random vectors that we don't control). In addition, we
# use a large number of neighbors to capture the large-scale structure.
lle = manifold.LocallyLinearEmbedding(out_dim=2, eigen_solver='dense',
            n_neighbors=6)
embedding = lle.fit_transform(X.T).T

###############################################################################
# Visualization
pl.figure(1, facecolor='w', figsize=(10, 8))
pl.clf()
ax = pl.axes([0., 0., 1., 1.])
pl.axis('off')

# Display a graph of the partial correlations
adjacency = model.precision_.copy()
d = 1/np.sqrt(np.diag(adjacency))
adjacency *= d
adjacency *= d[:, np.newaxis]
non_zero = (np.abs(np.triu(adjacency, k=1)) > 0.02)

# Plot the nodes using the coordinnates of our embedding
pl.scatter(embedding[0], embedding[1], s=100*d**2, c=labels,
           cmap=pl.cm.spectral)

# Plot the edges
start_idx, end_idx = np.where(non_zero)
#a sequence of (*line0*, *line1*, *line2*), where::
#
#            linen = (x0, y0), (x1, y1), ... (xm, ym)
segments = [[embedding[:, start], embedding[:, stop]]
            for start, stop in zip(start_idx, end_idx)]
values = np.abs(adjacency[non_zero])
lc = LineCollection(segments,
                     zorder=0, cmap=pl.cm.hot_r,
                     norm=pl.Normalize(0, .7*values.max()),
                     )
lc.set_array(values)
lc.set_linewidths(15*values)
ax.add_collection(lc)

# Add a label to each node. The challenge here is that we want to
# position the labels to avoid overlap with other labels
for index, (name, label, (x, y)) in enumerate(zip(names,
                                              labels, embedding.T)):
    dx = x - embedding[0]
    dx[index] = 1
    dy = y - embedding[1]
    dy[index] = 1
    this_dx = dx[np.argmin(np.abs(dy))]
    this_dy = dy[np.argmin(np.abs(dx))]
    if this_dx > 0:
        horizontalalignment = 'left'
        x = x + .002
    else:
        horizontalalignment = 'right'
        x = x - .002
    if this_dy > 0:
        verticalalignment = 'bottom'
        y = y + .002
    else:
        verticalalignment = 'top'
        y = y - .002
    pl.text(x, y, name, size=10,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            bbox=dict(facecolor='w',
                      edgecolor=pl.cm.spectral(label/float(n_labels)),
                      alpha=.6))

pl.xlim(embedding[0].min() - .15*embedding[0].ptp(),
        embedding[0].max() + .10*embedding[0].ptp(),)
pl.ylim(embedding[1].min() - .03*embedding[1].ptp(),
        embedding[1].max() + .03*embedding[1].ptp())


