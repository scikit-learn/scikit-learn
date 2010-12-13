"""
===========================================
Finding structure in the stock market
===========================================

An example of playing with stock market data to try and find some
structure in it.
"""
print __doc__

# Author: Gael Varoquaux gael.varoquaux@normalesup.org
# License: BSD

import datetime
from matplotlib import finance
import numpy as np

from scikits.learn import cluster

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
        'HPQ'  : 'Hewlett-Packard',
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
        'NWS'  : 'News Corporation',
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

#volumes = np.array([q.volume for q in quotes]).astype(np.float)
open    = np.array([q.open   for q in quotes]).astype(np.float)
close   = np.array([q.close  for q in quotes]).astype(np.float)
variation = close - open
correlations = np.corrcoef(variation)

_, labels = cluster.affinity_propagation(correlations)

for i in range(labels.max()+1):
    print 'Cluster %i: %s' % ((i+1),
                              ', '.join(names[labels==i]))
