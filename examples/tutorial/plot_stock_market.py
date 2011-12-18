import datetime
from matplotlib import finance
import numpy as np

################################################################################
if 1:
    # Choose a time period reasonnably calm: before the 2008 crash)
    # and after Google's start
    d1 = datetime.datetime(2004, 8, 19)
    d2 = datetime.datetime(2008, 01, 01)

    symbol_dict = {
            'TOT'  : 'Total',
            'XOM'  : 'Exxon',
            'CVX'  : 'Chevron',
            'COP'  : 'ConocoPhillips',
            'VLO'  : 'Valero Energy',
            'MSFT' : 'Microsoft',
            'B'    : 'Barnes group',
            'EK'   : 'Eastman kodak',
            'IBM'  : 'IBM',
            'TWX'  : 'Time Warner',
            'CMCSA': 'Comcast',
            'CVC'  : 'Cablevision',
            'YHOO' : 'Yahoo',
            'DELL' : 'Dell',
            'DE'   : 'Deere and company',
            'HPQ'  : 'Hewlett-Packard',
            'AMZN' : 'Amazon',
            'TM'   : 'Toyota',
            'CAJ'  : 'Canon',
            'MTU'  : 'Mitsubishi',
            'SNE'  : 'Sony',
            'EMR'  : 'Emerson electric',
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
            'VOD'  : 'Vodaphone',
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
            'PMI'  : 'PMI group',
            'AAPL' : 'Apple',
            'SAP'  : 'SAP',
            'CSCO' : 'Cisco',
            'QCOM' : 'Qualcomm',
            'HAL'  : 'Haliburton',
            'HTCH' : 'Hutchinson',
            'JDSU' : 'JDS uniphase',
            'TXN'  : 'Texas instruments',
            'O'    : 'Reality income',
            'UPS'  : 'UPS',
            'BP'   : 'BP',
            'L'    : 'Loews corporation',
            'M'    : "Macy's",
            'S'    : 'Sprint nextel',
            'XRX'  : 'Xerox',
            'WYNN' : 'Wynn resorts',
            'DIS'  : 'Walt disney',
            'WFR'  : 'MEMC electronic materials',
            'UTX'  : 'United Technology corp',
            'X'    : 'United States Steel corp',
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
            'MON'  : 'Monsanto',
            'CLF'  : 'Cliffs natural ressources',
            'BTU'  : 'Peabody energy',
            'ACI'  : 'Arch Coal',
            'BTU'  : 'Patriot coal corp',
            'PPG'  : 'PPC',
            'CMI'  : 'Cummins common stock',
            'JNJ'  : 'Johnson and johnson',
            'ABT'  : 'Abbott laboratories',
            'MRK'  : 'Merck and co',
            'T'    : 'AT and T',
            'VZ'   : 'Verizon',
            'FTR'  : 'Frontiers communication',
            'CTL'  : 'Centurylink',
            'MO'   : 'Altria group',
            'NLY'  : 'Annaly capital management',
            'QQQ'  : 'Powershares',
            'BMY'  : 'Bristal-myers squibb',
            'LLY'  : 'Eli lilly and co',
            'C'    : 'Citigroup',
            'MS'   : 'Morgan Stanley',
            'TGT'  : 'Target corporation',
            'SCHW' : 'Charles schwad',
            'ETFC' : 'E*Trade',
            'AMTD' : 'TD ameritrade holding',
            'INTC' : 'Intel',
            'AMD'  : 'AMD',
            'NOK'  : 'Nokia',
            'MU'   : 'Micron technologies',
            'NVDA' : 'Nvidia',
            'MRVL' : 'Marvel technology group',
            'SNDK' : 'Sandisk',
            'RIMM' : 'Research in mention',
            'TXN'  : 'Texas instruments',
            'EMC'  : 'EMC',
            'ORCL' : 'Oracle',
            'LOW'  : "Lowe's",
            'BBY'  : 'Best buy',
            'FDX'  : 'Fedex',
            'FE'   : 'First energy',
            'JNPR' : 'Juniper',
            'GOOG' : 'Google',
            'AXP'  : 'American express',
            'AMAT' : 'Applied material',
            '^DJI' : 'Dow Jones Industrial average',
            '^DJA' : 'Dow Jones Composite average',
            '^DJT' : 'Dow Jones Transportation average',
            '^DJU' : 'Dow Jones Utility average',
            '^IXIC': 'Nasdaq composite',
            #'^FCHI': 'CAC40',
        }

    symbols, names = np.array(symbol_dict.items()).T

    quotes = [finance.quotes_historical_yahoo(symbol, d1, d2, asobject=True)
                    for symbol in symbols]

    #volumes = np.array([q.volume for q in quotes]).astype(np.float)
    open    = np.array([q.open   for q in quotes]).astype(np.float)
    close   = np.array([q.close  for q in quotes]).astype(np.float)
    variation = close - open
    np.save('variation.npy', variation)
    np.save('names.npy', names)
else:
    names = np.load('names.npy')
    variation = np.load('variation.npy')

################################################################################
# Get our X and y variables
X = variation[names != 'Google'].T
y = variation[names == 'Google'].squeeze()
n = names[names != 'Google']

X -= X.mean(axis=0)
X /= X.std(axis=0)

################################################################################
# Shuffle the data because of non-stationarity
np.random.seed(0)
order = np.random.permutation(len(X))
X = X[order]
y = y[order]

from scikits.learn import linear_model, cross_val

lasso = linear_model.LassoCV()
scores = cross_val.cross_val_score(lasso, X, y, n_jobs=-1)
lasso.fit(X, y)
print n[lasso.coef_ != 0]

################################################################################
# Plot the time series
import pylab as pl
pl.plot(variation[names == 'Google'].squeeze())
pl.title('Google stock value')
pl.xlabel('Time')
pl.ylabel('Google quote daily increment')

#_, labels, _ = cluster.k_means(X.T, 10)
#
#for i in range(labels.max()+1):
#    print 'Cluster %i: %s' % ((i+1),
#                              ', '.join(n[labels==i]))
