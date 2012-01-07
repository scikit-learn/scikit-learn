# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 15:47:17 2012

@author: du
"""

"""
This script shows how to use Gaussian HMM.
It uses stock price data, which can be obtained from yahoo finance.
For more information on how to get stock price by matplitlib, please refere
to date_demo1.py of matplotlib
"""

import sys
import numpy as np
from pylab import figure, show
from matplotlib.finance import quotes_historical_yahoo
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import datetime
import sklearn.hmm

def getQuotes(tick):
    
    # start date
    date1 = datetime.date( 1995, 1, 1 )
    # end date    
    date2 = datetime.date( 2012, 1, 6 )

    # get quotes from yahoo finance    
    quotes = quotes_historical_yahoo(tick, date1, date2)
    if len(quotes) == 0:
        raise SystemExit
    
    # quotes[t] = (date, open_v, close_v, highest_v, lowest_v, volume)
    return quotes
    
def testHMM(n_components, quotes):
    # unpack quotes 
    dates = np.array([q[0] for q in quotes], dtype=int)
    close_v = np.array([q[2] for q in quotes])
    volume = np.array([q[2] for q in quotes])[1:] 

    # take diff of close value
    # this makes len(diff) = len(close_t) - 1
    # therefore, others quontity also need to be shifted
    diff = close_v[1:] - close_v[:-1]
    dates = dates[1:]
    close_v = close_v[1:]
    
    # pack diff and volume for training
    X = np.column_stack([diff, volume])
    
    # make an HMM instance and execute fit
    model = sklearn.hmm.GaussianHMM(n_components, "diag")
    model.fit([X], n_iter=1000)

    # predict the optimal sequence of internal hidden state
    hidden_states = model.decode(X)[1]
    
    return dates, close_v, volume, hidden_states, model


tick = "INTC"
n_components = 5

print "downloading quotes data"
sys.stdout.flush()
quotes = getQuotes(tick)

print "fittiing to HMM and decoding ...",
sys.stdout.flush()
dates, close_v, volume, hidden_states, model = testHMM(n_components, quotes)

print "done\n"
sys.stdout.flush()

# print trained parameters

print "transition matrix"
print model.transmat
print ""

print "means and vars of each hidden state"
for i in xrange(n_components):
    print "%dth hidden state"%i
    print "mean = ", model.means[i]
    print "var = ", np.diag(model.covars[i])
    print ""

years    = YearLocator()   # every year
months   = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y')
fig = figure()
ax = fig.add_subplot(111)

for i in xrange(n_components):
    # use fancy indexing to plot data in each state
    idx = (hidden_states == i)
    ax.plot_date(dates[idx], close_v[idx], 'o', label="%dth hidden state"%i)
ax.legend()

# format the ticks
ax.xaxis.set_major_locator(years)
ax.xaxis.set_major_formatter(yearsFmt)
ax.xaxis.set_minor_locator(months)
ax.autoscale_view()

#format the coords message box
def price(x): return '$%1.2f'%x
ax.fmt_xdata = DateFormatter('%Y-%m-%d')
ax.fmt_ydata = price
ax.grid(True)

fig.autofmt_xdate()
show()
