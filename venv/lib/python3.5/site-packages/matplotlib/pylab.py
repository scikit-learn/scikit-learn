"""
This is a procedural interface to the matplotlib object-oriented
plotting library.

The following plotting commands are provided; the majority have
MATLAB |reg| [*]_ analogs and similar arguments.

.. |reg| unicode:: 0xAE

_Plotting commands
  acorr     - plot the autocorrelation function
  annotate  - annotate something in the figure
  arrow     - add an arrow to the axes
  axes      - Create a new axes
  axhline   - draw a horizontal line across axes
  axvline   - draw a vertical line across axes
  axhspan   - draw a horizontal bar across axes
  axvspan   - draw a vertical bar across axes
  axis      - Set or return the current axis limits
  autoscale - turn axis autoscaling on or off, and apply it
  bar       - make a bar chart
  barh      - a horizontal bar chart
  broken_barh - a set of horizontal bars with gaps
  box       - set the axes frame on/off state
  boxplot   - make a box and whisker plot
  violinplot - make a violin plot
  cla       - clear current axes
  clabel    - label a contour plot
  clf       - clear a figure window
  clim      - adjust the color limits of the current image
  close     - close a figure window
  colorbar  - add a colorbar to the current figure
  cohere    - make a plot of coherence
  contour   - make a contour plot
  contourf  - make a filled contour plot
  csd       - make a plot of cross spectral density
  delaxes   - delete an axes from the current figure
  draw      - Force a redraw of the current figure
  errorbar  - make an errorbar graph
  figlegend - make legend on the figure rather than the axes
  figimage  - make a figure image
  figtext   - add text in figure coords
  figure   - create or change active figure
  fill     - make filled polygons
  findobj  - recursively find all objects matching some criteria
  gca      - return the current axes
  gcf      - return the current figure
  gci      - get the current image, or None
  getp      - get a graphics property
  grid     - set whether gridding is on
  hist     - make a histogram
  ioff     - turn interaction mode off
  ion      - turn interaction mode on
  isinteractive - return True if interaction mode is on
  imread   - load image file into array
  imsave   - save array as an image file
  imshow   - plot image data
  legend   - make an axes legend
  locator_params - adjust parameters used in locating axis ticks
  loglog   - a log log plot
  matshow  - display a matrix in a new figure preserving aspect
  margins  - set margins used in autoscaling
  pause    - pause for a specified interval
  pcolor   - make a pseudocolor plot
  pcolormesh - make a pseudocolor plot using a quadrilateral mesh
  pie      - make a pie chart
  plot     - make a line plot
  plot_date - plot dates
  plotfile  - plot column data from an ASCII tab/space/comma delimited file
  pie      - pie charts
  polar    - make a polar plot on a PolarAxes
  psd      - make a plot of power spectral density
  quiver   - make a direction field (arrows) plot
  rc       - control the default params
  rgrids   - customize the radial grids and labels for polar
  savefig  - save the current figure
  scatter  - make a scatter plot
  setp      - set a graphics property
  semilogx - log x axis
  semilogy - log y axis
  show     - show the figures
  specgram - a spectrogram plot
  spy      - plot sparsity pattern using markers or image
  stem     - make a stem plot
  subplot  - make one subplot (numrows, numcols, axesnum)
  subplots - make a figure with a set of (numrows, numcols) subplots
  subplots_adjust - change the params controlling the subplot positions of current figure
  subplot_tool - launch the subplot configuration tool
  suptitle   - add a figure title
  table    - add a table to the plot
  text     - add some text at location x,y to the current axes
  thetagrids - customize the radial theta grids and labels for polar
  tick_params - control the appearance of ticks and tick labels
  ticklabel_format - control the format of tick labels
  title    - add a title to the current axes
  tricontour - make a contour plot on a triangular grid
  tricontourf - make a filled contour plot on a triangular grid
  tripcolor - make a pseudocolor plot on a triangular grid
  triplot - plot a triangular grid
  xcorr   - plot the autocorrelation function of x and y
  xlim     - set/get the xlimits
  ylim     - set/get the ylimits
  xticks   - set/get the xticks
  yticks   - set/get the yticks
  xlabel   - add an xlabel to the current axes
  ylabel   - add a ylabel to the current axes

  autumn - set the default colormap to autumn
  bone   - set the default colormap to bone
  cool   - set the default colormap to cool
  copper - set the default colormap to copper
  flag   - set the default colormap to flag
  gray   - set the default colormap to gray
  hot    - set the default colormap to hot
  hsv    - set the default colormap to hsv
  jet    - set the default colormap to jet
  pink   - set the default colormap to pink
  prism  - set the default colormap to prism
  spring - set the default colormap to spring
  summer - set the default colormap to summer
  winter - set the default colormap to winter

_Event handling

  connect - register an event handler
  disconnect - remove a connected event handler

_Matrix commands

  cumprod   - the cumulative product along a dimension
  cumsum    - the cumulative sum along a dimension
  detrend   - remove the mean or besdt fit line from an array
  diag      - the k-th diagonal of matrix
  diff      - the n-th differnce of an array
  eig       - the eigenvalues and eigen vectors of v
  eye       - a matrix where the k-th diagonal is ones, else zero
  find      - return the indices where a condition is nonzero
  fliplr    - flip the rows of a matrix up/down
  flipud    - flip the columns of a matrix left/right
  linspace  - a linear spaced vector of N values from min to max inclusive
  logspace  - a log spaced vector of N values from min to max inclusive
  meshgrid  - repeat x and y to make regular matrices
  ones      - an array of ones
  rand      - an array from the uniform distribution [0,1]
  randn     - an array from the normal distribution
  rot90     - rotate matrix k*90 degress counterclockwise
  squeeze   - squeeze an array removing any dimensions of length 1
  tri       - a triangular matrix
  tril      - a lower triangular matrix
  triu      - an upper triangular matrix
  vander    - the Vandermonde matrix of vector x
  svd       - singular value decomposition
  zeros     - a matrix of zeros

_Probability

  normpdf   - The Gaussian probability density function
  rand      - random numbers from the uniform distribution
  randn     - random numbers from the normal distribution

_Statistics

  amax      - the maximum along dimension m
  amin      - the minimum along dimension m
  corrcoef  - correlation coefficient
  cov       - covariance matrix
  mean      - the mean along dimension m
  median    - the median along dimension m
  norm      - the norm of vector x
  prod      - the product along dimension m
  ptp       - the max-min along dimension m
  std       - the standard deviation along dimension m
  asum      - the sum along dimension m
  ksdensity - the kernel density estimate

_Time series analysis

  bartlett  - M-point Bartlett window
  blackman  - M-point Blackman window
  cohere    - the coherence using average periodiogram
  csd       - the cross spectral density using average periodiogram
  fft       - the fast Fourier transform of vector x
  hamming   - M-point Hamming window
  hanning   - M-point Hanning window
  hist      - compute the histogram of x
  kaiser    - M length Kaiser window
  psd       - the power spectral density using average periodiogram
  sinc      - the sinc function of array x

_Dates

  date2num  - convert python datetimes to numeric representation
  drange    - create an array of numbers for date plots
  num2date  - convert numeric type (float days since 0001) to datetime

_Other

  angle     - the angle of a complex array
  griddata  - interpolate irregularly distributed data to a regular grid
  load      - Deprecated--please use loadtxt.
  loadtxt   - load ASCII data into array.
  polyfit   - fit x, y to an n-th order polynomial
  polyval   - evaluate an n-th order polynomial
  roots     - the roots of the polynomial coefficients in p
  save      - Deprecated--please use savetxt.
  savetxt   - save an array to an ASCII file.
  trapz     - trapezoidal integration

__end

.. [*] MATLAB is a registered trademark of The MathWorks, Inc.


"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import warnings

from matplotlib.cbook import (
    flatten, exception_to_str, silent_list, iterable, dedent)

import matplotlib as mpl

from matplotlib.dates import (
    date2num, num2date, datestr2num, strpdate2num, drange, epoch2num,
    num2epoch, mx2num, DateFormatter, IndexDateFormatter, DateLocator,
    RRuleLocator, YearLocator, MonthLocator, WeekdayLocator, DayLocator,
    HourLocator, MinuteLocator, SecondLocator, rrule, MO, TU, WE, TH, FR,
    SA, SU, YEARLY, MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY, SECONDLY,
    relativedelta)

# bring all the symbols in so folks can import them from
# pylab in one fell swoop

## We are still importing too many things from mlab; more cleanup is needed.

from matplotlib.mlab import (
    amap, base_repr, binary_repr, bivariate_normal, center_matrix, csv2rec,
    demean, detrend, detrend_linear, detrend_mean, detrend_none, dist,
    dist_point_to_segment, distances_along_curve, entropy, exp_safe,
    fftsurr, find, frange, get_sparse_matrix, get_xyz_where, griddata,
    identity, inside_poly, is_closed_polygon, ispower2, isvector, l1norm,
    l2norm, log2, longest_contiguous_ones, longest_ones, movavg, norm_flat,
    normpdf, path_length, poly_below, poly_between, prctile, prctile_rank,
    rec2csv, rec_append_fields, rec_drop_fields, rec_join, rk4, rms_flat,
    segments_intersect, slopes, stineman_interp, vector_lengths,
    window_hanning, window_none)

from matplotlib import cbook, mlab, pyplot as plt
from matplotlib.pyplot import *

from numpy import *
from numpy.fft import *
from numpy.random import *
from numpy.linalg import *

import numpy as np
import numpy.ma as ma

# don't let numpy's datetime hide stdlib
import datetime

# This is needed, or bytes will be numpy.random.bytes from
# "from numpy.random import *" above
bytes = six.moves.builtins.bytes
