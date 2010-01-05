# Last Change: Sat Jun 09 08:00 PM 2007 J

#========================================================
# Constants used throughout the module (def args, etc...)
#========================================================
# This is the default dimension for representing confidence ellipses
DEF_VIS_DIM = (0, 1)
DEF_ELL_NP = 100
DEF_LEVEL = 0.39
#=====================================================================
# "magic number", that is number used to control regularization and co
# Change them at your risk !
#=====================================================================

# max deviation allowed when comparing double (this is actually stupid,
# I should actually use a number of decimals)
_MAX_DBL_DEV    = 1e-10

# max conditional number allowed
_MAX_COND       = 1e8
_MIN_INV_COND   = 1/_MAX_COND

# Default alpha for regularization
_DEF_ALPHA  = 1e-1

# Default min delta for regularization
_MIN_DBL_DELTA  = 1e-5

