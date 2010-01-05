# Last Change: Fri Nov 10 10:00 AM 2006 J

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

