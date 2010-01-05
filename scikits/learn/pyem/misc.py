# Last Change: Mon Jul 02 06:00 PM 2007 J

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
MAX_DBL_DEV    = 1e-10

## # max conditional number allowed
## _MAX_COND       = 1e8
## _MIN_INV_COND   = 1/_MAX_COND
## 
## # Default alpha for regularization
## _DEF_ALPHA  = 1e-1
## 
## # Default min delta for regularization
## _MIN_DBL_DELTA  = 1e-5
## 

class curry:
    def __init__(self, fun, *args, **kwargs):
        self.fun = fun
        self.pending = args[:]
        self.kwargs = kwargs.copy()

    def __call__(self, *args, **kwargs):
        if kwargs and self.kwargs:
            kw = self.kwargs.copy()
            kw.update(kwargs)
        else:
            kw = kwargs or self.kwargs

        return self.fun(*(self.pending + args), **kw)
