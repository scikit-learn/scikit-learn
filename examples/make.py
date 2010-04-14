#!/usr/bin/env python
import sys, os, glob
import matplotlib
import IPython.Shell
#matplotlib.rcdefaults()
matplotlib.use('Agg')

mplshell = IPython.Shell.MatplotlibShell('mpl')

formats = [('png', 100),
           ('hires.png', 200),
           ('pdf', 72)]

def figs():
    print 'making figs'
    import matplotlib.pyplot as plt
    for fname in glob.glob('*.py'):
        if fname.split('/')[-1] == __file__.split('/')[-1]: continue
        basename, ext = os.path.splitext(fname)
        imagefiles = dict([('%s.%s'%(basename, format), dpi)
                           for format, dpi in formats])
        all_exists = True
        for imagefile in imagefiles:
            if not os.path.exists(imagefile):
                all_exists = False
                break

        if all_exists:
            print '    already have %s'%fname
        else:
            print '    building %s'%fname
            plt.close('all')    # we need to clear between runs
            mplshell.magic_run(basename)
            for imagefile, dpi in imagefiles.iteritems():
                # todo: this will get called even if the run script
                # fails and exits, thus creating a stub pdf and png
                # iles preventing them from getting built successfully
                # later
                plt.savefig(imagefile, dpi=dpi)
    print 'all figures made'


def clean():
    patterns = (['#*', '*~', '*pyc'] +
                ['*.%s' % format for format, dpi in formats])
    for pattern in patterns:
        for fname in glob.glob(pattern):
            os.remove(fname)
    print 'all clean'



def all():
    figs()

funcd = {'figs':figs,
         'clean':clean,
         'all':all,
         }

if len(sys.argv)>1:
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are'%(
                    arg, funcd.keys()))
        func()
else:
    all()




