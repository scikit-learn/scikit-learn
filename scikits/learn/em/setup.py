#! /usr/bin/env python
# Last Change: Mon Jul 02 02:00 PM 2007 J
# TODO:
#   - check how to handle cmd line build options with distutils and use
#   it in the building process

"""This is a small python package to estimate Gaussian Mixtures Models
from data, using Expectation Maximization.

Maximum likelihood EM for mixture of Gaussian is implemented, with BIC computation
for number of cluster assessment.

There is also an experimental online EM version (the EM is updated for each new
sample), and I plan to add Variational Bayes and/or MCMC support for Bayesian approach
for estimating meta parameters of mixtures. """
from os.path import join

def configuration(parent_package='',top_path=None, package_name='em'):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('em',parent_package, top_path)
    config.add_data_dir('examples')
    config.add_data_dir('tests')
    config.add_data_dir('profile_data')
    config.add_extension('c_gden',
                         sources=[join('src', 'c_gden.c')])
    config.add_extension('_rawden',
                         sources=[join('src', 'pure_den.c')])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
