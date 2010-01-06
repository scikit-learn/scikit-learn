#! /usr/bin/env python
# Last Change: Fri Jan 23 12:00 PM 2009 J
"""This is a small python package to estimate Gaussian Mixtures Models
from data, using Expectation Maximization.

Maximum likelihood EM for mixture of Gaussian is implemented, with BIC computation
for number of cluster assessment."""

from os.path import join

DISTNAME    = 'em2'
VERSION     = '0.1'
DESCRIPTION = 'A python module for Expectation Maximization learning of mixtures pdf',
AUTHOR      = 'David Cournapeau',
AUTHOR_EMAIL= 'david@ar.media.kyoto-u.ac.jp',
#URL         ='http://ar.media.kyoto-u.ac.jp/members/david/softwares/em',

def configuration(parent_package='',top_path=None, package_name=DISTNAME):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name,parent_package,top_path, version=VERSION)

    config.add_data_dir('tests')
    config.add_extension('_lk',
                         sources=['_lk.c'])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
