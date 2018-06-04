# -*- coding: utf-8 -*-
"""
Testing the Sphinx compatibility shims
"""
from __future__ import division, absolute_import, print_function

import sphinx.util.console

from sphinx_gallery import sphinx_compatibility


def test_status_iterator(fakesphinxapp):
    for _ in sphinx_compatibility._app_status_iterator([1, 2, 3],
                                                       'summary',
                                                       length=3):
        pass

    assert len(fakesphinxapp.calls['status_iterator']) == 1
    call = fakesphinxapp.calls['status_iterator'][0]
    assert call.args == ([1, 2, 3], 'summary')
    assert 'color' not in call.kwargs
    assert 'colorfunc' not in call.kwargs
    assert call.kwargs['length'] == 3


def test_status_iterator_color(fakesphinxapp):
    for _ in sphinx_compatibility._app_status_iterator([1, 2, 3],
                                                       'summary',
                                                       color='green',
                                                       length=3):
        pass

    assert len(fakesphinxapp.calls['status_iterator']) == 1
    call = fakesphinxapp.calls['status_iterator'][0]
    assert call.args == ([1, 2, 3], 'summary')
    assert 'color' not in call.kwargs
    assert call.kwargs['colorfunc'] == sphinx.util.console.green
    assert call.kwargs['length'] == 3


def test_get_logger(fakesphinxapp):
    logger = sphinx_compatibility._app_get_logger('sphinx-gallery-tests')
    logger.error('error')
    logger.critical('critical')
    logger.warning('warning 1')
    logger.warning('warning 2', color='green')
    logger.info('info 1')
    logger.info('info 2', color='green')
    logger.verbose('verbose')
    logger.debug('debug')

    # Error + critical both go through warning:
    assert len(fakesphinxapp.calls['warn']) == 4
    error, critical, warning1, warning2 = fakesphinxapp.calls['warn']
    assert error.args == (sphinx.util.console.red('error'), )
    assert error.kwargs == {}
    assert critical.args == (sphinx.util.console.red('critical'), )
    assert critical.kwargs == {}
    assert warning1.args == ('warning 1', )
    assert warning1.kwargs == {}
    assert warning2.args == (sphinx.util.console.green('warning 2'), )
    assert warning2.kwargs == {}

    assert len(fakesphinxapp.calls['info']) == 2
    info1, info2 = fakesphinxapp.calls['info']
    assert info1.args == ('info 1', )
    assert info1.kwargs == {}
    assert info2.args == (sphinx.util.console.green('info 2'), )
    assert info2.kwargs == {}

    assert len(fakesphinxapp.calls['verbose']) == 1
    assert fakesphinxapp.calls['verbose'][0].args == ('verbose', )
    assert fakesphinxapp.calls['verbose'][0].kwargs == {}

    assert len(fakesphinxapp.calls['debug']) == 1
    assert fakesphinxapp.calls['debug'][0].args == ('debug', )
    assert fakesphinxapp.calls['debug'][0].kwargs == {}
