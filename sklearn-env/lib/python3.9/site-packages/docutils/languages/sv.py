# -*- coding: utf-8 -*-
# $Id: sv.py 8006 2016-12-22 23:02:44Z milde $
# Author: Adam Chodorowski <chodorowski@users.sourceforge.net>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Swedish language mappings for language-dependent features of Docutils.
"""

__docformat__ = 'reStructuredText'

labels = {
    'author':       u'Författare',
    'authors':      u'Författare',
    'organization': u'Organisation',
    'address':      u'Adress',
    'contact':      u'Kontakt',
    'version':      u'Version',
    'revision':     u'Revision',
    'status':       u'Status',
    'date':         u'Datum',
    'copyright':    u'Copyright',
    'dedication':   u'Dedikation',
    'abstract':     u'Sammanfattning',
    'attention':    u'Observera!',
    'caution':      u'Akta!', # 'Varning' already used for 'warning'
    'danger':       u'FARA!',
    'error':        u'Fel',
    'hint':         u'Vink',
    'important':    u'Viktigt',
    'note':         u'Notera',
    'tip':          u'Tips',
    'warning':      u'Varning',
    'contents':     u'Innehåll' }
"""Mapping of node class name to label text."""

bibliographic_fields = {
    # 'Author' and 'Authors' identical in Swedish; assume the plural:
    u'författare': 'authors',
    u' n/a':            'author',
    u'organisation':    'organization',
    u'adress':          'address',
    u'kontakt':         'contact',
    u'version':         'version',
    u'revision':        'revision',
    u'status':          'status',
    u'datum':           'date',
    u'copyright':       'copyright',
    u'dedikation':      'dedication', 
    u'sammanfattning':  'abstract' }
"""Swedish (lowcased) to canonical name mapping for bibliographic fields."""

author_separators = [';', ',']
"""List of separator strings for the 'Authors' bibliographic field. Tried in
order."""
