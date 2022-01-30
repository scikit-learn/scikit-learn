# -*- coding: utf-8 -*-
# $Id: sv.py 8012 2017-01-03 23:08:19Z milde $
# Author: Adam Chodorowski <chodorowski@users.sourceforge.net>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Swedish language mappings for language-dependent features of reStructuredText.
"""

__docformat__ = 'reStructuredText'

directives = {
      u'observera': 'attention',
      u'akta': 'caution', # also 'försiktigt'
      u'kod': 'code',
      u'fara': 'danger',
      u'fel': 'error',
      u'vink': 'hint', # also 'hint'
      u'viktigt': 'important',
      u'notera': 'note',
      u'tips': 'tip',
      u'varning': 'warning',
      u'anmärkning': 'admonition', # literal 'tillrättavisning', 'förmaning'
      u'sidorad': 'sidebar',
      u'ämne': 'topic',
      u'tema': 'topic',
      u'rad-block': 'line-block',
      u'parsed-literal (translation required)': 'parsed-literal', # 'tolkad-bokstavlig'?
      u'rubrik': 'rubric',
      u'epigraf': 'epigraph',
      u'höjdpunkter': 'highlights',
      u'pull-quote (translation required)': 'pull-quote',
      u'sammansatt': 'compound',
      u'container': 'container',
      # u'frågor': 'questions',
      # NOTE: A bit long, but recommended by http://www.nada.kth.se/dataterm/:
      # u'frågor-och-svar': 'questions',
      # u'vanliga-frågor': 'questions',
      u'tabell': 'table',
      u'csv-tabell': 'csv-table',
      u'list-tabell': 'list-table',
      u'meta': 'meta',
      u'matematik': 'math',
      # u'bildkarta': 'imagemap',   # FIXME: Translation might be too literal.
      u'bild': 'image',
      u'figur': 'figure',
      u'inkludera': 'include',
      u'rå': 'raw',
      u'ersätta': 'replace',
      u'unicode': 'unicode',
      u'datum': 'date',
      u'klass': 'class',
      u'roll': 'role',
      u'standardroll': 'default-role',
      u'titel': 'title',
      u'innehåll': 'contents',
      u'sektionsnumrering': 'sectnum',
      u'target-notes (translation required)': 'target-notes',
      u'sidhuvud': 'header',
      u'sidfot': 'footer',
      # u'fotnoter': 'footnotes',
      # u'citeringar': 'citations',
      }
"""Swedish name to registered (in directives/__init__.py) directive name
mapping."""

roles = {
      u'förkortning': 'abbreviation',
      u'akronym': 'acronym',
      u'kod': 'code',
      u'index': 'index',
      u'nedsänkt': 'subscript',
      u'upphöjd': 'superscript',
      u'titel-referens': 'title-reference',
      u'pep-referens': 'pep-reference',
      u'rfc-referens': 'rfc-reference',
      u'betoning': 'emphasis',
      u'stark': 'strong',
      u'bokstavlig': 'literal', # also 'ordagranna'
      u'matematik': 'math',
      u'namngiven-referens': 'named-reference',
      u'anonym-referens': 'anonymous-reference',
      u'fotnot-referens': 'footnote-reference',
      u'citat-referens': 'citation-reference',
      u'ersättnings-referens': 'substitution-reference',
      u'mål': 'target',
      u'uri-referens': 'uri-reference',
      u'rå': 'raw',}
"""Mapping of Swedish role names to canonical role names for interpreted text.
"""
