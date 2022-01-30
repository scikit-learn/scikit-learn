# -*- coding: utf-8 -*-
# $Id: ko.py 8541 2020-08-22 22:16:25Z milde $
# Author: Thomas SJ Kang <thomas.kangsj@ujuc.kr>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Korean-language mappings for language-dependent features of Docutils.
"""

__docformat__ = 'reStructuredText'

labels = {
      # fixed: language-dependent
      'author': u'저자',
      'authors': u'저자들',
      'organization': u'조직',
      'address': u'주소',
      'contact': u'연락처',
      'version': u'버전',
      'revision': u'리비전',
      'status': u'상태',
      'date': u'날짜',
      'copyright': u'저작권',
      'dedication': u'헌정',
      'abstract': u'요약',
      'attention': u'집중!',
      'caution': u'주의!',
      'danger': u'!위험!',
      'error': u'오류',
      'hint': u'실마리',
      'important': u'중요한',
      'note': u'비고',
      'tip': u'팁',
      'warning': u'경고',
      'contents': u'목차'}
"""Mapping of node class name to label text."""

bibliographic_fields = {
      # language-dependent: fixed
      u'저자': 'author',
      u'저자들': 'authors',
      u'조직': 'organization',
      u'주소': 'address',
      u'연락처': 'contact',
      u'버전': 'version',
      u'리비전': 'revision',
      u'상태': 'status',
      u'날짜': 'date',
      u'저작권': 'copyright',
      u'헌정': 'dedication',
      u'요약': 'abstract'}
"""Korean to canonical name mapping for bibliographic fields."""

author_separators = [';', ',']
"""List of separator strings for the 'Authors' bibliographic field. Tried in
order."""
