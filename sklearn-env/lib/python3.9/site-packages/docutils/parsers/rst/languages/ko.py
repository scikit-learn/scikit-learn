# -*- coding: utf-8 -*-
# $Id: ko.py 8541 2020-08-22 22:16:25Z milde $
# Author: Thomas SJ Kang <thomas.kangsj@ujuc.kr>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Korean-language mappings for language-dependent features of
reStructuredText.
"""

__docformat__ = 'reStructuredText'


directives = {
      # language-dependent: fixed
      u'집중': 'attention',
      u'주의': 'caution',
      u'코드': 'code',
      u'코드-블록': 'code',
      u'소스코드': 'code',
      u'위험': 'danger',
      u'오류': 'error',
      u'실마리': 'hint',
      u'중요한': 'important',
      u'비고': 'note',
      u'팁': 'tip',
      u'경고': 'warning',
      u'권고': 'admonition',
      u'사이드바': 'sidebar',
      u'주제': 'topic',
      u'라인-블록': 'line-block',
      u'파싱된-리터럴': 'parsed-literal',
      u'지시문': 'rubric',
      u'제명': 'epigraph',
      u'하이라이트': 'highlights',
      u'발췌문': 'pull-quote',
      u'합성어': 'compound',
      u'컨테이너': 'container',
      #u'질문': 'questions',
      u'표': 'table',
      u'csv-표': 'csv-table',
      u'list-표': 'list-table',
      #u'qa': 'questions',
      #u'faq': 'questions',
      u'메타': 'meta',
      u'수학': 'math',
      #u'이미지맵': 'imagemap',
      u'이미지': 'image',
      u'도표': 'figure',
      u'포함': 'include',
      'raw': 'raw',
      u'대신하다': 'replace',
      u'유니코드': 'unicode',
      u'날짜': 'date',
      u'클래스': 'class',
      u'역할': 'role',
      u'기본-역할': 'default-role',
      u'제목': 'title',
      u'내용': 'contents',
      'sectnum': 'sectnum',
      u'섹션-번호-매기기': 'sectnum',
      u'머리말': 'header',
      u'꼬리말': 'footer',
      #u'긱주': 'footnotes',
      #u'인용구': 'citations',
      u'목표-노트': 'target-notes',
      u'restructuredtext 테스트 지시어': 'restructuredtext-test-directive'}
"""Korean name to registered (in directives/__init__.py) directive name
mapping."""

roles = {
    # language-dependent: fixed
    u'약어': 'abbreviation',
    u'ab': 'abbreviation',
    u'두음문자': 'acronym',
    u'ac': 'acronym',
    u'코드': 'code',
    u'색인': 'index',
    u'i': 'index',
    u'다리-글자': 'subscript',
    u'sub': 'subscript',
    u'어깨-글자': 'superscript',
    u'sup': 'superscript',
    u'제목-참조': 'title-reference',
    u'제목': 'title-reference',
    u't': 'title-reference',
    u'pep-참조': 'pep-reference',
    u'pep': 'pep-reference',
    u'rfc-참조': 'rfc-reference',
    u'rfc': 'rfc-reference',
    u'강조': 'emphasis',
    u'굵게': 'strong',
    u'기울기': 'literal',
    u'수학': 'math',
    u'명명된-참조': 'named-reference',
    u'익명-참조': 'anonymous-reference',
    u'각주-참조': 'footnote-reference',
    u'인용-참조': 'citation-reference',
    u'대리-참조': 'substitution-reference',
    u'대상': 'target',
    u'uri-참조': 'uri-reference',
    u'uri': 'uri-reference',
    u'url': 'uri-reference',
    'raw': 'raw',}
"""Mapping of Korean role names to canonical role names for interpreted text.
"""
