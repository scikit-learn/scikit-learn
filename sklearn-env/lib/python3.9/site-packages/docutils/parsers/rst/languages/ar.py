# -*- coding: utf-8 -*-
# $Id: fa.py 4564 2016-08-10 11:48:42Z
# Author: Shahin <me@5hah.in>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Arabic-language mappings for language-dependent features of
reStructuredText.
"""

__docformat__ = 'reStructuredText'


directives = {
      # language-dependent: fixed
      u'تنبيه': u'attention',
      u'احتیاط': u'caution',
      u'كود': u'code',
      u'كود': u'code',
      u'كود': u'code',
      u'خطر': u'danger',
      u'خطأ': u'error',
      u'تلميح': u'hint',
      u'مهم': u'important',
      u'ملاحظة': u'note',
      u'نصيحة': u'tip',
      u'تحذير': u'warning',
      u'تذكير': u'admonition',
      u'شريط-جانبي': u'sidebar',
      u'موضوع': u'topic',
      u'قالب-سطري': u'line-block',
      u'لفظ-حرفي': u'parsed-literal',
      u'معيار': u'rubric',
      u'فكرة-الكتاب': u'epigraph',
      u'تمييز': u'highlights',
      u'نقل-قول': u'pull-quote',
      u'ترکیب': u'compound',
      u'وعاء': u'container',
      #'questions': u'questions',
      u'جدول': u'table',
      u'جدول-csv': u'csv-table',
      u'جدول-قوائم': u'list-table',
      #'qa': u'questions',
      #'faq': u'questions',
      u'ميتا': u'meta',
      u'رياضيات': u'math',
      #'imagemap': u'imagemap',
      u'صورة': u'image',
      u'رسم-توضيحي': u'figure',
      u'تضمين': u'include',
      u'خام': u'raw',
      u'تبديل': u'replace',
      u'یونیکد': u'unicode',
      u'تاریخ': u'date',
      u'كائن': u'class',
      u'قانون': u'role',
      u'قانون-افتراضي': u'default-role',
      u'عنوان': u'title',
      u'المحتوى': u'contents',
      u'رقم-الفصل': u'sectnum',
      u'رقم-القسم': u'sectnum',
      u'رأس-الصفحة': u'header',
      u'هامش': u'footer',
      #'footnotes': u'footnotes',
      #'citations': u'citations',
      u'': u'target-notes',
    }
"""Arabic name to registered (in directives/__init__.py) directive name
mapping."""

roles = {
    # language-dependent: fixed
    u'اختصار': u'abbreviation',
    u'اختزال': u'acronym',
    u'كود': u'code',
    u'فهرس': u'index',
    u'خفض': u'subscript',
    u'رفع': u'superscript',
    u'عنوان-مرجع': u'title-reference',
    u'مرجع-pep': u'pep-reference',
    u'rfc-مرجع': u'rfc-reference',
    u'تأكيد': u'emphasis',
    u'عريض': u'strong',
    u'لفظی': u'literal',
    u'رياضيات': u'math',
    u'مرجع-مسمى': u'named-reference',
    u'مرجع-مجهول': u'anonymous-reference',
    u'مرجع-هامشي': u'footnote-reference',
    u'مرجع-منقول': u'citation-reference',
    u'مرجع-معوض': u'substitution-reference',
    u'هدف': u'target',
    u'منبع-uri': u'uri-reference',
    u'uri': u'uri-reference',
    u'url': u'uri-reference',
    u'خام': u'raw',}
"""Mapping of Arabic role names to canonical role names for interpreted text.
"""
