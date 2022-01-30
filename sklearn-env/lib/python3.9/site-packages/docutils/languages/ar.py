# -*- coding: utf-8 -*-
# $Id: fa.py 4564 2016-08-10 11:48:42Z
# Author: Shahin <me@5hah.in>
# Copyright: This module has been placed in the public domain.

# New language mappings are welcome.  Before doing a new translation, please
# read <http://docutils.sf.net/docs/howto/i18n.html>.  Two files must be
# translated for each language: one in docutils/languages, the other in
# docutils/parsers/rst/languages.

"""
Arabic-language mappings for language-dependent features of Docutils.
"""

__docformat__ = 'reStructuredText'

labels = {
      # fixed: language-dependent
      u'author': u'المؤلف',
      u'authors': u'المؤلفون',
      u'organization': u'التنظيم',
      u'address': u'العنوان',
      u'contact': u'اتصل',
      u'version': u'نسخة',
      u'revision': u'مراجعة',
      u'status': u'الحالة',
      u'date': u'تاریخ',
      u'copyright': u'الحقوق',
      u'dedication': u'إهداء',
      u'abstract': u'ملخص',
      u'attention': u'تنبيه',
      u'caution': u'احتیاط',
      u'danger': u'خطر',
      u'error': u'خطأ',
      u'hint': u'تلميح',
      u'important': u'مهم',
      u'note': u'ملاحظة',
      u'tip': u'نصيحة',
      u'warning': u'تحذير',
      u'contents': u'المحتوى'}
"""Mapping of node class name to label text."""

bibliographic_fields = {
      # language-dependent: fixed
      u'مؤلف': u'author',
      u'مؤلفون': u'authors',
      u'التنظيم': u'organization',
      u'العنوان': u'address',
      u'اتصل': u'contact',
      u'نسخة': u'version',
      u'مراجعة': u'revision',
      u'الحالة': u'status',
      u'تاریخ': u'date',
      u'الحقوق': u'copyright',
      u'إهداء': u'dedication',
      u'ملخص': u'abstract'}
"""Arabic (lowcased) to canonical name mapping for bibliographic fields."""

author_separators = [u'؛', u'،']
"""List of separator strings for the 'Authors' bibliographic field. Tried in
order."""
