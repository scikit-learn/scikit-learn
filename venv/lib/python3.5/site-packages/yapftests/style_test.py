# -*- coding: utf-8 -*-
# Copyright 2015 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for yapf.style."""

import shutil
import tempfile
import textwrap
import unittest

from yapf.yapflib import style

from yapftests import utils


class UtilsTest(unittest.TestCase):

  def testContinuationAlignStyleStringConverter(self):
    self.assertEqual(style._ContinuationAlignStyleStringConverter(''), 'SPACE')
    self.assertEqual(
        style._ContinuationAlignStyleStringConverter('space'), 'SPACE')
    self.assertEqual(
        style._ContinuationAlignStyleStringConverter('fixed'), 'FIXED')
    self.assertEqual(
        style._ContinuationAlignStyleStringConverter('valign-right'),
        'VALIGN-RIGHT')
    with self.assertRaises(ValueError) as ctx:
      style._ContinuationAlignStyleStringConverter('blahblah')
    self.assertIn("unknown continuation align style: 'blahblah'",
                  str(ctx.exception))

  def testStringListConverter(self):
    self.assertEqual(style._StringListConverter('foo, bar'), ['foo', 'bar'])
    self.assertEqual(style._StringListConverter('foo,bar'), ['foo', 'bar'])
    self.assertEqual(style._StringListConverter('  foo'), ['foo'])
    self.assertEqual(
        style._StringListConverter('joe  ,foo,  bar'), ['joe', 'foo', 'bar'])

  def testBoolConverter(self):
    self.assertEqual(style._BoolConverter('true'), True)
    self.assertEqual(style._BoolConverter('1'), True)
    self.assertEqual(style._BoolConverter('false'), False)
    self.assertEqual(style._BoolConverter('0'), False)


def _LooksLikeChromiumStyle(cfg):
  return (cfg['INDENT_WIDTH'] == 2 and
          cfg['BLANK_LINE_BEFORE_NESTED_CLASS_OR_DEF'])


def _LooksLikeGoogleStyle(cfg):
  return (cfg['INDENT_WIDTH'] == 4 and
          cfg['BLANK_LINE_BEFORE_NESTED_CLASS_OR_DEF'])


def _LooksLikePEP8Style(cfg):
  return (cfg['INDENT_WIDTH'] == 4 and
          not cfg['BLANK_LINE_BEFORE_NESTED_CLASS_OR_DEF'])


def _LooksLikeFacebookStyle(cfg):
  return cfg['INDENT_WIDTH'] == 4 and cfg['DEDENT_CLOSING_BRACKETS']


class PredefinedStylesByNameTest(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    style.SetGlobalStyle(style.CreatePEP8Style())

  def testDefault(self):
    # default is PEP8
    cfg = style.CreateStyleFromConfig(None)
    self.assertTrue(_LooksLikePEP8Style(cfg))

  def testPEP8ByName(self):
    for pep8_name in ('PEP8', 'pep8', 'Pep8'):
      cfg = style.CreateStyleFromConfig(pep8_name)
      self.assertTrue(_LooksLikePEP8Style(cfg))

  def testGoogleByName(self):
    for google_name in ('google', 'Google', 'GOOGLE'):
      cfg = style.CreateStyleFromConfig(google_name)
      self.assertTrue(_LooksLikeGoogleStyle(cfg))

  def testChromiumByName(self):
    for chromium_name in ('chromium', 'Chromium', 'CHROMIUM'):
      cfg = style.CreateStyleFromConfig(chromium_name)
      self.assertTrue(_LooksLikeChromiumStyle(cfg))

  def testFacebookByName(self):
    for fb_name in ('facebook', 'FACEBOOK', 'Facebook'):
      cfg = style.CreateStyleFromConfig(fb_name)
      self.assertTrue(_LooksLikeFacebookStyle(cfg))


class StyleFromFileTest(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.test_tmpdir = tempfile.mkdtemp()
    style.SetGlobalStyle(style.CreatePEP8Style())

  @classmethod
  def tearDownClass(cls):
    shutil.rmtree(cls.test_tmpdir)

  def testDefaultBasedOnStyle(self):
    cfg = textwrap.dedent(u'''\
        [style]
        continuation_indent_width = 20
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikePEP8Style(cfg))
      self.assertEqual(cfg['CONTINUATION_INDENT_WIDTH'], 20)

  def testDefaultBasedOnPEP8Style(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = pep8
        continuation_indent_width = 40
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikePEP8Style(cfg))
      self.assertEqual(cfg['CONTINUATION_INDENT_WIDTH'], 40)

  def testDefaultBasedOnChromiumStyle(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = chromium
        continuation_indent_width = 30
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikeChromiumStyle(cfg))
      self.assertEqual(cfg['CONTINUATION_INDENT_WIDTH'], 30)

  def testDefaultBasedOnGoogleStyle(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = google
        continuation_indent_width = 20
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikeGoogleStyle(cfg))
      self.assertEqual(cfg['CONTINUATION_INDENT_WIDTH'], 20)

  def testDefaultBasedOnFacebookStyle(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = facebook
        continuation_indent_width = 20
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikeFacebookStyle(cfg))
      self.assertEqual(cfg['CONTINUATION_INDENT_WIDTH'], 20)

  def testBoolOptionValue(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = chromium
        SPLIT_BEFORE_NAMED_ASSIGNS=False
        split_before_logical_operator = true
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikeChromiumStyle(cfg))
      self.assertEqual(cfg['SPLIT_BEFORE_NAMED_ASSIGNS'], False)
      self.assertEqual(cfg['SPLIT_BEFORE_LOGICAL_OPERATOR'], True)

  def testStringListOptionValue(self):
    cfg = textwrap.dedent(u'''\
        [style]
        based_on_style = chromium
        I18N_FUNCTION_CALL = N_, V_, T_
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      cfg = style.CreateStyleFromConfig(filepath)
      self.assertTrue(_LooksLikeChromiumStyle(cfg))
      self.assertEqual(cfg['I18N_FUNCTION_CALL'], ['N_', 'V_', 'T_'])

  def testErrorNoStyleFile(self):
    with self.assertRaisesRegexp(style.StyleConfigError,
                                 'is not a valid style or file path'):
      style.CreateStyleFromConfig('/8822/xyznosuchfile')

  def testErrorNoStyleSection(self):
    cfg = textwrap.dedent(u'''\
        [s]
        indent_width=2
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      with self.assertRaisesRegexp(style.StyleConfigError,
                                   'Unable to find section'):
        style.CreateStyleFromConfig(filepath)

  def testErrorUnknownStyleOption(self):
    cfg = textwrap.dedent(u'''\
        [style]
        indent_width=2
        hummus=2
        ''')
    with utils.TempFileContents(self.test_tmpdir, cfg) as filepath:
      with self.assertRaisesRegexp(style.StyleConfigError,
                                   'Unknown style option'):
        style.CreateStyleFromConfig(filepath)


class StyleFromDict(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    style.SetGlobalStyle(style.CreatePEP8Style())

  def testDefaultBasedOnStyle(self):
    config_dict = {
        'based_on_style': 'pep8',
        'indent_width': 2,
        'blank_line_before_nested_class_or_def': True
    }
    cfg = style.CreateStyleFromConfig(config_dict)
    self.assertTrue(_LooksLikeChromiumStyle(cfg))
    self.assertEqual(cfg['INDENT_WIDTH'], 2)

  def testDefaultBasedOnStyleBadDict(self):
    self.assertRaisesRegexp(style.StyleConfigError, 'Unknown style option',
                            style.CreateStyleFromConfig,
                            {'based_on_styl': 'pep8'})
    self.assertRaisesRegexp(style.StyleConfigError, 'not a valid',
                            style.CreateStyleFromConfig,
                            {'INDENT_WIDTH': 'FOUR'})


class StyleFromCommandLine(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    style.SetGlobalStyle(style.CreatePEP8Style())

  def testDefaultBasedOnStyle(self):
    cfg = style.CreateStyleFromConfig(
        '{based_on_style: pep8,'
        ' indent_width: 2,'
        ' blank_line_before_nested_class_or_def: True}')
    self.assertTrue(_LooksLikeChromiumStyle(cfg))
    self.assertEqual(cfg['INDENT_WIDTH'], 2)

  def testDefaultBasedOnStyleNotStrict(self):
    cfg = style.CreateStyleFromConfig(
        '{based_on_style : pep8'
        ' ,indent_width=2'
        ' blank_line_before_nested_class_or_def:True}')
    self.assertTrue(_LooksLikeChromiumStyle(cfg))
    self.assertEqual(cfg['INDENT_WIDTH'], 2)

  def testDefaultBasedOnExplicitlyUnicodeTypeString(self):
    cfg = style.CreateStyleFromConfig(u'{}')
    self.assertIsInstance(cfg, dict)

  def testDefaultBasedOnDetaultTypeString(self):
    cfg = style.CreateStyleFromConfig('{}')
    self.assertIsInstance(cfg, dict)

  def testDefaultBasedOnStyleBadString(self):
    self.assertRaisesRegexp(style.StyleConfigError, 'Unknown style option',
                            style.CreateStyleFromConfig,
                            '{based_on_styl: pep8}')
    self.assertRaisesRegexp(style.StyleConfigError, 'not a valid',
                            style.CreateStyleFromConfig, '{INDENT_WIDTH: FOUR}')
    self.assertRaisesRegexp(style.StyleConfigError, 'Invalid style dict',
                            style.CreateStyleFromConfig,
                            '{based_on_style: pep8')


class StyleHelp(unittest.TestCase):

  def testHelpKeys(self):
    settings = sorted(style.Help())
    expected = sorted(style._style)
    self.assertListEqual(settings, expected)


if __name__ == '__main__':
  unittest.main()
