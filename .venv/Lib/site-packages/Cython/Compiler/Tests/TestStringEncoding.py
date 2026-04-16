import unittest

import Cython.Compiler.StringEncoding as StringEncoding


class StringEncodingTest(unittest.TestCase):
    """
    Test the StringEncoding module.
    """
    def test_string_contains_lone_surrogates(self):
        self.assertFalse(StringEncoding.string_contains_lone_surrogates("abc"))
        self.assertFalse(StringEncoding.string_contains_lone_surrogates("\uABCD"))
        self.assertFalse(StringEncoding.string_contains_lone_surrogates("\N{SNOWMAN}"))

        self.assertTrue(StringEncoding.string_contains_lone_surrogates("\uD800\uDFFF"))
        obfuscated_surrogate_pair = ("\uDFFF" + "\uD800")[::-1]
        self.assertTrue(StringEncoding.string_contains_lone_surrogates(obfuscated_surrogate_pair))
        self.assertTrue(StringEncoding.string_contains_lone_surrogates("\uD800"))
        self.assertTrue(StringEncoding.string_contains_lone_surrogates("\uDFFF"))
        self.assertTrue(StringEncoding.string_contains_lone_surrogates("\uDFFF\uD800"))
        self.assertTrue(StringEncoding.string_contains_lone_surrogates("\uD800x\uDFFF"))
