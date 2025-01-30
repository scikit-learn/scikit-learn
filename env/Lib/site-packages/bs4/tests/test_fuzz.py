"""This file contains test cases reported by third parties using
fuzzing tools, primarily from Google's oss-fuzz project. Some of these
represent real problems with Beautiful Soup, but many are problems in
libraries that Beautiful Soup depends on, and many of the test cases
represent different ways of triggering the same problem.

Grouping these test cases together makes it easy to see which test
cases represent the same problem, and puts the test cases in close
proximity to code that can trigger the problems.
"""
import os
import pytest
from bs4 import (
    BeautifulSoup,
    ParserRejectedMarkup,
)
try:
    from soupsieve.util import SelectorSyntaxError
    import lxml
    import html5lib
    fully_fuzzable = True
except ImportError:
    fully_fuzzable = False
    

@pytest.mark.skipif(not fully_fuzzable, reason="Prerequisites for fuzz tests are not installed.")
class TestFuzz(object):

    # Test case markup files from fuzzers are given this extension so
    # they can be included in builds.
    TESTCASE_SUFFIX = ".testcase"

    # Copied 20230512 from
    # https://github.com/google/oss-fuzz/blob/4ac6a645a197a695fe76532251feb5067076b3f3/projects/bs4/bs4_fuzzer.py
    #
    # Copying the code lets us precisely duplicate the behavior of
    # oss-fuzz.  The downside is that this code changes over time, so
    # multiple copies of the code must be kept around to run against
    # older tests. I'm not sure what to do about this, but I may
    # retire old tests after a time.
    def fuzz_test_with_css(self, filename):
        data = self.__markup(filename)
        parsers = ['lxml-xml', 'html5lib', 'html.parser', 'lxml']
        try:
            idx = int(data[0]) % len(parsers)
        except ValueError:
            return

        css_selector, data = data[1:10], data[10:]

        try:
            soup = BeautifulSoup(data[1:], features=parsers[idx])
        except ParserRejectedMarkup:
            return
        except ValueError:
            return

        list(soup.find_all(True))
        try:
            soup.css.select(css_selector.decode('utf-8', 'replace'))
        except SelectorSyntaxError:
            return
        soup.prettify()
    
    # This class of error has been fixed by catching a less helpful
    # exception from html.parser and raising ParserRejectedMarkup
    # instead.
    @pytest.mark.parametrize(
        "filename", [
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5703933063462912",
            "crash-ffbdfa8a2b26f13537b68d3794b0478a4090ee4a",
        ]
    )
    def test_rejected_markup(self, filename):
        markup = self.__markup(filename)
        with pytest.raises(ParserRejectedMarkup):
            BeautifulSoup(markup, 'html.parser')
            
    # This class of error has to do with very deeply nested documents
    # which overflow the Python call stack when the tree is converted
    # to a string. This is an issue with Beautiful Soup which was fixed
    # as part of [bug=1471755].
    #
    # These test cases are in the older format that doesn't specify
    # which parser to use or give a CSS selector.
    @pytest.mark.parametrize(
        "filename", [
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5984173902397440",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5167584867909632",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-6124268085182464",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-6450958476902400",
        ]
    )
    def test_deeply_nested_document_without_css(self, filename):
        # Parsing the document and encoding it back to a string is
        # sufficient to demonstrate that the overflow problem has
        # been fixed.
        markup = self.__markup(filename)
        BeautifulSoup(markup, 'html.parser').encode()

    # This class of error has to do with very deeply nested documents
    # which overflow the Python call stack when the tree is converted
    # to a string. This is an issue with Beautiful Soup which was fixed
    # as part of [bug=1471755].
    @pytest.mark.parametrize(
        "filename", [
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5000587759190016",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5375146639360000",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5492400320282624",
        ]
    )
    def test_deeply_nested_document(self, filename): 
       self.fuzz_test_with_css(filename)
        
    @pytest.mark.parametrize(
        "filename", [
            "clusterfuzz-testcase-minimized-bs4_fuzzer-4670634698080256",
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5270998950477824",
        ]
    )
    def test_soupsieve_errors(self, filename):
        self.fuzz_test_with_css(filename)
        
    # This class of error represents problems with html5lib's parser,
    # not Beautiful Soup. I use
    # https://github.com/html5lib/html5lib-python/issues/568 to notify
    # the html5lib developers of these issues.
    #
    # These test cases are in the older format that doesn't specify
    # which parser to use or give a CSS selector.
    @pytest.mark.skip(reason="html5lib-specific problems")
    @pytest.mark.parametrize(
        "filename", [
            # b"""ÿ<!DOCTyPEV PUBLIC'''Ð'"""
            "clusterfuzz-testcase-minimized-bs4_fuzzer-4818336571064320",

            # b')<a><math><TR><a><mI><a><p><a>'
            "clusterfuzz-testcase-minimized-bs4_fuzzer-4999465949331456",

            # b'-<math><sElect><mi><sElect><sElect>'
            "clusterfuzz-testcase-minimized-bs4_fuzzer-5843991618256896",
           
            # b'ñ<table><svg><html>'
            "clusterfuzz-testcase-minimized-bs4_fuzzer-6241471367348224",

            # <TABLE>, some ^@ characters, some <math> tags.
            "clusterfuzz-testcase-minimized-bs4_fuzzer-6600557255327744",

            # Nested table
            "crash-0d306a50c8ed8bcd0785b67000fcd5dea1d33f08"
        ]
    )
    def test_html5lib_parse_errors_without_css(self, filename):
        markup = self.__markup(filename)
        print(BeautifulSoup(markup, 'html5lib').encode())

    # This class of error represents problems with html5lib's parser,
    # not Beautiful Soup. I use
    # https://github.com/html5lib/html5lib-python/issues/568 to notify
    # the html5lib developers of these issues.
    @pytest.mark.skip(reason="html5lib-specific problems")
    @pytest.mark.parametrize(
        "filename", [
            # b'-      \xff\xff  <math>\x10<select><mi><select><select>t'
            "clusterfuzz-testcase-minimized-bs4_fuzzer-6306874195312640",
        ]
    )
    def test_html5lib_parse_errors(self, filename):
        self.fuzz_test_with_css(filename)
        
    def __markup(self, filename):
        if not filename.endswith(self.TESTCASE_SUFFIX):
            filename += self.TESTCASE_SUFFIX
        this_dir = os.path.split(__file__)[0]
        path = os.path.join(this_dir, 'fuzz', filename)
        return open(path, 'rb').read()
