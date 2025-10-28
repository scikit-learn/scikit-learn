import difflib
import glob
import gzip
import os
import sys
import tempfile
import unittest

import Cython.Build.Dependencies
import Cython.Compiler.Main
import Cython.Utils
from Cython.TestUtils import CythonTest


class TestCyCache(CythonTest):

    def setUp(self):
        CythonTest.setUp(self)
        self.temp_dir = tempfile.mkdtemp(
            prefix='cycache-test',
            dir='TEST_TMP' if os.path.isdir('TEST_TMP') else None)
        self.src_dir = tempfile.mkdtemp(prefix='src', dir=self.temp_dir)
        self.cache_dir = tempfile.mkdtemp(prefix='cache', dir=self.temp_dir)

    def cache_files(self, file_glob):
        return glob.glob(os.path.join(self.cache_dir, file_glob))

    def fresh_cythonize(self, *args, **kwargs):
        Cython.Utils.clear_function_caches()
        Cython.Build.Dependencies._dep_tree = None  # discard method caches
        Cython.Build.Dependencies.cythonize(*args, **kwargs)

    def fresh_compile(self, *args, **kwargs):
        Cython.Utils.clear_function_caches()
        Cython.Compiler.Main.compile(*args, **kwargs)


    def _test_cycache_switch(self, compilation_method):
        content1 = 'value = 1\n'
        content2 = 'value = 2\n'
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'

        with open(a_pyx, 'w') as f:
            f.write(content1)

        compilation_method(a_pyx, cache=self.cache_dir)
        compilation_method(a_pyx, cache=self.cache_dir)

        self.assertEqual(1, len(self.cache_files('a.c*')))
        with open(a_c) as f:
            a_contents1 = f.read()
        os.unlink(a_c)

        with open(a_pyx, 'w') as f:
            f.write(content2)

        compilation_method(a_pyx, cache=self.cache_dir)

        with open(a_c) as f:
            a_contents2 = f.read()
        os.unlink(a_c)

        self.assertNotEqual(a_contents1, a_contents2, 'C file not changed!')
        self.assertEqual(2, len(self.cache_files('a.c*')))

        with open(a_pyx, 'w') as f:
            f.write(content1)

        compilation_method(a_pyx, cache=self.cache_dir)

        self.assertEqual(2, len(self.cache_files('a.c*')))
        with open(a_c) as f:
            a_contents = f.read()
        self.assertEqual(
            a_contents, a_contents1,
            msg='\n'.join(list(difflib.unified_diff(
                a_contents.split('\n'), a_contents1.split('\n')))[:10]))

    def test_cycache_switch_cythonize(self):
        self._test_cycache_switch(self.fresh_cythonize)

    def test_cycache_switch_compile(self):
        self._test_cycache_switch(self.fresh_compile)

    def _test_cycache_uses_cache(self, compilation_method):
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'
        with open(a_pyx, 'w') as f:
            f.write('pass')

        compilation_method(a_pyx, cache=self.cache_dir)

        a_cache = os.path.join(self.cache_dir, os.listdir(self.cache_dir)[0])
        with gzip.GzipFile(a_cache, 'wb') as gzipfile:
            gzipfile.write(b'fake stuff')
        os.unlink(a_c)

        compilation_method(a_pyx, cache=self.cache_dir)

        with open(a_c) as f:
            a_contents = f.read()
        self.assertEqual(a_contents, 'fake stuff',
                         'Unexpected contents: %s...' % a_contents[:100])


    def test_cycache_uses_cache_cythonize(self):
        self._test_cycache_uses_cache(self.fresh_cythonize)

    def test_cycache_uses_cache_compile(self):
        self._test_cycache_uses_cache(self.fresh_compile)

    def _test_cycache_annotation(self, compilation_method):
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'
        a_html = a_pyx[:-4] + '.html'
        with open(a_pyx, 'w') as f:
            f.write('pass')

        compilation_method(a_pyx, cache=self.cache_dir, annotate='default')
        self.assertTrue(os.path.exists(a_html), a_html)
        os.unlink(a_html)
        os.unlink(a_c)
        compilation_method(a_pyx, cache=self.cache_dir, annotate='default')
        self.assertTrue(os.path.exists(a_html), a_html)

    def test_cycache_annotation_cythonize(self):
        self._test_cycache_annotation(self.fresh_cythonize)

    def test_cycache_annotation_compile(self):
        self._test_cycache_annotation(self.fresh_compile)

    def _test_multi_file_output(self, compilation_method):
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'
        a_h = a_pyx[:-4] + '.h'
        a_api_h = a_pyx[:-4] + '_api.h'
        with open(a_pyx, 'w') as f:
            f.write('cdef public api int foo(int x): return x\n')

        compilation_method(a_pyx, cache=self.cache_dir)

        expected = [a_c, a_h, a_api_h]
        for output in expected:
            self.assertTrue(os.path.exists(output), output)
            os.unlink(output)

        compilation_method(a_pyx, cache=self.cache_dir)

        for output in expected:
            self.assertTrue(os.path.exists(output), output)

    def test_multi_file_output_cythonize(self):
        self._test_multi_file_output(self.fresh_cythonize)

    def test_multi_file_output_compile(self):
        self._test_multi_file_output(self.fresh_compile)

    def _test_options_invalidation(self, compilation_method):
        hash_pyx = os.path.join(self.src_dir, 'options.pyx')
        hash_c = hash_pyx[:-len('.pyx')] + '.c'
        hash_cpp = hash_pyx[:-len('.pyx')] + '.cpp'

        with open(hash_pyx, 'w') as f:
            f.write('pass')
        compilation_method(hash_pyx, cache=self.cache_dir, cplus=False)
        self.assertEqual(1, len(self.cache_files('options.c*')))

        os.unlink(hash_c)

        compilation_method(hash_pyx, cache=self.cache_dir, cplus=True)

        self.assertEqual(2, len(self.cache_files('options.c*')))

        try:
            os.unlink(hash_c)
        except FileNotFoundError:
            # fresh_cythonize() produces .c file, fresh_compile produces .cpp file
            os.unlink(hash_cpp)

        compilation_method(hash_pyx, cache=self.cache_dir, cplus=False, show_version=False)

        self.assertEqual(2, len(self.cache_files('options.c*')))

        os.unlink(hash_c)

        compilation_method(hash_pyx, cache=self.cache_dir, cplus=False, show_version=True)

        self.assertEqual(2, len(self.cache_files('options.c*')))
    def test_options_invalidation_cythonize(self):
        self._test_options_invalidation(self.fresh_cythonize)

    def test_options_invalidation_compile(self):
        self._test_options_invalidation(self.fresh_compile)
