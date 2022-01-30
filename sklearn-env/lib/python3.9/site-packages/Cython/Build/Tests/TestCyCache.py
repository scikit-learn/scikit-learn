import difflib
import glob
import gzip
import os
import tempfile

import Cython.Build.Dependencies
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

    def test_cycache_switch(self):
        content1 = 'value = 1\n'
        content2 = 'value = 2\n'
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'

        open(a_pyx, 'w').write(content1)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        self.assertEqual(1, len(self.cache_files('a.c*')))
        a_contents1 = open(a_c).read()
        os.unlink(a_c)

        open(a_pyx, 'w').write(content2)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        a_contents2 = open(a_c).read()
        os.unlink(a_c)

        self.assertNotEqual(a_contents1, a_contents2, 'C file not changed!')
        self.assertEqual(2, len(self.cache_files('a.c*')))

        open(a_pyx, 'w').write(content1)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        self.assertEqual(2, len(self.cache_files('a.c*')))
        a_contents = open(a_c).read()
        self.assertEqual(
            a_contents, a_contents1,
            msg='\n'.join(list(difflib.unified_diff(
                a_contents.split('\n'), a_contents1.split('\n')))[:10]))

    def test_cycache_uses_cache(self):
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'
        open(a_pyx, 'w').write('pass')
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        a_cache = os.path.join(self.cache_dir, os.listdir(self.cache_dir)[0])
        gzip.GzipFile(a_cache, 'wb').write('fake stuff'.encode('ascii'))
        os.unlink(a_c)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        a_contents = open(a_c).read()
        self.assertEqual(a_contents, 'fake stuff',
                         'Unexpected contents: %s...' % a_contents[:100])

    def test_multi_file_output(self):
        a_pyx = os.path.join(self.src_dir, 'a.pyx')
        a_c = a_pyx[:-4] + '.c'
        a_h = a_pyx[:-4] + '.h'
        a_api_h = a_pyx[:-4] + '_api.h'
        open(a_pyx, 'w').write('cdef public api int foo(int x): return x\n')
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        expected = [a_c, a_h, a_api_h]
        for output in expected:
            self.assertTrue(os.path.exists(output), output)
            os.unlink(output)
        self.fresh_cythonize(a_pyx, cache=self.cache_dir)
        for output in expected:
            self.assertTrue(os.path.exists(output), output)

    def test_options_invalidation(self):
        hash_pyx = os.path.join(self.src_dir, 'options.pyx')
        hash_c = hash_pyx[:-len('.pyx')] + '.c'

        open(hash_pyx, 'w').write('pass')
        self.fresh_cythonize(hash_pyx, cache=self.cache_dir, cplus=False)
        self.assertEqual(1, len(self.cache_files('options.c*')))

        os.unlink(hash_c)
        self.fresh_cythonize(hash_pyx, cache=self.cache_dir, cplus=True)
        self.assertEqual(2, len(self.cache_files('options.c*')))

        os.unlink(hash_c)
        self.fresh_cythonize(hash_pyx, cache=self.cache_dir, cplus=False, show_version=False)
        self.assertEqual(2, len(self.cache_files('options.c*')))

        os.unlink(hash_c)
        self.fresh_cythonize(hash_pyx, cache=self.cache_dir, cplus=False, show_version=True)
        self.assertEqual(2, len(self.cache_files('options.c*')))
