"""
Provides a collection of utilities for comparing (image) results.

"""
from __future__ import absolute_import, division, print_function

import six

import atexit
import functools
import hashlib
import itertools
import os
import re
import shutil
import sys
from tempfile import TemporaryFile

import numpy as np

import matplotlib
from matplotlib.compat import subprocess
from matplotlib.testing.exceptions import ImageComparisonFailure
from matplotlib import _png
from matplotlib import _get_cachedir
from matplotlib import cbook

__all__ = ['compare_float', 'compare_images', 'comparable_formats']


def make_test_filename(fname, purpose):
    """
    Make a new filename by inserting `purpose` before the file's
    extension.
    """
    base, ext = os.path.splitext(fname)
    return '%s-%s%s' % (base, purpose, ext)


def compare_float(expected, actual, relTol=None, absTol=None):
    """
    Fail if the floating point values are not close enough, with
    the given message.

    You can specify a relative tolerance, absolute tolerance, or both.

    """
    if relTol is None and absTol is None:
        raise ValueError("You haven't specified a 'relTol' relative "
                         "tolerance or a 'absTol' absolute tolerance "
                         "function argument. You must specify one.")
    msg = ""

    if absTol is not None:
        absDiff = abs(expected - actual)
        if absTol < absDiff:
            template = ['',
                        'Expected: {expected}',
                        'Actual:   {actual}',
                        'Abs diff: {absDiff}',
                        'Abs tol:  {absTol}']
            msg += '\n  '.join([line.format(**locals()) for line in template])

    if relTol is not None:
        # The relative difference of the two values.  If the expected value is
        # zero, then return the absolute value of the difference.
        relDiff = abs(expected - actual)
        if expected:
            relDiff = relDiff / abs(expected)

        if relTol < relDiff:
            # The relative difference is a ratio, so it's always unit-less.
            template = ['',
                        'Expected: {expected}',
                        'Actual:   {actual}',
                        'Rel diff: {relDiff}',
                        'Rel tol:  {relTol}']
            msg += '\n  '.join([line.format(**locals()) for line in template])

    return msg or None


def get_cache_dir():
    cachedir = _get_cachedir()
    if cachedir is None:
        raise RuntimeError('Could not find a suitable configuration directory')
    cache_dir = os.path.join(cachedir, 'test_cache')
    if not os.path.exists(cache_dir):
        try:
            cbook.mkdirs(cache_dir)
        except IOError:
            return None
    if not os.access(cache_dir, os.W_OK):
        return None
    return cache_dir


def get_file_hash(path, block_size=2 ** 20):
    md5 = hashlib.md5()
    with open(path, 'rb') as fd:
        while True:
            data = fd.read(block_size)
            if not data:
                break
            md5.update(data)

    if path.endswith('.pdf'):
        from matplotlib import checkdep_ghostscript
        md5.update(checkdep_ghostscript()[1].encode('utf-8'))
    elif path.endswith('.svg'):
        from matplotlib import checkdep_inkscape
        md5.update(checkdep_inkscape().encode('utf-8'))

    return md5.hexdigest()


def make_external_conversion_command(cmd):
    def convert(old, new):
        cmdline = cmd(old, new)
        pipe = subprocess.Popen(cmdline, universal_newlines=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate()
        errcode = pipe.wait()
        if not os.path.exists(new) or errcode:
            msg = "Conversion command failed:\n%s\n" % ' '.join(cmdline)
            if stdout:
                msg += "Standard output:\n%s\n" % stdout
            if stderr:
                msg += "Standard error:\n%s\n" % stderr
            raise IOError(msg)

    return convert


# Modified from https://bugs.python.org/issue25567.
_find_unsafe_bytes = re.compile(br'[^a-zA-Z0-9_@%+=:,./-]').search


def _shlex_quote_bytes(b):
    return (b if _find_unsafe_bytes(b) is None
            else b"'" + b.replace(b"'", b"'\"'\"'") + b"'")


class _SVGConverter(object):
    def __init__(self):
        self._proc = None
        # We cannot rely on the GC to trigger `__del__` at exit because
        # other modules (e.g. `subprocess`) may already have their globals
        # set to `None`, which make `proc.communicate` or `proc.terminate`
        # fail.  By relying on `atexit` we ensure the destructor runs before
        # `None`-setting occurs.
        atexit.register(self.__del__)

    def _read_to_prompt(self):
        """Did Inkscape reach the prompt without crashing?
        """
        stream = iter(functools.partial(self._proc.stdout.read, 1), b"")
        prompt = (b"\n", b">")
        n = len(prompt)
        its = itertools.tee(stream, n)
        for i, it in enumerate(its):
            next(itertools.islice(it, i, i), None)  # Advance `it` by `i`.
        while True:
            window = tuple(map(next, its))
            if len(window) != n:
                # Ran out of data -- one of the `next(it)` raised
                # StopIteration, so the tuple is shorter.
                return False
            if self._proc.poll() is not None:
                # Inkscape exited.
                return False
            if window == prompt:
                # Successfully read until prompt.
                return True

    def __call__(self, orig, dest):
        if (not self._proc  # First run.
                or self._proc.poll() is not None):  # Inkscape terminated.
            env = os.environ.copy()
            # If one passes e.g. a png file to Inkscape, it will try to
            # query the user for conversion options via a GUI (even with
            # `--without-gui`).  Unsetting `DISPLAY` prevents this (and causes
            # GTK to crash and Inkscape to terminate, but that'll just be
            # reported as a regular exception below).
            env.pop("DISPLAY", None)  # May already be unset.
            # Do not load any user options.
            # `os.environ` needs native strings on Py2+Windows.
            env[str("INKSCAPE_PROFILE_DIR")] = os.devnull
            # Old versions of Inkscape (0.48.3.1, used on Travis as of now)
            # seem to sometimes deadlock when stderr is redirected to a pipe,
            # so we redirect it to a temporary file instead.  This is not
            # necessary anymore as of Inkscape 0.92.1.
            self._stderr = TemporaryFile()
            self._proc = subprocess.Popen(
                [str("inkscape"), "--without-gui", "--shell"],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=self._stderr, env=env)
            if not self._read_to_prompt():
                raise OSError("Failed to start Inkscape")

        try:
            fsencode = os.fsencode
        except AttributeError:  # Py2.
            def fsencode(s):
                return s.encode(sys.getfilesystemencoding())

        # Inkscape uses glib's `g_shell_parse_argv`, which has a consistent
        # behavior across platforms, so we can just use `shlex.quote`.
        orig_b, dest_b = map(_shlex_quote_bytes, map(fsencode, [orig, dest]))
        if b"\n" in orig_b or b"\n" in dest_b:
            # Who knows whether the current folder name has a newline, or if
            # our encoding is even ASCII compatible...  Just fall back on the
            # slow solution (Inkscape uses `fgets` so it will always stop at a
            # newline).
            return make_external_conversion_command(lambda old, new: [
                str('inkscape'), '-z', old, '--export-png', new])(orig, dest)
        self._proc.stdin.write(orig_b + b" --export-png=" + dest_b + b"\n")
        self._proc.stdin.flush()
        if not self._read_to_prompt():
            # Inkscape's output is not localized but gtk's is, so the
            # output stream probably has a mixed encoding.  Using
            # `getfilesystemencoding` should at least get the filenames
            # right...
            self._stderr.seek(0)
            raise ImageComparisonFailure(
                self._stderr.read().decode(
                    sys.getfilesystemencoding(), "replace"))

    def __del__(self):
        if self._proc:
            if self._proc.poll() is None:  # Not exited yet.
                self._proc.communicate(b"quit\n")
                self._proc.wait()
            self._proc.stdin.close()
            self._proc.stdout.close()
            self._stderr.close()


def _update_converter():
    gs, gs_v = matplotlib.checkdep_ghostscript()
    if gs_v is not None:
        def cmd(old, new):
            return [str(gs), '-q', '-sDEVICE=png16m', '-dNOPAUSE', '-dBATCH',
             '-sOutputFile=' + new, old]
        converter['pdf'] = make_external_conversion_command(cmd)
        converter['eps'] = make_external_conversion_command(cmd)

    if matplotlib.checkdep_inkscape() is not None:
        converter['svg'] = _SVGConverter()


#: A dictionary that maps filename extensions to functions which
#: themselves map arguments `old` and `new` (filenames) to a list of strings.
#: The list can then be passed to Popen to convert files with that
#: extension to png format.
converter = {}
_update_converter()


def comparable_formats():
    """
    Returns the list of file formats that compare_images can compare
    on this system.

    """
    return ['png'] + list(converter)


def convert(filename, cache):
    """
    Convert the named file into a png file.  Returns the name of the
    created file.

    If *cache* is True, the result of the conversion is cached in
    `matplotlib._get_cachedir() + '/test_cache/'`.  The caching is based
    on a hash of the exact contents of the input file.  The is no limit
    on the size of the cache, so it may need to be manually cleared
    periodically.

    """
    base, extension = filename.rsplit('.', 1)
    if extension not in converter:
        reason = "Don't know how to convert %s files to png" % extension
        from . import is_called_from_pytest
        if is_called_from_pytest():
            import pytest
            pytest.skip(reason)
        else:
            from nose import SkipTest
            raise SkipTest(reason)
    newname = base + '_' + extension + '.png'
    if not os.path.exists(filename):
        raise IOError("'%s' does not exist" % filename)

    # Only convert the file if the destination doesn't already exist or
    # is out of date.
    if (not os.path.exists(newname) or
            os.stat(newname).st_mtime < os.stat(filename).st_mtime):
        if cache:
            cache_dir = get_cache_dir()
        else:
            cache_dir = None

        if cache_dir is not None:
            hash_value = get_file_hash(filename)
            new_ext = os.path.splitext(newname)[1]
            cached_file = os.path.join(cache_dir, hash_value + new_ext)
            if os.path.exists(cached_file):
                shutil.copyfile(cached_file, newname)
                return newname

        converter[extension](filename, newname)

        if cache_dir is not None:
            shutil.copyfile(newname, cached_file)

    return newname

#: Maps file extensions to a function which takes a filename as its
#: only argument to return a list suitable for execution with Popen.
#: The purpose of this is so that the result file (with the given
#: extension) can be verified with tools such as xmllint for svg.
verifiers = {}

# Turning this off, because it seems to cause multiprocessing issues
if False and matplotlib.checkdep_xmllint():
    verifiers['svg'] = lambda filename: [
        'xmllint', '--valid', '--nowarning', '--noout', filename]


@cbook.deprecated("2.1")
def verify(filename):
    """Verify the file through some sort of verification tool."""
    if not os.path.exists(filename):
        raise IOError("'%s' does not exist" % filename)
    base, extension = filename.rsplit('.', 1)
    verifier = verifiers.get(extension, None)
    if verifier is not None:
        cmd = verifier(filename)
        pipe = subprocess.Popen(cmd, universal_newlines=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate()
        errcode = pipe.wait()
        if errcode != 0:
            msg = "File verification command failed:\n%s\n" % ' '.join(cmd)
            if stdout:
                msg += "Standard output:\n%s\n" % stdout
            if stderr:
                msg += "Standard error:\n%s\n" % stderr
            raise IOError(msg)


def crop_to_same(actual_path, actual_image, expected_path, expected_image):
    # clip the images to the same size -- this is useful only when
    # comparing eps to pdf
    if actual_path[-7:-4] == 'eps' and expected_path[-7:-4] == 'pdf':
        aw, ah, ad = actual_image.shape
        ew, eh, ed = expected_image.shape
        actual_image = actual_image[int(aw / 2 - ew / 2):int(
            aw / 2 + ew / 2), int(ah / 2 - eh / 2):int(ah / 2 + eh / 2)]
    return actual_image, expected_image


def calculate_rms(expectedImage, actualImage):
    "Calculate the per-pixel errors, then compute the root mean square error."
    if expectedImage.shape != actualImage.shape:
        raise ImageComparisonFailure(
            "Image sizes do not match expected size: {0} "
            "actual size {1}".format(expectedImage.shape, actualImage.shape))
    # Convert to float to avoid overflowing finite integer types.
    return np.sqrt(((expectedImage - actualImage).astype(float) ** 2).mean())


def compare_images(expected, actual, tol, in_decorator=False):
    """
    Compare two "image" files checking differences within a tolerance.

    The two given filenames may point to files which are convertible to
    PNG via the `.converter` dictionary. The underlying RMS is calculated
    with the `.calculate_rms` function.

    Parameters
    ----------
    expected : str
        The filename of the expected image.
    actual :str
        The filename of the actual image.
    tol : float
        The tolerance (a color value difference, where 255 is the
        maximal difference).  The test fails if the average pixel
        difference is greater than this value.
    in_decorator : bool
        If called from image_comparison decorator, this should be
        True. (default=False)

    Examples
    --------
    img1 = "./baseline/plot.png"
    img2 = "./output/plot.png"
    compare_images( img1, img2, 0.001 ):

    """
    if not os.path.exists(actual):
        raise Exception("Output image %s does not exist." % actual)

    if os.stat(actual).st_size == 0:
        raise Exception("Output image file %s is empty." % actual)

    # Convert the image to png
    extension = expected.split('.')[-1]

    if not os.path.exists(expected):
        raise IOError('Baseline image %r does not exist.' % expected)

    if extension != 'png':
        actual = convert(actual, False)
        expected = convert(expected, True)

    # open the image files and remove the alpha channel (if it exists)
    expectedImage = _png.read_png_int(expected)
    actualImage = _png.read_png_int(actual)
    expectedImage = expectedImage[:, :, :3]
    actualImage = actualImage[:, :, :3]

    actualImage, expectedImage = crop_to_same(
        actual, actualImage, expected, expectedImage)

    diff_image = make_test_filename(actual, 'failed-diff')

    if tol <= 0.0:
        if np.array_equal(expectedImage, actualImage):
            return None

    # convert to signed integers, so that the images can be subtracted without
    # overflow
    expectedImage = expectedImage.astype(np.int16)
    actualImage = actualImage.astype(np.int16)

    rms = calculate_rms(expectedImage, actualImage)

    if rms <= tol:
        return None

    save_diff_image(expected, actual, diff_image)

    results = dict(rms=rms, expected=str(expected),
                   actual=str(actual), diff=str(diff_image), tol=tol)

    if not in_decorator:
        # Then the results should be a string suitable for stdout.
        template = ['Error: Image files did not match.',
                    'RMS Value: {rms}',
                    'Expected:  \n    {expected}',
                    'Actual:    \n    {actual}',
                    'Difference:\n    {diff}',
                    'Tolerance: \n    {tol}', ]
        results = '\n  '.join([line.format(**results) for line in template])
    return results


def save_diff_image(expected, actual, output):
    expectedImage = _png.read_png(expected)
    actualImage = _png.read_png(actual)
    actualImage, expectedImage = crop_to_same(
        actual, actualImage, expected, expectedImage)
    expectedImage = np.array(expectedImage).astype(float)
    actualImage = np.array(actualImage).astype(float)
    if expectedImage.shape != actualImage.shape:
        raise ImageComparisonFailure(
            "Image sizes do not match expected size: {0} "
            "actual size {1}".format(expectedImage.shape, actualImage.shape))
    absDiffImage = np.abs(expectedImage - actualImage)

    # expand differences in luminance domain
    absDiffImage *= 255 * 10
    save_image_np = np.clip(absDiffImage, 0, 255).astype(np.uint8)
    height, width, depth = save_image_np.shape

    # The PDF renderer doesn't produce an alpha channel, but the
    # matplotlib PNG writer requires one, so expand the array
    if depth == 3:
        with_alpha = np.empty((height, width, 4), dtype=np.uint8)
        with_alpha[:, :, 0:3] = save_image_np
        save_image_np = with_alpha

    # Hard-code the alpha channel to fully solid
    save_image_np[:, :, 3] = 255

    _png.write_png(save_image_np, output)
