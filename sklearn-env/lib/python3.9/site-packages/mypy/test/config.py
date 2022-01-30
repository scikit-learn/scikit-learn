import os.path

provided_prefix = os.getenv('MYPY_TEST_PREFIX', None)
if provided_prefix:
    PREFIX = provided_prefix
else:
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    PREFIX = os.path.dirname(os.path.dirname(this_file_dir))

# Location of test data files such as test case descriptions.
test_data_prefix = os.path.join(PREFIX, 'test-data', 'unit')
package_path = os.path.join(PREFIX, 'test-data', 'packages')

# Temp directory used for the temp files created when running test cases.
# This is *within* the tempfile.TemporaryDirectory that is chroot'ed per testcase.
# It is also hard-coded in numerous places, so don't change it.
test_temp_dir = 'tmp'
