# coding: utf-8
"""Tests for the IPython tab-completion machinery.
"""
#-----------------------------------------------------------------------------
# Module imports
#-----------------------------------------------------------------------------

# stdlib
import io
import os
import sys
import tempfile
from datetime import datetime

# third party
import nose.tools as nt

# our own packages
from traitlets.config.loader import Config
from IPython.utils.tempdir import TemporaryDirectory
from IPython.core.history import HistoryManager, extract_hist_ranges

def setUp():
    nt.assert_equal(sys.getdefaultencoding(), "utf-8")

def test_history():
    ip = get_ipython()
    with TemporaryDirectory() as tmpdir:
        hist_manager_ori = ip.history_manager
        hist_file = os.path.join(tmpdir, 'history.sqlite')
        try:
            ip.history_manager = HistoryManager(shell=ip, hist_file=hist_file)
            hist = [u'a=1', u'def f():\n    test = 1\n    return test', u"b='€Æ¾÷ß'"]
            for i, h in enumerate(hist, start=1):
                ip.history_manager.store_inputs(i, h)

            ip.history_manager.db_log_output = True
            # Doesn't match the input, but we'll just check it's stored.
            ip.history_manager.output_hist_reprs[3] = "spam"
            ip.history_manager.store_output(3)

            nt.assert_equal(ip.history_manager.input_hist_raw, [''] + hist)
            
            # Detailed tests for _get_range_session
            grs = ip.history_manager._get_range_session
            nt.assert_equal(list(grs(start=2,stop=-1)), list(zip([0], [2], hist[1:-1])))
            nt.assert_equal(list(grs(start=-2)), list(zip([0,0], [2,3], hist[-2:])))
            nt.assert_equal(list(grs(output=True)), list(zip([0,0,0], [1,2,3], zip(hist, [None,None,'spam']))))

            # Check whether specifying a range beyond the end of the current
            # session results in an error (gh-804)
            ip.magic('%hist 2-500')
            
            # Check that we can write non-ascii characters to a file
            ip.magic("%%hist -f %s" % os.path.join(tmpdir, "test1"))
            ip.magic("%%hist -pf %s" % os.path.join(tmpdir, "test2"))
            ip.magic("%%hist -nf %s" % os.path.join(tmpdir, "test3"))
            ip.magic("%%save %s 1-10" % os.path.join(tmpdir, "test4"))

            # New session
            ip.history_manager.reset()
            newcmds = [u"z=5",
                       u"class X(object):\n    pass",
                       u"k='p'",
                       u"z=5"]
            for i, cmd in enumerate(newcmds, start=1):
                ip.history_manager.store_inputs(i, cmd)
            gothist = ip.history_manager.get_range(start=1, stop=4)
            nt.assert_equal(list(gothist), list(zip([0,0,0],[1,2,3], newcmds)))
            # Previous session:
            gothist = ip.history_manager.get_range(-1, 1, 4)
            nt.assert_equal(list(gothist), list(zip([1,1,1],[1,2,3], hist)))

            newhist = [(2, i, c) for (i, c) in enumerate(newcmds, 1)]

            # Check get_hist_tail
            gothist = ip.history_manager.get_tail(5, output=True,
                                                    include_latest=True)
            expected = [(1, 3, (hist[-1], "spam"))] \
                + [(s, n, (c, None)) for (s, n, c) in newhist]
            nt.assert_equal(list(gothist), expected)

            gothist = ip.history_manager.get_tail(2)
            expected = newhist[-3:-1]
            nt.assert_equal(list(gothist), expected)

            # Check get_hist_search
            gothist = ip.history_manager.search("*test*")
            nt.assert_equal(list(gothist), [(1,2,hist[1])] )

            gothist = ip.history_manager.search("*=*")
            nt.assert_equal(list(gothist),
                            [(1, 1, hist[0]),
                             (1, 2, hist[1]),
                             (1, 3, hist[2]),
                             newhist[0],
                             newhist[2],
                             newhist[3]])

            gothist = ip.history_manager.search("*=*", n=4)
            nt.assert_equal(list(gothist),
                            [(1, 3, hist[2]),
                             newhist[0],
                             newhist[2],
                             newhist[3]])

            gothist = ip.history_manager.search("*=*", unique=True)
            nt.assert_equal(list(gothist),
                            [(1, 1, hist[0]),
                             (1, 2, hist[1]),
                             (1, 3, hist[2]),
                             newhist[2],
                             newhist[3]])

            gothist = ip.history_manager.search("*=*", unique=True, n=3)
            nt.assert_equal(list(gothist),
                            [(1, 3, hist[2]),
                             newhist[2],
                             newhist[3]])

            gothist = ip.history_manager.search("b*", output=True)
            nt.assert_equal(list(gothist), [(1,3,(hist[2],"spam"))] )

            # Cross testing: check that magic %save can get previous session.
            testfilename = os.path.realpath(os.path.join(tmpdir, "test.py"))
            ip.magic("save " + testfilename + " ~1/1-3")
            with io.open(testfilename, encoding='utf-8') as testfile:
                nt.assert_equal(testfile.read(),
                                        u"# coding: utf-8\n" + u"\n".join(hist)+u"\n")

            # Duplicate line numbers - check that it doesn't crash, and
            # gets a new session
            ip.history_manager.store_inputs(1, "rogue")
            ip.history_manager.writeout_cache()
            nt.assert_equal(ip.history_manager.session_number, 3)
        finally:
            # Ensure saving thread is shut down before we try to clean up the files
            ip.history_manager.save_thread.stop()
            # Forcibly close database rather than relying on garbage collection
            ip.history_manager.db.close()
            # Restore history manager
            ip.history_manager = hist_manager_ori


def test_extract_hist_ranges():
    instr = "1 2/3 ~4/5-6 ~4/7-~4/9 ~9/2-~7/5 ~10/"
    expected = [(0, 1, 2),  # 0 == current session
                (2, 3, 4),
                (-4, 5, 7),
                (-4, 7, 10),
                (-9, 2, None),  # None == to end
                (-8, 1, None),
                (-7, 1, 6),
                (-10, 1, None)]
    actual = list(extract_hist_ranges(instr))
    nt.assert_equal(actual, expected)

def test_magic_rerun():
    """Simple test for %rerun (no args -> rerun last line)"""
    ip = get_ipython()
    ip.run_cell("a = 10", store_history=True)
    ip.run_cell("a += 1", store_history=True)
    nt.assert_equal(ip.user_ns["a"], 11)
    ip.run_cell("%rerun", store_history=True)
    nt.assert_equal(ip.user_ns["a"], 12)

def test_timestamp_type():
    ip = get_ipython()
    info = ip.history_manager.get_session_info()
    nt.assert_true(isinstance(info[1], datetime))

def test_hist_file_config():
    cfg = Config()
    tfile = tempfile.NamedTemporaryFile(delete=False)
    cfg.HistoryManager.hist_file = tfile.name
    try:
        hm = HistoryManager(shell=get_ipython(), config=cfg)
        nt.assert_equal(hm.hist_file, cfg.HistoryManager.hist_file)
    finally:
        try:
            os.remove(tfile.name)
        except OSError:
            # same catch as in testing.tools.TempFileMixin
            # On Windows, even though we close the file, we still can't
            # delete it.  I have no clue why
            pass

def test_histmanager_disabled():
    """Ensure that disabling the history manager doesn't create a database."""
    cfg = Config()
    cfg.HistoryAccessor.enabled = False

    ip = get_ipython()
    with TemporaryDirectory() as tmpdir:
        hist_manager_ori = ip.history_manager
        hist_file = os.path.join(tmpdir, 'history.sqlite')
        cfg.HistoryManager.hist_file = hist_file
        try:
            ip.history_manager = HistoryManager(shell=ip, config=cfg)
            hist = [u'a=1', u'def f():\n    test = 1\n    return test', u"b='€Æ¾÷ß'"]
            for i, h in enumerate(hist, start=1):
                ip.history_manager.store_inputs(i, h)
            nt.assert_equal(ip.history_manager.input_hist_raw, [''] + hist)
            ip.history_manager.reset()
            ip.history_manager.end_session()
        finally:
            ip.history_manager = hist_manager_ori

    # hist_file should not be created
    nt.assert_false(os.path.exists(hist_file))
