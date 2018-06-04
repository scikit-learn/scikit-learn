import sys
from IPython.external.qt_for_kernel import QtCore, QtGui

# If we create a QApplication, keep a reference to it so that it doesn't get
# garbage collected.
_appref = None


def inputhook(context):
    global _appref
    app = QtCore.QCoreApplication.instance()
    if not app:
        _appref = app = QtGui.QApplication([" "])
    event_loop = QtCore.QEventLoop(app)

    if sys.platform == 'win32':
        # The QSocketNotifier method doesn't appear to work on Windows.
        # Use polling instead.
        timer = QtCore.QTimer()
        timer.timeout.connect(event_loop.quit)
        while not context.input_is_ready():
            timer.start(50)  # 50 ms
            event_loop.exec_()
            timer.stop()
    else:
        # On POSIX platforms, we can use a file descriptor to quit the event
        # loop when there is input ready to read.
        notifier = QtCore.QSocketNotifier(context.fileno(),
                                          QtCore.QSocketNotifier.Read)
        # connect the callback we care about before we turn it on
        notifier.activated.connect(event_loop.exit)
        notifier.setEnabled(True)
        # only start the event loop we are not already flipped
        if not context.input_is_ready():
            event_loop.exec_()
