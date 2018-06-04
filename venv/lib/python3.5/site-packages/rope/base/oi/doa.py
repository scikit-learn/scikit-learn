try:
    import pickle
except ImportError:
    import cPickle as pickle
import marshal
import os
import socket
import subprocess
import sys
import tempfile
import threading


class PythonFileRunner(object):
    """A class for running python project files"""

    def __init__(self, pycore, file_, args=None, stdin=None,
                 stdout=None, analyze_data=None):
        self.pycore = pycore
        self.file = file_
        self.analyze_data = analyze_data
        self.observers = []
        self.args = args
        self.stdin = stdin
        self.stdout = stdout

    def run(self):
        """Execute the process"""
        env = dict(os.environ)
        file_path = self.file.real_path
        path_folders = self.pycore.project.get_source_folders() + \
            self.pycore.project.get_python_path_folders()
        env['PYTHONPATH'] = os.pathsep.join(folder.real_path
                                            for folder in path_folders)
        runmod_path = self.pycore.project.find_module('rope.base.oi.runmod').real_path
        self.receiver = None
        self._init_data_receiving()
        send_info = '-'
        if self.receiver:
            send_info = self.receiver.get_send_info()
        args = [sys.executable, runmod_path, send_info,
                self.pycore.project.address, self.file.real_path]
        if self.analyze_data is None:
            del args[1:4]
        if self.args is not None:
            args.extend(self.args)
        self.process = subprocess.Popen(
            executable=sys.executable, args=args, env=env,
            cwd=os.path.split(file_path)[0], stdin=self.stdin,
            stdout=self.stdout, stderr=self.stdout, close_fds=os.name != 'nt')

    def _init_data_receiving(self):
        if self.analyze_data is None:
            return
        # Disabling FIFO data transfer due to blocking when running
        # unittests in the GUI.
        # XXX: Handle FIFO data transfer for `rope.ui.testview`
        if True or os.name == 'nt':
            self.receiver = _SocketReceiver()
        else:
            self.receiver = _FIFOReceiver()
        self.receiving_thread = threading.Thread(
            target=self._receive_information)
        self.receiving_thread.setDaemon(True)
        self.receiving_thread.start()

    def _receive_information(self):
        #temp = open('/dev/shm/info', 'wb')
        for data in self.receiver.receive_data():
            self.analyze_data(data)
            #temp.write(str(data) + '\n')
        #temp.close()
        for observer in self.observers:
            observer()

    def wait_process(self):
        """Wait for the process to finish"""
        self.process.wait()
        if self.analyze_data:
            self.receiving_thread.join()

    def kill_process(self):
        """Stop the process"""
        if self.process.poll() is not None:
            return
        try:
            if hasattr(self.process, 'terminate'):
                self.process.terminate()
            elif os.name != 'nt':
                os.kill(self.process.pid, 9)
            else:
                import ctypes
                handle = int(self.process._handle)
                ctypes.windll.kernel32.TerminateProcess(handle, -1)
        except OSError:
            pass

    def add_finishing_observer(self, observer):
        """Notify this observer when execution finishes"""
        self.observers.append(observer)


class _MessageReceiver(object):

    def receive_data(self):
        pass

    def get_send_info(self):
        pass


class _SocketReceiver(_MessageReceiver):

    def __init__(self):
        self.server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.data_port = 3037
        while self.data_port < 4000:
            try:
                self.server_socket.bind(('', self.data_port))
                break
            except socket.error:
                self.data_port += 1
        self.server_socket.listen(1)

    def get_send_info(self):
        return str(self.data_port)

    def receive_data(self):
        conn, addr = self.server_socket.accept()
        self.server_socket.close()
        my_file = conn.makefile('rb')
        while True:
            try:
                yield pickle.load(my_file)
            except EOFError:
                break
        my_file.close()
        conn.close()


class _FIFOReceiver(_MessageReceiver):

    def __init__(self):
        # XXX: this is insecure and might cause race conditions
        self.file_name = self._get_file_name()
        os.mkfifo(self.file_name)

    def _get_file_name(self):
        prefix = tempfile.gettempdir() + '/__rope_'
        i = 0
        while os.path.exists(prefix + str(i).rjust(4, '0')):
            i += 1
        return prefix + str(i).rjust(4, '0')

    def get_send_info(self):
        return self.file_name

    def receive_data(self):
        my_file = open(self.file_name, 'rb')
        while True:
            try:
                yield marshal.load(my_file)
            except EOFError:
                break
        my_file.close()
        os.remove(self.file_name)
