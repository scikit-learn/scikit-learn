import threading
import time

import lazy_loader as lazy


def import_np():
    time.sleep(0.5)
    lazy.load("numpy")


for _ in range(10):
    threading.Thread(target=import_np).start()
