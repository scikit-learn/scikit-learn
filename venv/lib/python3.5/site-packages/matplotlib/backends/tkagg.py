from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import tkinter as Tk

import numpy as np

from matplotlib.backends import _tkagg

def blit(photoimage, aggimage, bbox=None, colormode=1):
    tk = photoimage.tk

    if bbox is not None:
        bbox_array = bbox.__array__()
        # x1, x2, y1, y2
        bboxptr = (bbox_array[0, 0], bbox_array[1, 0],
                   bbox_array[0, 1], bbox_array[1, 1])
    else:
        bboxptr = 0
    data = np.asarray(aggimage)
    dataptr = (data.shape[0], data.shape[1], data.ctypes.data)
    try:
        tk.call(
            "PyAggImagePhoto", photoimage,
            dataptr, colormode, bboxptr)
    except Tk.TclError:
        if hasattr(tk, 'interpaddr'):
            _tkagg.tkinit(tk.interpaddr(), 1)
        else:
            # very old python?
            _tkagg.tkinit(tk, 0)
        tk.call("PyAggImagePhoto", photoimage,
                dataptr, colormode, bboxptr)

def test(aggimage):
    r = Tk.Tk()
    c = Tk.Canvas(r, width=aggimage.width, height=aggimage.height)
    c.pack()
    p = Tk.PhotoImage(width=aggimage.width, height=aggimage.height)
    blit(p, aggimage)
    c.create_image(aggimage.width,aggimage.height,image=p)
    blit(p, aggimage)
    while True: r.update_idletasks()
