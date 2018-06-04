# -*- encoding: utf-8 -*-
"""
django-videothumbs
"""

from django.utils.six import StringIO
from django.conf.global_settings import FILE_UPLOAD_TEMP_DIR
import math
import os
import time
from PIL import Image

TMP_DIR = FILE_UPLOAD_TEMP_DIR or '/tmp/'

def generate_thumb(storage, video_path, thumb_size=None, format='jpg', frames=100):
    histogram = []

    http_status = 200
    
    name = video_path
    path = storage.path(video_path)

    if not storage.exists(video_path):
        return "", '404'

    framemask = "%s%s%s%s" % (TMP_DIR,
                              name.split('/')[-1].split('.')[0] + str(time.time()),
                              '.%d.',
                              format)
    # ffmpeg command for grabbing N number of frames
    cmd = "/usr/bin/ffmpeg -y -t 00:00:05 -i '%s' '%s'" % (path, framemask)

    # make sure that this command worked or return.
    if os.system(cmd) != 0:
        return "", '500'

    # loop through the generated images, open, and generate the image histogram.
    for i in range(1, frames + 1):
        fname = framemask % i

        if not os.path.exists(os.path.join(TMP_DIR, fname)):
            break

        image = Image.open(fname)

        # Convert to RGB if necessary
        if image.mode not in ('L', 'RGB'):
            image = image.convert('RGB')

        # The list of image historgrams
        histogram.append(image.histogram())

    n = len(histogram)
    avg = []

    # Get the image average of the first image
    for c in range(len(histogram[0])):
        ac = 0.0
        for i in range(n):
            ac = ac + (float(histogram[i][c]) / n)
        avg.append(ac)

    minn = -1
    minRMSE = -1

    # Compute the mean squared error
    for i in range(1, n + 1):
        results = 0.0
        num = len(avg)

        for j in range(num):
            median_error = avg[j] - float(histogram[i - 1][j])
            results += (median_error * median_error) / num
        rmse = math.sqrt(results)

        if minn == -1 or rmse < minRMSE:
            minn = i
            minRMSE = rmse

    file_location = framemask % (minn)
    image = Image.open(file_location)

    # If you want to generate a square thumbnail
    if not thumb_size is None:
        thumb_w, thumb_h = thumb_size

        if thumb_w == thumb_h:
            # quad
            xsize, ysize = image.size
            # get minimum size
            minsize = min(xsize, ysize)
            # largest square possible in the image
            xnewsize = (xsize - minsize) / 2
            ynewsize = (ysize - minsize) / 2
            # crop it
            image2 = image.crop((xnewsize, ynewsize, xsize - xnewsize, ysize - ynewsize))
            # load is necessary after crop
            image2.load()
            # thumbnail of the cropped image (with ANTIALIAS to make it look better)
            image2.thumbnail(thumb_size, Image.ANTIALIAS)
        else:
            # not quad
            image2 = image
            image2.thumbnail(thumb_size, Image.ANTIALIAS)
    else:
        image2 = image

    io = StringIO()

    # PNG and GIF are the same, JPG is JPEG
    if format.upper() == 'JPG':
        format = 'JPEG'

    image2.save(io, format)

    # We don't know how many frames we capture. We just captured the first 5 seconds, so keep removing until not found
    for i in range(1, 9999):
        fname = framemask % i
        try:
            os.unlink(fname)
        except OSError:
            break

    return io.getvalue(), http_status
