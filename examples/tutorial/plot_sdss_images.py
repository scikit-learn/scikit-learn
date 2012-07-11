"""
SDSS Images
-----------

This script plots an example quasar, star, and galaxy image for use in
the tutorial.
"""
import os
import urllib2

import pylab as pl
from matplotlib import image

from sklearn.datasets import get_data_home


def _fetch(outfile, RA, DEC, scale=0.2, width=400, height=400):
    """Fetch the image at the given RA, DEC from the SDSS server"""
    url = ("http://casjobs.sdss.org/ImgCutoutDR7/"
           "getjpeg.aspx?ra=%.8f&dec=%.8f&scale=%.2f&width=%i&height=%i"
           % (RA, DEC, scale, width, height))
    print "downloading %s" % url
    print " -> %s" % outfile
    fhandle = urllib2.urlopen(url)
    open(outfile, 'w').write(fhandle.read())


def fetch_image(object_type, data_home=None):
    """Return the data array for the image of object type"""
    data_home = get_data_home(data_home)
    if not os.path.exists(data_home):
        os.makedirs(data_home)

    filename = os.path.join(data_home, '%s_image.jpg' % object_type)
    if not os.path.exists(filename):
        RA = image_locations[object_type]['RA']
        DEC = image_locations[object_type]['DEC']
        _fetch(filename, RA, DEC)
    
    return image.imread(filename)


image_locations = dict(star=dict(RA=180.63040108,
                                 DEC=64.96767375),
                       galaxy=dict(RA=197.51943983,
                                   DEC=0.94881436),
                       quasar=dict(RA=226.18451462,
                                   DEC=4.07456639))


# Plot the images
fig = pl.figure(figsize=(9, 3))

# Check that PIL is installed for jpg support
if 'jpg' not in fig.canvas.get_supported_filetypes():
    raise ValueError("PIL required to load SDSS jpeg images")

object_types = ['star', 'galaxy', 'quasar']

for i, object_type in enumerate(object_types):
    ax = pl.subplot(131 + i, xticks=[], yticks=[])
    I = fetch_image(object_type)
    ax.imshow(I)
    if object_type != 'galaxy':
        pl.arrow(0.65, 0.65, -0.1, -0.1, width=0.005, head_width=0.03,
                 length_includes_head=True,
                 color='w', transform=ax.transAxes)
    pl.text(0.99, 0.01, object_type, fontsize='large', color='w', ha='right',
            transform=ax.transAxes)

pl.subplots_adjust(bottom=0.04, top=0.94, left=0.02, right=0.98, wspace=0.04)

pl.show()
