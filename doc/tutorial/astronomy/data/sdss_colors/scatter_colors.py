import numpy as np
import pylab as pl

data = np.load('sdssdr6_colors_class_train.npy')

# only plot 10000 points: otherwise it takes too much memory
np.random.shuffle(data)
data = data[:10000]

redshift = data['redshift']

print "%i qsos" % np.sum(redshift > 0)
print "%i stars" % np.sum(redshift == 0)

kwargs = dict(s=1, c=(redshift > 0), lw=0)

pl.figure(figsize=(6, 8))

pl.subplot(311).scatter(data['u-g'], data['g-r'], **kwargs)

pl.subplot(312).scatter(data['g-r'], data['r-i'], **kwargs)

pl.subplot(313).scatter(data['r-i'], data['i-z'], **kwargs)

pl.show()
