import numpy as np
import matplotlib.pyplot as plt
from scikits.learn.machine.manifold_learning.regression.neighbors import \
     Neighbors, Kneighbors, Parzen


n = 100 # number of points
data1 = np.random.randn(n,2) + 3.0
data2 = np.random.randn(n, 2) + 5.0
data = np.concatenate((data1, data2))
#plt.subplot(211)

ax = plt.gca()

#plt.subplot(212)
# mesh
h = 0.1 # step size
x = np.arange(0, 12, h)
y = np.arange(0, 12, h)
X, Y = np.meshgrid(x, y)
Z = np.zeros(X.shape)


neigh = Neighbors(data, 3)


for i in range(len(x)):
    for j in range(len(y)):
        nn = np.array(neigh.kneighbors(np.array((x[i], y[j]))))
        _t = np.sign(nn[:,1]-100)
        Z[i,j] = np.sign(_t.sum())

ax = plt.subplot(111)
plt.pcolormesh(X, Y, Z)
## im = plt.imshow(Z, cmap=plt.cm.jet)
#im.set_interpolation('bilinear')
plt.scatter(data1[:,0], data1[:,1])
plt.scatter(data2[:,0], data2[:,1], c='red')
#plt.show()

plt.show()
