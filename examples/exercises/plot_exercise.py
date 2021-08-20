import matplotlib.pyplot as plt
import numpy as np

plt.style.use('default')
plt.style.use('seaborn-talk')

# class I
X01 = [(2, 9), (7, 19), (1, 10), (3, 19), (4, 16), (5, 18)]

# class II
X02 = [(4, 3), (6, 7), (1, -10), (3, -1), (9, 5), (5, -7)]

plt.xlim([0, 10])
plt.scatter(*(zip(*X01)))
plt.scatter(*(zip(*X02)))
plt.show()

plt.scatter(*(zip(*X01)))
plt.scatter(*(zip(*X02)))
plt.plot([0, 10], [0, 15], 'C2-');
plt.show()

plt.scatter(*(zip(*X01)))
plt.scatter(*(zip(*X02)))
plt.scatter(5, 9, c='C3')  # test point
plt.plot([0, 10], [2, 22], 'C4-');
plt.plot([0, 10], [0, 15], 'C2-');
plt.show()

plt.plot([0, 10], [0, 20], 'C2-', zorder=0)
plt.plot([0, 10], [5, 25], 'C2--', zorder=0)
plt.plot([0, 10], [-5, 15], 'C2--', zorder=0)
plt.scatter(5, 9, c='C3')  # test point
plt.annotate('', (0,5), (0.55, -4), arrowprops=dict(arrowstyle='<->'))
plt.text(0.5, -1.5, '$m$')
plt.text(0.3, 2.5, '$m$')
plt.annotate('', (3.05,18.1), (3.8, 7.75), arrowprops=dict(arrowstyle='<->'))
plt.text(2.9, 14, '$m_i$')
plt.scatter(*(zip(*X01)), zorder=1)
plt.scatter(*(zip(*X02)), zorder=1)
sv = X01[:2] + X02[:2]
plt.scatter(*(zip(*sv)), zorder=1, facecolors='none', edgecolors='r', s=500);
plt.show()

def flc(c, n=100):
  """Level curves of objective functions"""
  return np.array([(c * np.cos(ang), c * np.sin(ang))
                   for ang in np.linspace(0, 2*np.pi, n)])
def g(n=100):
  """Constraint"""
  return np.array([(x, 1./x) for x in np.linspace(0.1, 2.5, n)])
plt.figure(figsize=(8, 8))
plt.xlim([-2.5, 2.5])
plt.ylim([-2.5, 2.5])
plt.grid(color='0.5', linestyle='--', linewidth=0.5)
for c in (0.2, 0.6, 1, 1.8, 2.2):
  plt.plot(*flc(c).T, color='C0')
plt.plot(*flc(np.sqrt(2)).T, color='C0', label='$f(x, y) = c$')
plt.plot(*g().T, color='C1', label='$g(x, y) = 0$')
plt.plot(*-g().T, color='C1')
plt.scatter(1, 1, c='C2', zorder=4)
plt.scatter(np.sqrt(0.345), np.sqrt(1/0.345), c='C3', zorder=4)
plt.legend();
plt.show()

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
ax = plt.figure().add_subplot(111, projection=Axes3D.name)
X, Y = np.meshgrid(np.linspace(-1, 3, 50), np.linspace(-5, 1, 50))
L = X**2 + Y * (X - 1)
ax.plot_surface(X, Y, L, cmap=cm.hsv)
ax.view_init(elev=45, azim=120)
ax.set_xlabel('$x$', labelpad=20)
ax.set_ylabel('$\lambda$', labelpad=20)
ax.set_zlabel('$\mathcal{L}$', labelpad=10);
plt.show()

def generate_circle_data(R1=0, R2=1, N=500):
  """Generate N points in a circle for radius range (R1, R2)"""
  r = lambda: R1 + np.random.random() * (R2 - R1)
  return np.array([(r() * np.cos(ang), r() * np.sin(ang))
                   for ang in np.linspace(0, 2*np.pi, N)])
C01 = generate_circle_data()
C02 = generate_circle_data(1, 2)
plt.scatter(*C01.T, marker='.')
plt.scatter(*C02.T, marker='.');
plt.show()


ax = plt.figure().add_subplot(111, projection=Axes3D.name)
Z01 = np.array([x**2 + y**2 for x, y in C01])
Z02 = np.array([x**2 + y**2 for x, y in C02])
ax.scatter(*C01.T, Z01, cmap=cm.hsv)
ax.scatter(*C02.T, Z02, cmap=cm.hsv)
ax.view_init(elev=15, azim=60)
plt.show()

from sklearn import svm
clf = svm.SVC(kernel="linear")
clf_linear = svm.SVC(kernel="linear")
clf_rbf = svm.SVC(kernel="rbf")
clf_poly3 = svm.SVC(kernel="poly", degree=3)
clf_poly10 = svm.SVC(kernel="poly", degree=10)
titles = ("Linear", "RBF", "Polynomial, degree = 3", "Plynomial, degree = 10")
xx, yy = np.meshgrid(np.arange(-3, 3, 0.01), np.arange(-3, 3, 0.01))
for i, clf in enumerate((clf_linear, clf_rbf, clf_poly3, clf_poly10)):
  clf.fit(np.concatenate((C01, C02), axis=0), [-1]*len(C01) + [1]*len(C02))
  plt.subplot(2, 2, i + 1)
  Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
  Z = Z.reshape(xx.shape)
  plt.contourf(xx, yy, Z, cmap=plt.cm.Paired)
  plt.scatter(*C01.T, marker='.')
  plt.scatter(*C02.T, marker='.')
  plt.title(titles[i])
plt.tight_layout()
plt.show()
