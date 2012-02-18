import numpy as np
from scipy import stats

e = np.exp(1)

np.random.seed(4)


try:
    from enthought.mayavi import mlab
except ImportError:
    from mayavi import mlab


y = np.random.normal(scale=0.5, size=(200000))
x = np.random.normal(scale=0.5, size=(200000))
z = np.random.normal(scale=0.1, size=(len(x)), )
def pdf(x):
    return 0.5*(  stats.norm(scale=0.25/e).pdf(x) 
                + stats.norm(scale=4/e).pdf(x))
density = pdf(x) * pdf(y)
pdf_z = pdf(5*z)

density *= pdf_z

a = x+y
b = 2*y
c = a-b+z

norm = np.sqrt(a.var() + b.var())
a /= norm
b /= norm

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
mlab.clf()
pts = mlab.points3d(a, b, c, density, colormap='jet', mode='2dvertex')
mlab.outline(extent=[-3*a.std(), 3*a.std(), -3*b.std(), 3*b.std(),
                     -3*c.std(), 3*c.std()], line_width=2)

Y = np.c_[a, b, c]
U, pca_score, V = np.linalg.svd(Y, full_matrices=False)
x_pca_axis, y_pca_axis, z_pca_axis = V.T*pca_score/pca_score.min()

mlab.view(-20.8, 83, 9, [0.18, 0.2, -0.24])
#mlab.savefig('pca_3d.jpg')
mlab.quiver3d(0.1*x_pca_axis, 0.1*y_pca_axis, 0.1*z_pca_axis,
                2*x_pca_axis, 2*y_pca_axis, 2*z_pca_axis,
                color=(0.6, 0, 0), line_width=2)

x_pca_axis, y_pca_axis, z_pca_axis = 3*V.T
x_pca_plane = np.r_[x_pca_axis[:2], - x_pca_axis[1::-1]]
y_pca_plane = np.r_[y_pca_axis[:2], - y_pca_axis[1::-1]]
z_pca_plane = np.r_[z_pca_axis[:2], - z_pca_axis[1::-1]]

x_pca_plane.shape = (2, 2)
y_pca_plane.shape = (2, 2)
z_pca_plane.shape = (2, 2)

mlab.mesh(x_pca_plane, y_pca_plane, z_pca_plane, color=(0.6, 0, 0),
            opacity=0.1)
mlab.mesh(x_pca_plane, y_pca_plane, z_pca_plane, color=(0.6, 0, 0),
            representation='wireframe', line_width=3, opacity=0.3)

#mlab.title('PCA axis')
mlab.view(-20.8, 83, 9, [0.18, 0.2, -0.24])
#mlab.savefig('pca_3d_axis.jpg')

# A view
mlab.view(3.3, 43.8, 9.2, [0.04, -0.11, -0.17])
#mlab.savefig('pca_3d_aligned.jpg')

