def generate_toy_data(n_atoms, n_samples, image_size):
    n_features = image_size[0] * image_size[1]

    np.random.seed(0)
    U = np.random.randn(n_samples, n_atoms)
    V = np.random.randn(n_atoms, n_features)

    centers = [(3, 3), (6, 7), (8, 1)]
    sz = [1, 2, 1]
    for k in range(n_atoms):
        img = np.zeros(image_size)
        xmin, xmax = centers[k][0] - sz[k], centers[k][0] + sz[k]
        ymin, ymax = centers[k][1] - sz[k], centers[k][1] + sz[k]
        img[xmin:xmax][:, ymin:ymax] = 1.0
        V[k, :] = img.ravel()

    # Y is defined by : Y = UV + noise
    Y = np.dot(U, V)
    Y += 0.1 * np.random.randn(Y.shape[0], Y.shape[1])  # Add noise
    return Y, U, V


# Generate toy data
n_atoms = 3
n_samples = 100
img_sz = (10, 10)
Y, U, V = generate_toy_data(n_atoms, n_samples, img_sz)

# Estimate U,V
alpha = 0.5
SPCA = SparsePCA(n_atoms, alpha, max_iter=100, method='lasso', n_jobs=1)
SPCA.fit(Y)
V_estimated = SPCA.components_

# View results
import pylab as pl
pl.close('all')

for k in range(n_atoms):
    pl.matshow(np.reshape(V_estimated[k, :], img_sz))
    pl.title('Atom %d' % k)
    pl.colorbar()

pl.figure()
pl.plot(SPCA.error_)
pl.xlabel('Iteration')
pl.ylabel('Cost function')
pl.show()