import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils.multiclass import unique_labels

from sklearn.glvq.gmlvq import GmlvqModel
from sklearn.glvq.grlvq import GrlvqModel
from sklearn.glvq.lgmlvq import LgmlvqModel

from sklearn.glvq.glvq import GlvqModel

nb_ppc = 100
toy_data = np.append(np.random.multivariate_normal([0, 0], np.eye(2)/2, size=nb_ppc),
                     np.random.multivariate_normal([5, 0], np.eye(2)/2, size=nb_ppc), axis=0)
toy_label = np.append(np.zeros(nb_ppc), np.ones(nb_ppc), axis=0)

glvq = GlvqModel()
glvq.fit(toy_data,toy_label)
pred = glvq.predict(toy_data)

plt.subplot(421)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(glvq.w_[:, 0], glvq.w_[:, 1])
plt.axis('equal')

plt.subplot(422)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(glvq.w_[:, 0], glvq.w_[:, 1])
plt.axis('equal')

toy_data = np.append(np.random.multivariate_normal([0, 0], np.array([[0.3,0],[0,4]]), size=nb_ppc),
                     np.random.multivariate_normal([4, 4], np.array([[0.3,0],[0,4]]), size=nb_ppc), axis=0)
grlvq = GrlvqModel()
grlvq.fit(toy_data,toy_label)
pred = grlvq.predict(toy_data)

plt.subplot(423)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(grlvq.w_[:, 0], grlvq.w_[:, 1])
plt.axis('equal')

toy_data = toy_data.dot(np.diag(grlvq.lambda_))

plt.subplot(424)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(grlvq.w_.dot(np.diag(grlvq.lambda_))[:, 0], grlvq.w_.dot(np.diag(grlvq.lambda_))[:, 1])
plt.axis('equal')

print('grlvq:',grlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))

toy_data = np.append(np.random.multivariate_normal([0, 0], np.array([[5,4],[4,6]]), size=nb_ppc),
                     np.random.multivariate_normal([9, 0], np.array([[5,4],[4,6]]), size=nb_ppc), axis=0)
gmlvq = GmlvqModel()
gmlvq.fit(toy_data,toy_label)
pred = gmlvq.predict(toy_data)

plt.subplot(425)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(gmlvq.w_[:, 0], gmlvq.w_[:, 1])
plt.scatter(gmlvq.w_[:, 0], gmlvq.w_[:, 1], c=gmlvq.c_w_, marker='.')
plt.axis('equal')

toy_data = toy_data.dot(gmlvq.omega_)

plt.subplot(426)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(gmlvq.w_.dot(gmlvq.omega_)[:, 0], gmlvq.w_.dot(gmlvq.omega_)[:, 1])
plt.scatter(gmlvq.w_[:, 0], gmlvq.w_[:, 1], c=gmlvq.c_w_, marker='.')
plt.axis('equal')

print('gmlvq:',gmlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))

toy_data = np.append(np.random.multivariate_normal([0, 1], np.array([[5,-4],[-4,6]]), size=nb_ppc),
                     np.random.multivariate_normal([0, 0], np.array([[5,4],[4,6]]), size=nb_ppc), axis=0)
lgmlvq = LgmlvqModel(display=True)
lgmlvq.fit(toy_data,toy_label)
pred = lgmlvq.predict(toy_data)

plt.subplot(427)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(lgmlvq.w_[:, 0], lgmlvq.w_[:, 1])
plt.axis('equal')

classes = unique_labels(toy_label)
class_points = []
for i in range(len(classes)):
    class_points.append(toy_data[classes[i]==toy_label].dot(lgmlvq.psis_[i]))

plt.subplot(428)
for i in range(len(class_points)):
    plt.scatter(class_points[i][:, 0], class_points[i][:, 1])
    #plt.scatter(class_points[i][:, 0], class_points[i][:, 1], c=pred, marker='.')
    #plt.scatter(glvq.w_[:, 0], glvq.w_[:, 1])
plt.axis('equal')

print('lgmlvq:',lgmlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))
plt.show()