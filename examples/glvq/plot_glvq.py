import numpy as np
import matplotlib.pyplot as plt

from sklearn.glvq.plotting import project_plot2d
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

f = plt.figure(1)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=toy_label)
plt.scatter(toy_data[:, 0], toy_data[:, 1], c=pred, marker='.')
plt.scatter(glvq.w_[:, 0], glvq.w_[:, 1])
plt.axis('equal')
f.show()

toy_data = np.append(np.random.multivariate_normal([0, 0], np.array([[0.3,0],[0,4]]), size=nb_ppc),
                     np.random.multivariate_normal([4, 4], np.array([[0.3,0],[0,4]]), size=nb_ppc), axis=0)
grlvq = GrlvqModel()
grlvq.fit(toy_data,toy_label)
project_plot2d(grlvq,toy_data,toy_label,2)

print('grlvq:',grlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))

#toy_data = np.append(np.random.multivariate_normal([0, 0], np.array([[5,4],[4,6]]), size=nb_ppc),
#                     np.random.multivariate_normal([9, 0], np.array([[5,4],[4,6]]), size=nb_ppc), axis=0)
gmlvq = GmlvqModel()
gmlvq.fit(toy_data,toy_label)
project_plot2d(gmlvq,toy_data,toy_label,3)

print('gmlvq:',gmlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))

toy_data = np.append(np.random.multivariate_normal([0, 1], np.array([[5,-4],[-4,6]]), size=nb_ppc),
                     np.random.multivariate_normal([0, 0], np.array([[5,4],[4,6]]), size=nb_ppc), axis=0)
lgmlvq = LgmlvqModel()
lgmlvq.fit(toy_data,toy_label)
project_plot2d(lgmlvq, toy_data, toy_label, 4)

print('lgmlvq:',lgmlvq.score(toy_data,toy_label))
print('gvlq:', glvq.score(toy_data,toy_label))
plt.show()
