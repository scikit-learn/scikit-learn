
import numpy
from sklearn.earth import Earth
from matplotlib import pyplot

#Create some fake data
numpy.random.seed(2)
m = 10000
n = 10
X = 80*numpy.random.uniform(size=(m,n)) - 40
y = 100*numpy.abs(numpy.sin((X[:,6])/10) - 4.0) + 10*numpy.random.normal(size=m)

#Fit an Earth model
model = Earth(max_degree = 3,minspan_alpha=.5)
model.fit(X,y)

#Print the model
print model.trace()
print model.summary()

#Plot the model
y_hat = model.predict(X)
pyplot.figure()
pyplot.plot(X[:,6],y,'r.')
pyplot.plot(X[:,6],y_hat,'b.')
pyplot.show()

