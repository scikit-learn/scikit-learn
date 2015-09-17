# -*- coding: utf-8 -*-
"""
Concordance correlation coefficient VS Pearson's correlation coefficient.

The predicted values are equal to x+noise+bias. As the bias increased from
1 to 100 the Pearson correlation coefficient remains unaffected, but the
concordance correlation coefficient tends to zero.

This demonstrates the appropriateness of the concordance correlation coefficient
over Pearson's correlation coefficient in situations where the predicted values
might suffer from high bias.
"""
import matplotlib.pyplot as plt

from metrics import *
n_samples=10

pearson=[]
concor=[]

y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
y_pred = np.array([[0, 0, 0, 1], [1, 0, 1, 1], [0, 0, 0, 1]])
    
print(pearson_correlation_coefficient(y_true,y_pred))


for i in range(1,100):
    y_true = np.arange(n_samples)
    y_pred = y_true + i+np.random.rand(n_samples)*5
    
    pearson.append(pearson_correlation_coefficient(y_true,y_pred))
    concor.append(concordance_correlation_coefficient(y_true,y_pred))
    

    
plt.figure()
plt.plot(range(1,100),pearson,range(1,100),concor)
plt.legend(["pearson","concordance"],loc="best")