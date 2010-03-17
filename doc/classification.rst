Classification
==============


Classifying data is a common task in machine learning. Suppose some
given data points each belong to one of two classes, and the goal is
to decide which class a new data point will be in.

In the case of support vector machines, a data point is viewed as a
p-dimensional vector (a list of p numbers), and we want to know
whether we can separate such points with a p − 1-dimensional
hyperplane. This is called a linear classifier. There are many
hyperplanes that might classify the data. One reasonable choice as the
best hyperplane is the one that represents the largest separation, or
margin, between the two classes. So we choose the hyperplane so that
the distance from it to the nearest data point on each side is
maximized. If such a hyperplane exists, it is known as the
maximum-margin hyperplane and the linear classifier it defines is
known as a maximum margin classifier.

