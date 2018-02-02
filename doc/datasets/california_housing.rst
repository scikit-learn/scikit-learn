
.. _california_housing:

The California Housing dataset
================================


This data set contains a set of housing information taken in the 1990 Census (in the US). The
:func:`sklearn.datasets.fetch_california_housing` function is the data
fetching / caching function that downloads the data
archive from AT&T.

.. _This data set from the original source: http://www.dcc.fc.up.pt/~ltorgo/Regression/cal_housing.html

As described on the original website:

    We collected information on the variables using all the block groups in California from the 1990 Census. In this sample a block group on average includes 1425.5 individuals living in a geographically compact area. Naturally, the geographical area included varies inversely with the population density. We computed distances among the centroids of each block group as measured in latitude and longitude. We excluded all the block groups reporting zero entries for the independent and dependent variables. The final data contained 20,640 observations on 9 variables. The dependent variable is ln(median house value). The file contains all the the variables. Specifically, it contains median house value, median income, housing median age, total rooms, total bedrooms, population, households, latitude, and longitude in that order.

The data set can be accessed using Python by importing like this:


>>> from sklearn.datasets import fetch_california_housing 
>>> cal_housing_data = fetch_california_housing()

Then the user can do what they please with the data. Note that the data
has already been read into the compiler, and using a file importer after this step is not necessary.