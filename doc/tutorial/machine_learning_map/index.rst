.. _ml_map:


.. include:: ../../includes/big_toc_css.rst

Choosing the right estimator
=======================================================


Often the hardest part of solving a machine learning problem can
be finding the right estimator for the job.

Different estimators are better suited for different types of data
and different problems.

The flowchart below is designed to give users a bit of
a rough guide on how to approach problems with regard to
which estimators to try on your data.

Click on any estimator in the chart below to see its documentation.



.. raw:: html

        <img src="../../_static/ml_map.png" class="map" alt="Move mouse over image" usemap="#imgmap">
      	    <map name="imgmap">
	    	<area href="../../documentation.html" title="Back to Documentation" shape="poly" coords="97,1094, 76,1097, 56,1105, 40,1120, 35,1132, 34,1145, 35,1153, 40,1162, 46,1171, 54,1177, 62,1182, 72,1187, 81,1188, 100,1189, 118,1186, 127,1182, 136,1177, 146,1170, 152,1162, 155,1158, 158,1146, 158,1126, 143,1110, 138,1105, 127,1100, 97,1094"></area>
		<area href="../../modules/linear_model.html#elastic-net" title="Elastic Net Documentation" shape="poly" coords="1556,446, 1556,446, 1556,476, 1556,476, 1556,476, 1676,476, 1676,476, 1676,476, 1676,446, 1676,446, 1676,446, 1556,446, 1556,446" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/ensemble.html" title="Ensembe Methods Documentation" shape="poly" coords="209,200, 209,200, 209,252, 209,252, 209,252, 332,252, 332,252, 332,252, 332,200, 332,200, 332,200, 209,200, 209,200" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/ensemble.html" title="Ensembe Methods Documentation" shape="poly" coords="1828,506, 1828,506, 1828,544, 1828,544, 1828,544, 2054,544, 2054,544, 2054,544, 2054,506, 2054,506, 2054,506, 1828,506, 1828,506" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/mixture.html" title="Gaussian mixture models Documentation" shape="poly" coords="142,637, 142,637, 142,667, 142,667, 142,667, 265,667, 265,667, 265,667, 265,637, 265,637, 265,637, 142,637, 142,637" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/manifold.html#isomap" title="Isomap Documentation" shape="poly" coords="1500,799, 1500,799, 1500,844, 1500,844, 1500,844, 1618,844, 1618,844, 1618,844, 1618,800, 1618,800, 1618,800, 1500,799, 1500,799" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/kernel_approximation.html" title="Kernel Approximation Documentation" shape="poly" coords="1477,982, 1477,982, 1477,1055, 1477,1055, 1477,1055, 1638,1055, 1638,1055, 1638,1055, 1638,982, 1638,982, 1638,982, 1477,982, 1477,982" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/kernel_approximation.html" title="Kernel Approximation Documentation" shape="poly" coords="472,100, 472,100, 472,173, 472,173, 472,173, 634,173, 634,173, 634,173, 634,100, 634,100, 634,100, 472,100, 472,100" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/clustering.html#k-means" title="KMeans Documentation" shape="poly" coords="377,605, 377,605, 377,655, 377,655, 377,655, 476,655, 476,655, 476,655, 476,605, 476,605, 476,605, 377,605, 377,605" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/neighbors.html" title="Nearest Neighbors" shape="poly" coords="440,219, 440,219, 440,293, 440,293, 440,293, 574,293, 574,293, 574,293, 574,219, 574,219, 574,219, 440,219, 440,219" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/linear_model.html#lasso" title="Lasso Documentation" shape="poly" coords="1550,408, 1550,408, 1550,436, 1550,436, 1550,436, 1671,436, 1671,436, 1671,436, 1671,408, 1671,408, 1671,408, 1550,408, 1550,408" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/svm.html#classification" title="LinearSVC Documentation" shape="poly" coords="609,419, 609,419, 609,492, 609,492, 609,492, 693,492, 693,492, 693,492, 693,419, 693,419, 693,419, 609,419, 609,419" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/manifold.html#locally-linear-embedding" title="Locally Linear Embedding Documentation" shape="poly" coords="1719,888, 1719,888, 1719,945, 1719,945, 1719,945, 1819,945, 1819,945, 1819,945, 1819,888, 1819,888, 1819,888, 1719,888, 1719,888" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/clustering.html#mean-shift" title="Mean Shift Documentation" shape="poly" coords="562,949, 562,949, 562,981, 562,981, 562,981, 682,981, 682,981, 682,981, 682,949, 682,949, 682,949, 562,949, 562,949" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/clustering.html#mini-batch-k-means" title="Mini Batch K-means Documentation" shape="poly" coords="343,917, 343,917, 343,990, 343,990, 343,990, 461,990, 461,990, 461,990, 461,917, 461,917, 461,917, 343,917, 343,917" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/naive_bayes.html" title="Naive Bayes Documentation" shape="poly" coords="194,339, 194,339, 194,412, 194,412, 194,412, 294,412, 294,412, 294,412, 294,339, 294,339, 294,339, 194,339, 194,339" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/decomposition.html#principal-component-analysis-pca" title="Principal Component Analysis Documentation" shape="poly" coords="1208,778, 1208,778, 1208,851, 1208,851, 1208,851, 1350,851, 1350,851, 1350,851, 1350,778, 1350,778, 1350,778, 1208,778, 1208,778" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/linear_model.html#ridge-regression" title="Ridge Regression Documentation" shape="poly" coords="1696,648, 1696,648, 1696,687, 1696,687, 1696,687, 1890,687, 1890,687, 1890,687, 1890,648, 1890,648, 1890,648, 1696,648, 1696,648" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/sgd.html#classification" title="SGD Classifier Documentation" shape="poly" coords="691,205, 691,205, 691,278, 691,278, 691,278, 803,278, 803,278, 803,278, 803,205, 803,205, 803,205, 691,205, 691,205" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/sgd.html#regression" title="SGD Regression Documentation" shape="poly" coords="1317,425, 1317,425, 1317,498, 1317,498, 1317,498, 1436,498, 1436,498, 1436,498, 1436,425, 1436,425, 1436,425, 1317,425, 1317,425" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/clustering.html#spectral-clustering" title="Spectral Clustering Documentation" shape="poly" coords="145,572, 145,572, 145,631, 145,631, 145,631, 267,631, 267,631, 267,631, 267,572, 267,572, 267,572, 145,572, 145,572" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/manifold.html#spectral-embedding" title="Spectral Embedding Documentation" shape="poly" coords="1502,849, 1502,849, 1502,910, 1502,910, 1502,910, 1618,910, 1618,910, 1618,910, 1618,849, 1618,849, 1618,849, 1502,849, 1502,849" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/svm.html#classification" title="SVC Documentation" shape="poly" coords="210,157, 210,157, 210,194, 210,194, 210,194, 333,194, 333,194, 333,194, 333,157, 333,157, 333,157, 210,157, 210,157" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/svm.html#regression" title="SVR Documentation" shape="poly" coords="1696,692, 1696,692, 1696,732, 1696,732, 1696,732, 1890,732, 1890,732, 1890,732, 1890,692, 1890,692, 1890,692, 1696,692, 1696,692" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/svm.html#regression" title="SVR Documentation" shape="poly" coords="1831,458, 1831,458, 1831,496, 1831,496, 1831,496, 2052,496, 2052,496, 2052,496, 2052,458, 2052,458, 2052,458, 1831,458, 1831,458" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
		<area href="../../modules/mixture.html#bgmm" title=" Bayesian GMM Documentation" shape="poly" coords="562,994, 562,994, 562,1026, 562,1026, 562,1026, 682,1026, 682,1026, 682,1026, 682,994, 682,994, 682,994, 562,994, 562,994" data-maphilight='{"strokeColor":"0000ff","strokeWidth":5,"fillColor":"66FF66","fillOpacity":0.4}'></area>
	    </map>
	</img>
