from sklearn.datasets import load_digits
from sklearn.gp_search import GPSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.linear_model import SGDClassifier
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.pipeline import Pipeline

import logging
import matplotlib.pyplot as plt
import numpy as np

def extend_result(n_tests,tmp_res):
	res = np.zeros(n_tests)
	l = len(tmp_res) -1
	for i in range(n_tests):
		res[i] = tmp_res[min(i,l)]

	return res


def test1():
	iris = load_digits()
	X, y = iris.data, iris.target
	clf = RandomForestClassifier(n_estimators=20)

	# specify parameters and distributions to sample from
	parameters = {"max_depth": [3, 3],
					"max_features": [1,11],
					"min_samples_split": [1,11],
					"min_samples_leaf": [1,11],
					"bootstrap": [True, False],
					"criterion": ["gini", "entropy"]}

	parameters_details = {"max_depth": 'int',
					"max_features": 'int',
					"min_samples_split": 'int',
					"min_samples_leaf": 'int',
					"bootstrap": 'cat',
					"criterion": 'cat'}

	search = GPSearchCV(parameters,parameters_details,estimator=clf,X=X,y=y,n_iter=20)
	search._fit()


def test2():
	parameters = {'kernel' : ['rbf','poly'],'d' : [1,3],'C' : [1,10] }
	parameters_details = {'kernel' : 'cat','d' : 'int','C' : 'float'}
	def scoring_function(x):
		return 0.5

	search = GPSearchCV(parameters,parameters_details,estimator=scoring_function,n_iter=20)
	search._fit()


def test3():
	# Display progress logs on stdout
	logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

	# Load some categories from the training set
	categories = [
	    'alt.atheism',
	    'talk.religion.misc',
	]
	# Uncomment the following to do the analysis on all the categories
	#categories = None

	print("Loading 20 newsgroups dataset for categories:")
	print(categories)

	data = fetch_20newsgroups(subset='train', categories=categories)
	print("%d documents" % len(data.filenames))
	print("%d categories" % len(data.target_names))
	print()

	# define a pipeline combining a text feature extractor with a simple
	# classifier
	pipeline = Pipeline([
	    ('vect', CountVectorizer()),
	    ('tfidf', TfidfTransformer()),
	    ('clf', SGDClassifier()),
	])

	# uncommenting more parameters will give better exploring power but will
	# increase processing time in a combinatorial way
	parameters = {
	    'vect__max_df': [0.5,1.],
	    #'vect__max_features': (None, 5000, 10000, 50000),
	    'vect__ngram_range': [(1, 1), (1, 2)],  # unigrams or bigrams
	    #'tfidf__use_idf': (True, False),
	    #'tfidf__norm': ('l1', 'l2'),
	    'clf__alpha': [0.000001, 0.00001],
	    'clf__penalty': ['l2', 'elasticnet'],
	    #'clf__n_iter': (10, 50, 80),
	}

	parameters_details = {
	    'vect__max_df': 'float',
	    #'vect__max_features': (None, 5000, 10000, 50000),
	    'vect__ngram_range': 'cat' ,  # unigrams or bigrams
	    #'tfidf__use_idf': (True, False),
	    #'tfidf__norm': ('l1', 'l2'),
	    'clf__alpha': 'float',
	    'clf__penalty': 'cat',
	    #'clf__n_iter': (10, 50, 80),
	}

	search = GPSearchCV(parameters,parameters_details,estimator=pipeline,X=data.data, y=data.target,n_iter=20)
	search._fit()


def gp_vs_random_search(test_name,n_tests,search_lenght):
	"""
	Compare GP-based search vs a simple random one
	Choose test_name in {'iris','text'}
	"""

	n_iter_search = search_lenght

	if(test_name == 'iris'):
		iris = load_digits()
		X, y = iris.data, iris.target
		pipeline = RandomForestClassifier(n_estimators=20)

		# specify parameters and distributions to sample from
		parameters = {"max_depth": [2, 15],
						"max_features": [1,20],
						"min_samples_split": [1,20],
						"min_samples_leaf": [1,20],
						"bootstrap": [True, False],
						"criterion": ["gini", "entropy"]}

		parameters_details = {"max_depth": 'int',
						"max_features": 'int',
						"min_samples_split": 'int',
						"min_samples_leaf": 'int',
						"bootstrap": 'cat',
						"criterion": 'cat'}

	elif(test_name == 'text'):
		# Display progress logs on stdout
		logging.basicConfig(level=logging.INFO,
	                    format='%(asctime)s %(levelname)s %(message)s')

		# Load some categories from the training set
		categories = [
		    'alt.atheism',
		    'talk.religion.misc',
		]
		# Uncomment the following to do the analysis on all the categories
		#categories = None
		print("Loading 20 newsgroups dataset for categories:")
		print(categories)

		data = fetch_20newsgroups(subset='train', categories=categories)
		print("%d documents" % len(data.filenames))
		print("%d categories" % len(data.target_names))

		X = data.data
		y = data.target

		# define a pipeline combining a text feature extractor with a simple
		# classifier
		pipeline = Pipeline([
		    ('vect', CountVectorizer()),
		    ('tfidf', TfidfTransformer()),
		    ('clf', SGDClassifier()),
		])

		# uncommenting more parameters will give better exploring power but will
		# increase processing time in a combinatorial way
		parameters = {
		    'vect__max_df': [0.5,1.],
		    #'vect__max_features': (None, 5000, 10000, 50000),
		    'vect__ngram_range': [(1, 1), (1, 2)],  # unigrams or bigrams
		    #'tfidf__use_idf': (True, False),
		    #'tfidf__norm': ('l1', 'l2'),
		    'clf__alpha': [0.000001, 0.00001],
		    'clf__penalty': ['l2', 'elasticnet'],
		    #'clf__n_iter': (10, 50, 80),
		}

		parameters_details = {
		    'vect__max_df': 'float',
		    #'vect__max_features': (None, 5000, 10000, 50000),
		    'vect__ngram_range': 'cat' ,  # unigrams or bigrams
		    #'tfidf__use_idf': (True, False),
		    #'tfidf__norm': ('l1', 'l2'),
		    'clf__alpha': 'float',
		    'clf__penalty': 'cat',
		    #'clf__n_iter': (10, 50, 80),
		}

	else:
		print('Dataset not available for test')

	# GP search
	all_gp_results = []
	print 'GP search'
	for i in range(n_tests):
		search = GPSearchCV(parameters,parameters_details,estimator=pipeline,X=X,y=y,
							n_iter=n_iter_search, n_init=20, verbose=False)
		_,scores = search._fit()

		max_scores = [scores[0]]
		print 'Test',i,'-',len(scores),'parameters tested'

		for j in range(1,len(scores)):
			max_scores.append(max(max_scores[j-1],scores[j]))
		all_gp_results.append(extend_result(n_iter_search,max_scores))
	all_gp_results = np.asarray(all_gp_results)
	print all_gp_results.shape

	# Randomized search
	print 'Random search'
	all_random_results = []
	for i in range(n_tests):
		random_search = GPSearchCV(parameters,parameters_details,estimator=pipeline,X=X,y=y,
                                  	n_iter=n_iter_search, n_init=n_iter_search, verbose=False)
		_,scores = search._fit()

		max_scores = [scores[0]]
		print 'Test',i,'-',len(scores),'parameters tested'

		for j in range(1,len(scores)):
			max_scores.append(max(max_scores[j-1],scores[j]))
		all_random_results.append(extend_result(n_iter_search,max_scores))
	all_random_results = np.asarray(all_random_results)

	plt.figure()
	plt.plot(range(n_iter_search),np.mean(all_gp_results,axis=0),'r',label='GP')
	plt.plot(range(n_iter_search),np.mean(all_random_results,axis=0),'g',label='Random')
	plt.legend()
	plt.title('Test GP vs Random on ' + test_name +' dataset - Average on ' + str(n_tests) + ' trials')
	plt.xlabel('Iterations')
	plt.ylabel('Max CV performance')
	plt.show()



if __name__ == "__main__":
	
	# print 'Routine Test'
	# test2()
	test_name = 'iris'
	n_tests = 20
	search_lenght = 50
	print '\nTest GP vs Random on',test_name,'dataset - Average on',n_tests,'trials'
	gp_vs_random_search(test_name,n_tests,search_lenght)