# The _DEFAULT_PARAM_GRIDS was shamelessly stolen from
# https://github.com/EducationalTestingService/skll/blob/master/skll/learner.py

_DEFAULT_PARAM_GRIDS = {'AdaBoostClassifier':
                        [{'learning_rate': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'AdaBoostRegressor':
                        [{'learning_rate': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'DecisionTreeClassifier':
                        [{'max_features': ["auto", None]}],
                        'DecisionTreeRegressor':
                        [{'max_features': ["auto", None]}],
                        'ElasticNet':
                        [{'alpha': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'GradientBoostingClassifier':
                        [{'max_depth': [1, 3, 5]}],
                        'GradientBoostingRegressor':
                        [{'max_depth': [1, 3, 5]}],
                        'KNeighborsClassifier':
                        [{'n_neighbors': [1, 5, 10, 100],
                          'weights': ['uniform', 'distance']}],
                        'KNeighborsRegressor':
                        [{'n_neighbors': [1, 5, 10, 100],
                          'weights': ['uniform', 'distance']}],
                        'Lasso':
                        [{'alpha': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'LinearRegression':
                        [{}],
                        'LinearSVC':
                        [{'C': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'LogisticRegression':
                        [{'C': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'SVC': [{'C': [0.01, 0.1, 1.0, 10.0, 100.0],
                                 'gamma': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'MultinomialNB':
                        [{'alpha': [0.1, 0.25, 0.5, 0.75, 1.0]}],
                        'RandomForestClassifier':
                        [{'max_depth': [1, 5, 10, None]}],
                        'RandomForestRegressor':
                        [{'max_depth': [1, 5, 10, None]}],
                        'Ridge':
                        [{'alpha': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'SGDClassifier':
                        [{'alpha': [0.000001, 0.00001, 0.0001, 0.001, 0.01],
                          'penalty': ['l1', 'l2', 'elasticnet']}],
                        'SGDRegressor':
                        [{'alpha': [0.000001, 0.00001, 0.0001, 0.001, 0.01],
                          'penalty': ['l1', 'l2', 'elasticnet']}],
                        'LinearSVR':
                        [{'C': [0.01, 0.1, 1.0, 10.0, 100.0]}],
                        'SVR':
                        [{'C': [0.01, 0.1, 1.0, 10.0, 100.0],
                          'gamma': [0.01, 0.1, 1.0, 10.0, 100.0]}]}
