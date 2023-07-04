ការចាប់ផ្តើម
គោលបំណងនៃការណែនាំនេះគឺដើម្បីបង្ហាញពីលក្ខណៈសំខាន់ៗមួយចំនួនដែល scikit-learn ផ្តល់ជូន។ វាសន្មត់ថាជាចំណេះដឹងការងារជាមូលដ្ឋាននៃការអនុវត្តការរៀនម៉ាស៊ីន (ការបំពេញគំរូ ការទស្សន៍ទាយ សុពលភាពឆ្លងកាត់។ល។)។ សូមយោងទៅលើ :ref:`ការណែនាំអំពីការដំឡើង <installation-instructions>` របស់យើងសម្រាប់ការដំឡើង scikit-learn ។

Scikit-learn គឺជាបណ្ណាល័យការរៀនម៉ាស៊ីនប្រភពបើកចំហដែលគាំទ្រការរៀនសូត្រដែលមានការគ្រប់គ្រង និងគ្មានការត្រួតពិនិត្យ។ វាក៏ផ្តល់នូវឧបករណ៍ជាច្រើនសម្រាប់ការបំពេញគំរូ ដំណើរការទិន្នន័យជាមុន ការជ្រើសរើសគំរូ ការវាយតម្លៃគំរូ និងឧបករណ៍ប្រើប្រាស់ជាច្រើនទៀត។

សមនិងទស្សន៍ទាយ៖ មូលដ្ឋានគ្រឹះប៉ាន់ស្មាន
Scikit-learn ផ្ដល់ជូននូវក្បួនដោះស្រាយ និងគំរូនៃការរៀនម៉ាស៊ីនដែលភ្ជាប់មកជាមួយជាច្រើន ដែលហៅថា :term:`estimators`។ ការប៉ាន់ស្មាននីមួយៗអាចត្រូវបានសមទៅនឹងទិន្នន័យមួយចំនួនដោយប្រើវិធី :term:`fit` របស់វា។

នេះគឺជាឧទាហរណ៍ដ៏សាមញ្ញមួយដែលយើងសមនឹង :class:`~sklearn.ensemble.RandomForestClassifier` ទៅនឹងទិន្នន័យមូលដ្ឋានមួយចំនួន៖
>>> from sklearn.ensemble import RandomForestClassifier
>>> clf = RandomForestClassifier(random_state=0)
>>> X = [[ 1,  2,  3],  # 2 samples, 3 features
...      [11, 12, 13]]
>>> y = [0, 1]  # classes of each sample
>>> clf.fit(X, y)
RandomForestClassifier(random_state=0)

វិធីសាស្រ្ត :term:`fit` ជាទូទៅទទួលយក 2 ធាតុចូល៖

ម៉ាទ្រីសគំរូ (ឬម៉ាទ្រីសរចនា): ពាក្យ៖`X`។ ទំហំ X គឺជាធម្មតា (n_samples, n_features) ដែលមានន័យថាគំរូត្រូវបានតំណាងជាជួរ ហើយលក្ខណៈពិសេសត្រូវបានតំណាងជាជួរឈរ។
តម្លៃគោលដៅ៖ ពាក្យ៖`y` ដែលជាចំនួនពិតសម្រាប់កិច្ចការតំរែតំរង់ ឬចំនួនគត់សម្រាប់ចាត់ថ្នាក់ (ឬសំណុំតម្លៃដាច់ពីគ្នាផ្សេងទៀត)។ សម្រាប់កិច្ចការសិក្សាដែលមិនមានការត្រួតពិនិត្យ y មិនចាំបាច់បញ្ជាក់ទេ។ y ជាធម្មតាជាអារេ 1d ដែលធាតុ i th ត្រូវគ្នាទៅនឹងគោលដៅនៃគំរូ i th (ជួរដេក) នៃ X ។
ទាំង X និង y ជាធម្មតាត្រូវបានគេរំពឹងថាជាអារេ numpy ឬសមមូល :term:`array-like` data type ទោះបីជាការប៉ាន់ស្មានខ្លះដំណើរការជាមួយទម្រង់ផ្សេងទៀតដូចជា sparse matrices។

នៅពេលដែលឧបករណ៍ប៉ាន់ស្មានត្រូវបានបំពាក់ វាអាចត្រូវបានប្រើសម្រាប់ការទស្សន៍ទាយតម្លៃគោលដៅនៃទិន្នន័យថ្មី។ អ្នកមិនចាំបាច់បង្វឹកអ្នកប៉ាន់ស្មានឡើងវិញទេ៖
>>> clf.predict(X)  # predict classes of the training data
array([0, 1])
>>> clf.predict([[4, 5, 6], [14, 15, 16]])  # predict classes of new data
array([0, 1])
Transformers និង Pre-processors
លំហូរការងារនៃការរៀនម៉ាស៊ីនជាញឹកញាប់ត្រូវបានផ្សំឡើងដោយផ្នែកផ្សេងៗ។ បំពង់បង្ហូរប្រេងធម្មតាមានជំហានដំណើរការមុនដែលបំប្លែង ឬបញ្ចូលទិន្នន័យ និងអ្នកព្យាករណ៍ចុងក្រោយដែលព្យាករណ៍តម្លៃគោលដៅ។

នៅក្នុង scikit-learn, pre-processors និង transformers ដើរតាម API ដូចគ្នាទៅនឹង objects ប៉ាន់ស្មាន (តាមពិតពួកវាទាំងអស់ទទួលបានមរតកពី BaseEstimator class ដូចគ្នា)។ វត្ថុបំប្លែងមិនមានវិធីសាស្ត្រ :term:`predict` ទេ ប៉ុន្តែជាវិធីសាស្ត្រ :term:`transform` ដែលបញ្ចេញនូវម៉ាទ្រីសគំរូ X ដែលបានផ្លាស់ប្តូរថ្មី៖
>>> from sklearn.preprocessing import StandardScaler
>>> X = [[0, 15],
...      [1, -10]]
>>> # scale data according to computed scaling values
>>> StandardScaler().fit(X).transform(X)
array([[-1.,  1.],
       [ 1., -1.]])
ពេលខ្លះ អ្នកចង់អនុវត្តការបំប្លែងផ្សេងៗគ្នាចំពោះមុខងារផ្សេងៗគ្នា៖ :ref:`ColumnTransformer<column_transformer>` ត្រូវបានរចនាឡើងសម្រាប់ករណីប្រើប្រាស់ទាំងនេះ។

បំពង់៖ ច្រវាក់អ្នកកែច្នៃមុន និងអ្នកប៉ាន់ប្រមាណ
Transformers និងអ្នកប៉ាន់ប្រមាណ (ទស្សន៍ទាយ) អាចត្រូវបានផ្សំជាមួយគ្នាទៅជាវត្ថុបង្រួបបង្រួមតែមួយ៖ a :class:`~sklearn.pipeline.Pipeline`។ បំពង់ផ្តល់នូវ API ដូចគ្នាទៅនឹងការប៉ាន់ស្មានធម្មតា៖ វាអាចត្រូវបានបំពាក់ និងប្រើសម្រាប់ការទស្សន៍ទាយដោយសម និងទស្សន៍ទាយ។ ដូចដែលយើងនឹងឃើញនៅពេលក្រោយ ការប្រើប្រាស់បំពង់មួយក៏នឹងការពារអ្នកពីការលេចធ្លាយទិន្នន័យ ពោលគឺការបង្ហាញទិន្នន័យសាកល្បងមួយចំនួននៅក្នុងទិន្នន័យបណ្តុះបណ្តាលរបស់អ្នក។

ក្នុងឧទាហរណ៍ខាងក្រោម យើង៖ ref:`ផ្ទុកសំណុំទិន្នន័យ Iris <datasets>` បំបែកវាទៅជារថភ្លើង និងឈុតសាកល្បង ហើយគណនាពិន្ទុភាពត្រឹមត្រូវនៃបំពង់បង្ហូរប្រេងលើទិន្នន័យសាកល្បង៖
>>> from sklearn.preprocessing import StandardScaler
>>> from sklearn.linear_model import LogisticRegression
>>> from sklearn.pipeline import make_pipeline
>>> from sklearn.datasets import load_iris
>>> from sklearn.model_selection import train_test_split
>>> from sklearn.metrics import accuracy_score
...
>>> # create a pipeline object
>>> pipe = make_pipeline(
...     StandardScaler(),
...     LogisticRegression()
... )
...
>>> # load the iris dataset and split it into train and test sets
>>> X, y = load_iris(return_X_y=True)
>>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
...
>>> # fit the whole pipeline
>>> pipe.fit(X_train, y_train)
Pipeline(steps=[('standardscaler', StandardScaler()),
                ('logisticregression', LogisticRegression())])
>>> # we can now use it like any other estimator
>>> accuracy_score(pipe.predict(X_test), y_test)
0.97...
ការវាយតម្លៃគំរូ
ការបញ្ចូលគំរូទៅនឹងទិន្នន័យមួយចំនួនមិនមានន័យថាវានឹងអាចទស្សន៍ទាយបានល្អលើទិន្នន័យដែលមើលមិនឃើញនោះទេ។ នេះចាំបាច់ត្រូវវាយតម្លៃដោយផ្ទាល់។ យើងទើបតែបានឃើញ :func:`~sklearn.model_selection.train_test_split` ជំនួយដែលបំបែកសំណុំទិន្នន័យទៅជាសំណុំរថភ្លើង និងតេស្ត ប៉ុន្តែ scikit-learn ផ្តល់នូវឧបករណ៍ជាច្រើនទៀតសម្រាប់ការវាយតម្លៃគំរូ ជាពិសេសសម្រាប់ :ref:`cross-validation < cross_validation>`។

នៅទីនេះ យើងបង្ហាញយ៉ាងខ្លីពីរបៀបអនុវត្តនីតិវិធីឆ្លងកាត់សុពលភាព 5 ដង ដោយប្រើ :func:`~sklearn.model_selection.cross_validate` ជំនួយ។ ចំណាំថាវាក៏អាចធ្វើម្តងទៀតដោយដៃនៅលើផ្នត់ ប្រើយុទ្ធសាស្រ្តបំបែកទិន្នន័យផ្សេងគ្នា និងប្រើមុខងារដាក់ពិន្ទុផ្ទាល់ខ្លួន។ សូមមើល :ref:`User Guide <cross_validation>` របស់យើងសម្រាប់ព័ត៌មានលម្អិតបន្ថែម៖
>>> from sklearn.datasets import make_regression
>>> from sklearn.linear_model import LinearRegression
>>> from sklearn.model_selection import cross_validate
...
>>> X, y = make_regression(n_samples=1000, random_state=0)
>>> lr = LinearRegression()
...
>>> result = cross_validate(lr, X, y)  # defaults to 5-fold CV
>>> result['test_score']  # r_squared score is high because dataset is easy
array([1., 1., 1., 1., 1.])
ស្វែងរកប៉ារ៉ាម៉ែត្រដោយស្វ័យប្រវត្តិ
រាល់ការប៉ាន់ប្រមាណមានប៉ារ៉ាម៉ែត្រ (ជារឿយៗគេហៅថា hyper-parameters នៅក្នុងអក្សរសិល្ប៍) ដែលអាចកែតម្រូវបាន។ អំណាចទូទៅនៃការប៉ាន់ប្រមាណ ច្រើនតែអាស្រ័យទៅលើប៉ារ៉ាម៉ែត្រមួយចំនួន។ ឧទាហរណ៍ :class:`~sklearn.ensemble.RandomForestRegressor` មានប៉ារ៉ាម៉ែត្រ n_estimators ដែលកំណត់ចំនួនដើមឈើនៅក្នុងព្រៃ និងប៉ារ៉ាម៉ែត្រ max_depth ដែលកំណត់ជម្រៅអតិបរមានៃដើមឈើនីមួយៗ។ ជាញឹកញយ វាមិនច្បាស់ថាតម្លៃពិតប្រាកដនៃប៉ារ៉ាម៉ែត្រទាំងនេះគួរតែជាអ្វីទេ ដោយសារវាអាស្រ័យលើទិន្នន័យនៅនឹងដៃ។

Scikit-learn ផ្តល់ឧបករណ៍ដើម្បីស្វែងរកបន្សំប៉ារ៉ាម៉ែត្រល្អបំផុតដោយស្វ័យប្រវត្តិ (តាមរយៈសុពលភាពឆ្លង)។ ក្នុងឧទាហរណ៍ខាងក្រោម យើងស្វែងរកដោយចៃដន្យលើចន្លោះប៉ារ៉ាម៉ែត្រនៃព្រៃចៃដន្យជាមួយនឹងវត្ថុ :class:`~sklearn.model_selection.RandomizedSearchCV`។ នៅពេលការស្វែងរកបានបញ្ចប់ :class:`~sklearn.model_selection.RandomizedSearchCV` ដើរតួជា :class:`~sklearn.ensemble.RandomForestRegressor` ដែលត្រូវបានបំពាក់ជាមួយនឹងសំណុំប៉ារ៉ាម៉ែត្រល្អបំផុត។ អានបន្ថែមនៅក្នុង :ref:`User Guide <grid_search>`:

>>> from sklearn.datasets import fetch_california_housing
>>> from sklearn.ensemble import RandomForestRegressor
>>> from sklearn.model_selection import RandomizedSearchCV
>>> from sklearn.model_selection import train_test_split
>>> from scipy.stats import randint
...
>>> X, y = fetch_california_housing(return_X_y=True)
>>> X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
...
>>> # define the parameter space that will be searched over
>>> param_distributions = {'n_estimators': randint(1, 5),
...                        'max_depth': randint(5, 10)}
...
>>> # now create a searchCV object and fit it to the data
>>> search = RandomizedSearchCV(estimator=RandomForestRegressor(random_state=0),
...                             n_iter=5,
...                             param_distributions=param_distributions,
...                             random_state=0)
>>> search.fit(X_train, y_train)
RandomizedSearchCV(estimator=RandomForestRegressor(random_state=0), n_iter=5,
                   param_distributions={'max_depth': ...,
                                        'n_estimators': ...},
                   random_state=0)
>>> search.best_params_
{'max_depth': 9, 'n_estimators': 4}

>>> # the search object now acts like a normal random forest estimator
>>> # with max_depth=9 and n_estimators=4
>>> search.score(X_test, y_test)
0.73...

ចំណាំ

នៅក្នុងការអនុវត្ត អ្នកស្ទើរតែតែងតែចង់ :ref:`ស្វែងរកតាមបំពង់ <composite_grid_search>` ជំនួសឱ្យការប៉ាន់ស្មានតែមួយ។ មូលហេតុចម្បងមួយគឺថា ប្រសិនបើអ្នកអនុវត្តជំហានមុនដំណើរការទៅសំណុំទិន្នន័យទាំងមូលដោយមិនប្រើបំពង់បង្ហូរប្រេង ហើយបន្ទាប់មកអនុវត្តប្រភេទនៃសុពលភាពឆ្លងណាមួយ នោះអ្នកនឹងបំបែកការសន្មតជាមូលដ្ឋាននៃភាពឯករាជ្យរវាងទិន្នន័យបណ្តុះបណ្តាល និងការធ្វើតេស្ត។ ជាការពិតណាស់ ចាប់តាំងពីអ្នកបានដំណើរការទិន្នន័យជាមុនដោយប្រើសំណុំទិន្នន័យទាំងមូល ព័ត៌មានមួយចំនួនអំពីសំណុំតេស្តមានសម្រាប់ឈុតរថភ្លើង។ នេះនឹងនាំទៅដល់ការប៉ាន់ប្រមាណលើសអំណាចទូទៅរបស់អ្នកប៉ាន់ស្មាន (អ្នកអាចអានបន្ថែមនៅក្នុងប្រកាស Kaggle នេះ)។

ការប្រើប្រាស់បំពង់បង្ហូរប្រេងសម្រាប់ការផ្ទៀងផ្ទាត់ឆ្លងដែន និងការស្វែងរកយ៉ាងទូលំទូលាយនឹងរារាំងអ្នកពីបញ្ហាទូទៅនេះ។

ជំហានបន្ទាប់
យើងបានគ្របដណ្តប់ដោយសង្ខេបអំពីការសមនឹងការប៉ាន់ស្មាន និងការទស្សន៍ទាយ ជំហានមុនដំណើរការ បំពង់បង្ហូរ ឧបករណ៍ដែលមានសុពលភាព និងការស្វែងរកប៉ារ៉ាម៉ែត្រខ្ពស់ដោយស្វ័យប្រវត្តិ។ មគ្គុទ្ទេសក៍នេះគួរតែផ្តល់ឱ្យអ្នកនូវទិដ្ឋភាពទូទៅនៃលក្ខណៈសំខាន់ៗមួយចំនួននៃបណ្ណាល័យ ប៉ុន្តែនៅមានច្រើនទៀតសម្រាប់ scikit-learn!

សូមមើល :ref:`user_guide` របស់យើងសម្រាប់ព័ត៌មានលម្អិតអំពីឧបករណ៍ទាំងអស់ដែលយើងផ្តល់។ អ្នកក៏អាចស្វែងរកបញ្ជីពេញលេញនៃ API សាធារណៈនៅក្នុង :ref:`api_ref`។

អ្នកក៏អាចមើល :ref:`examples <general_examples>` ជាច្រើនរបស់យើង ដែលបង្ហាញពីការប្រើប្រាស់ scikit-learn ក្នុងបរិបទផ្សេងៗគ្នាជាច្រើន។

ឯកសារ :ref:`tutorials <tutorial_menu>` ក៏មានធនធានសិក្សាបន្ថែមផងដែរ។
