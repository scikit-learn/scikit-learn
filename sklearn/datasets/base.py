"""
Base IO code for all datasets
"""

# Copyright (c) 2007 David Cournapeau <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
#               2010 Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause
import sys
import shutil
from collections import namedtuple
from os import environ, listdir, makedirs
from os.path import dirname, exists, expanduser, isdir, join, splitext
import hashlib

from ..utils import Bunch
from ..utils import check_random_state

import numpy as np
from numpy import ndarray

from urllib.request import urlretrieve

from typing import Dict, Tuple, Union, Any
from abc import ABCMeta, abstractmethod
from ..utils import deprecated

RemoteFileMetadata = namedtuple('RemoteFileMetadata',
                                ['filename', 'url', 'checksum'])
Dataset = Union[Bunch, Tuple[Any, Any]]


def get_data_home(data_home=None):
    """Return the path of the scikit-learn data dir.

    This folder is used by some large dataset loaders to avoid downloading the
    data several times.

    By default the data dir is set to a folder named 'scikit_learn_data' in the
    user home folder.

    Alternatively, it can be set by the 'SCIKIT_LEARN_DATA' environment
    variable or programmatically by giving an explicit folder path. The '~'
    symbol is expanded to the user home folder.

    If the folder does not already exist, it is automatically created.

    Parameters
    ----------
    data_home : str | None
        The path to scikit-learn data dir.
    """
    if data_home is None:
        data_home = environ.get('SCIKIT_LEARN_DATA',
                                join('~', 'scikit_learn_data'))
    data_home = expanduser(data_home)
    if not exists(data_home):
        makedirs(data_home)
    return data_home


def clear_data_home(data_home=None):
    """Delete all the content of the data home cache.

    Parameters
    ----------
    data_home : str | None
        The path to scikit-learn data dir.
    """
    data_home = get_data_home(data_home)
    shutil.rmtree(data_home)


def load_files(container_path, description=None, categories=None,
               load_content=True, shuffle=True, encoding=None,
               decode_error='strict', random_state=0):
    """Load text files with categories as subfolder names.

    Individual samples are assumed to be files stored a two levels folder
    structure such as the following:

        container_folder/
            category_1_folder/
                file_1.txt
                file_2.txt
                ...
                file_42.txt
            category_2_folder/
                file_43.txt
                file_44.txt
                ...

    The folder names are used as supervised signal label names. The individual
    file names are not important.

    This function does not try to extract features into a numpy array or scipy
    sparse matrix. In addition, if load_content is false it does not try to
    load the files in memory.

    To use text files in a scikit-learn classification or clustering algorithm,
    you will need to use the `sklearn.feature_extraction.text` module to build
    a feature extraction transformer that suits your problem.

    If you set load_content=True, you should also specify the encoding of the
    text using the 'encoding' parameter. For many modern text files, 'utf-8'
    will be the correct encoding. If you leave encoding equal to None, then the
    content will be made of bytes instead of Unicode, and you will not be able
    to use most functions in `sklearn.feature_extraction.text`.

    Similar feature extractors should be built for other kind of unstructured
    data input such as images, audio, video, ...

    Read more in the :ref:`User Guide <datasets>`.

    Parameters
    ----------
    container_path : string or unicode
        Path to the main folder holding one subfolder per category

    description : string or unicode, optional (default=None)
        A paragraph describing the characteristic of the dataset: its source,
        reference, etc.

    categories : A collection of strings or None, optional (default=None)
        If None (default), load all the categories. If not None, list of
        category names to load (other categories ignored).

    load_content : boolean, optional (default=True)
        Whether to load or not the content of the different files. If true a
        'data' attribute containing the text information is present in the data
        structure returned. If not, a filenames attribute gives the path to the
        files.

    shuffle : bool, optional (default=True)
        Whether or not to shuffle the data: might be important for models that
        make the assumption that the samples are independent and identically
        distributed (i.i.d.), such as stochastic gradient descent.

    encoding : string or None (default is None)
        If None, do not try to decode the content of the files (e.g. for images
        or other non-text content). If not None, encoding to use to decode text
        files to Unicode if load_content is True.

    decode_error : {'strict', 'ignore', 'replace'}, optional
        Instruction on what to do if a byte sequence is given to analyze that
        contains characters not of the given `encoding`. Passed as keyword
        argument 'errors' to bytes.decode.

    random_state : int, RandomState instance or None (default=0)
        Determines random number generation for dataset shuffling. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are: either
        data, the raw text data to learn, or 'filenames', the files
        holding it, 'target', the classification labels (integer index),
        'target_names', the meaning of the labels, and 'DESCR', the full
        description of the dataset.
    """
    target = []
    target_names = []
    filenames = []

    folders = [f for f in sorted(listdir(container_path))
               if isdir(join(container_path, f))]

    if categories is not None:
        folders = [f for f in folders if f in categories]

    for label, folder in enumerate(folders):
        target_names.append(folder)
        folder_path = join(container_path, folder)
        documents = [join(folder_path, d)
                     for d in sorted(listdir(folder_path))]
        target.extend(len(documents) * [label])
        filenames.extend(documents)

    # convert to array for fancy indexing
    filenames = np.array(filenames)
    target = np.array(target)

    if shuffle:
        random_state = check_random_state(random_state)
        indices = np.arange(filenames.shape[0])
        random_state.shuffle(indices)
        filenames = filenames[indices]
        target = target[indices]

    if load_content:
        data = []
        for filename in filenames:
            with open(filename, 'rb') as f:
                data.append(f.read())
        if encoding is not None:
            data = [d.decode(encoding, decode_error) for d in data]
        return Bunch(data=data,
                     filenames=filenames,
                     target_names=target_names,
                     target=target,
                     DESCR=description)

    return Bunch(filenames=filenames,
                 target_names=target_names,
                 target=target,
                 DESCR=description)


@deprecated("'load_data' was renamed to"
            "'SimpleCSVLocalDatasetLoader.read_X_y_csv'"
            "in version 0.21 and will be removed in 0.23.")
def load_data(path: str) -> Tuple[ndarray, ndarray, ndarray]:
    """Loads data from module_path/data/data_file_name.

    Parameters
    ----------
    path : string
        The data file path.

    Returns
    -------
    data : Numpy array
        A 2D array with each row representing one sample and each column
        representing the features of a given sample.

    target : Numpy array
        A 1D array holding target variables for all the samples in `data.
        For example target[0] is the target varible for data[0].

    target_names : Numpy array
        A 1D array containing the names of the classifications. For example
        target_names[0] is the name of the target[0] class.
    """
    return SimpleCSVLocalDatasetLoader.read_X_y_csv(path=path)


def _attempt_cast_to_int(arr: ndarray) -> ndarray:
    arri = arr.astype('int', casting='unsafe')
    if (arr == arri).all():
        return arri
    return arr


@deprecated("'load_wine' was renamed to 'Wine().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_wine(return_X_y=False):
    """Load and return the wine dataset (classification).

    .. versionadded:: 0.18

    The wine dataset is a classic and very easy multi-class classification
    dataset.

    =================   ==============
    Classes                          3
    Samples per class        [59,71,48]
    Samples total                  178
    Dimensionality                  13
    Features            real, positive
    =================   ==============

    Read more in the :ref:`User Guide <wine_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are: 'data', the
        data to learn, 'target', the classification labels, 'target_names', the
        meaning of the labels, 'feature_names', the meaning of the features,
        and 'DESCR', the full description of the dataset.

    (data, target) : tuple if ``return_X_y`` is True

    The copy of UCI ML Wine Data Set dataset is downloaded and modified to fit
    standard format from:
    https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data

    Examples
    --------
    Let's say you are interested in the samples 10, 80, and 140, and want to
    know their class name.

    >>> from sklearn.datasets import Wine
    >>> data = Wine().load()
    >>> data.target[[10, 80, 140]]
    array([0, 1, 2])
    >>> list(data.target_names)
    ['class_0', 'class_1', 'class_2']
    """
    return Wine().load(return_X_y=return_X_y)


@deprecated("'load_iris' was renamed to 'Iris().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_iris(return_X_y=False):
    """Load and return the iris dataset (classification).

    The iris dataset is a classic and very easy multi-class classification
    dataset.

    =================   ==============
    Classes                          3
    Samples per class               50
    Samples total                  150
    Dimensionality                   4
    Features            real, positive
    =================   ==============

    Read more in the :ref:`User Guide <iris_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object. See
        below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        'target_names', the meaning of the labels, 'feature_names', the
        meaning of the features, 'DESCR', the full description of
        the dataset, 'filename', the physical location of
        iris csv dataset (added in version `0.20`).

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18

    Notes
    -----
        .. versionchanged:: 0.20
            Fixed two wrong data points according to Fisher's paper.
            The new version is the same as in R, but not as in the UCI
            Machine Learning Repository.

    Examples
    --------
    Let's say you are interested in the samples 10, 25, and 50, and want to
    know their class name.

    >>> from sklearn.datasets import Iris
    >>> data = Iris().load()
    >>> data.target[[10, 25, 50]]
    array([0, 0, 1])
    >>> list(data.target_names)
    ['setosa', 'versicolor', 'virginica']
    """
    return Iris().load(return_X_y=return_X_y)


@deprecated("'load_breast_cancer' was renamed to 'BreastCancer().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_breast_cancer(return_X_y=False):
    """Load and return the breast cancer wisconsin dataset (classification).

    The breast cancer dataset is a classic and very easy binary classification
    dataset.

    =================   ==============
    Classes                          2
    Samples per class    212(M),357(B)
    Samples total                  569
    Dimensionality                  30
    Features            real, positive
    =================   ==============

    Read more in the :ref:`User Guide <breast_cancer_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        'target_names', the meaning of the labels, 'feature_names', the
        meaning of the features, and 'DESCR', the full description of
        the dataset, 'filename', the physical location of
        breast cancer csv dataset (added in version `0.20`).

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18

    The copy of UCI ML Breast Cancer Wisconsin (Diagnostic) dataset is
    downloaded from:
    https://goo.gl/U2Uwz2

    Examples
    --------
    Let's say you are interested in the samples 10, 50, and 85, and want to
    know their class name.

    >>> from sklearn.datasets import BreastCancer
    >>> data = BreastCancer().load()
    >>> data.target[[10, 50, 85]]
    array([0, 1, 0])
    >>> list(data.target_names)
    ['malignant', 'benign']
    """
    return BreastCancer().load(return_X_y=return_X_y)


@deprecated("'load_digits' was renamed to 'Digits().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_digits(n_class=10, return_X_y=False):
    """Load and return the digits dataset (classification).

    Each datapoint is a 8x8 image of a digit.

    =================   ==============
    Classes                         10
    Samples per class             ~180
    Samples total                 1797
    Dimensionality                  64
    Features             integers 0-16
    =================   ==============

    Read more in the :ref:`User Guide <digits_dataset>`.

    Parameters
    ----------
    n_class : integer, between 0 and 10, optional (default=10)
        The number of classes to return.

    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'images', the images corresponding
        to each sample, 'target', the classification labels for each
        sample, 'target_names', the meaning of the labels, and 'DESCR',
        the full description of the dataset.

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18

    This is a copy of the test set of the UCI ML hand-written digits datasets
    https://archive.ics.uci.edu/ml/datasets/Optical+Recognition+of+Handwritten+Digits

    Examples
    --------
    To load the data and visualize the images::

        >>> from sklearn.datasets import Digits
        >>> digits = Digits().load()
        >>> print(digits.data.shape)
        (1797, 64)
        >>> import matplotlib.pyplot as plt #doctest: +SKIP
        >>> plt.gray() #doctest: +SKIP
        >>> plt.matshow(digits.images[0]) #doctest: +SKIP
        >>> plt.show() #doctest: +SKIP
    """
    return Digits(n_class=n_class).load(return_X_y=return_X_y)


@deprecated("'load_diabetes' was renamed to 'Diabetes().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_diabetes(return_X_y=False):
    """Load and return the diabetes dataset (regression).

    ==============      ==================
    Samples total       442
    Dimensionality      10
    Features            real, -.2 < x < .2
    Targets             integer 25 - 346
    ==============      ==================

    Read more in the :ref:`User Guide <diabetes_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the regression target for each
        sample, 'data_filename', the physical location
        of diabetes data csv dataset, and 'target_filename', the physical
        location of diabetes targets csv datataset (added in version `0.20`).

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18
    """
    return Diabetes().load(return_X_y=return_X_y)


@deprecated("'load_linnerud' was renamed to 'Linnerud().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_linnerud(return_X_y=False):
    """Load and return the linnerud dataset (multivariate regression).

    ==============    ============================
    Samples total     20
    Dimensionality    3 (for both data and target)
    Features          integer
    Targets           integer
    ==============    ============================

    Read more in the :ref:`User Guide <linnerrud_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are: 'data' and
        'targets', the two multivariate datasets, with 'data' corresponding to
        the exercise and 'targets' corresponding to the physiological
        measurements, as well as 'feature_names' and 'target_names'.
        In addition, you will also have access to 'data_filename',
        the physical location of linnerud data csv dataset, and
        'target_filename', the physical location of
        linnerud targets csv datataset (added in version `0.20`).

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18
    """
    return Linnerud().load(return_X_y=return_X_y)


@deprecated("'load_boston' was renamed to 'Boston().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_boston(return_X_y=False):
    """Load and return the boston house-prices dataset (regression).

    ==============     ==============
    Samples total                 506
    Dimensionality                 13
    Features           real, positive
    Targets             real 5. - 50.
    ==============     ==============

    Read more in the :ref:`User Guide <boston_dataset>`.

    Parameters
    ----------
    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object.
        See below for more information about the `data` and `target` object.

        .. versionadded:: 0.18

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the regression targets,
        'DESCR', the full description of the dataset,
        and 'filename', the physical location of boston
        csv dataset (added in version `0.20`).

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.18

    Notes
    -----
        .. versionchanged:: 0.20
            Fixed a wrong data point at [445, 0].

    Examples
    --------
    >>> from sklearn.datasets import Boston
    >>> boston = Boston().load()
    >>> print(boston.data.shape)
    (506, 13)
    """
    return Boston().load(return_X_y=return_X_y)


@deprecated("'load_sample_images' was renamed to 'SampleImages().load'"
            "in version 0.21 and will be removed in 0.23.")
def load_sample_images():
    """Load sample images for image manipulation.

    Loads both, ``china`` and ``flower``.

    Read more in the :ref:`User Guide <sample_images>`.

    Returns
    -------
    data : Bunch
        Dictionary-like object with the following attributes : 'images', the
        two sample images, 'filenames', the file names for the images, and
        'DESCR' the full description of the dataset.

    Examples
    --------
    To load the data and visualize the images:

    >>> from sklearn.datasets import SampleImages
    >>> dataset = SampleImages().load()     #doctest: +SKIP
    >>> len(dataset.images)                #doctest: +SKIP
    2
    >>> first_img_data = dataset.images[0] #doctest: +SKIP
    >>> first_img_data.shape               #doctest: +SKIP
    (427, 640, 3)
    >>> first_img_data.dtype               #doctest: +SKIP
    dtype('uint8')
    """
    return SampleImages().load()


@deprecated("'load_sample_image' was refactored"
            "in version 0.21 and will be removed in 0.23."
            "Use 'SampleImages().load' instead")
def load_sample_image(image_name):
    """Load the numpy array of a single sample image

    Read more in the :ref:`User Guide <sample_images>`.

    Parameters
    -----------
    image_name : {`china.jpg`, `flower.jpg`}
        The name of the sample image loaded

    Returns
    -------
    img : 3D array
        The image as a numpy array: height x width x color

    Examples
    ---------
    >>> from sklearn.datasets import SampleImages
    >>> china = SampleImages('china.jpg').load().images[0]   # doctest: +SKIP
    >>> china.dtype                              # doctest: +SKIP
    dtype('uint8')
    >>> china.shape                              # doctest: +SKIP
    (427, 640, 3)
    >>> flower = SampleImages('flower.jpg').load().images[0] # doctest: +SKIP
    >>> flower.dtype                             # doctest: +SKIP
    dtype('uint8')
    >>> flower.shape                             # doctest: +SKIP
    (427, 640, 3)
    """
    return SampleImages(image_name).load().images[0]


class DatasetLoader(object, metaclass=ABCMeta):
    """Abstract class for all dataset loaders in scikit-learn."""

    def load(self, return_X_y=False) -> Dataset:
        bunch = self._raw_data_to_bunch()
        if return_X_y:
            return bunch.data, bunch.target
        return bunch

    @abstractmethod
    def _raw_data_to_bunch(self) -> Bunch:
        raise NotImplementedError


class LocalDatasetLoader(DatasetLoader):
    _module_path = dirname(__file__)
    _data_dir = join(_module_path, 'data')
    _descr_dir = join(_module_path, 'descr')
    _images_dir = join(_module_path, 'images')

    @property
    def X_file(self) -> Union[str, ndarray]:
        raise NotImplementedError

    @property
    def y_file(self) -> str:
        return self.X_file

    @property
    def descr_file(self):
        return self.X_file.split('.', maxsplit=1)[0] + '.rst'

    @property
    def local_data_paths(self) -> Dict[str, str]:
        try:
            return {
                'X': join(self._data_dir, self.X_file),
                'y': join(self._data_dir, self.y_file),
                'descr': join(self._descr_dir, self.descr_file),
            }
        except TypeError:
            return {
                'X': np.array([join(self._data_dir, f) for f in self.X_file]),
                'y': np.array([join(self._data_dir, f) for f in self.y_file]),
                'descr': join(self._descr_dir, self.descr_file),
            }

    @property
    def feature_names(self) -> ndarray:
        raise NotImplementedError

    @property
    def target_names(self) -> ndarray:
        raise NotImplementedError

    _attempt_cast_to_int = staticmethod(_attempt_cast_to_int)

    def _raw_data_to_bunch(self) -> Bunch:
        bunch = self.read_data()
        bunch.DESCR = self._read_description()
        return self.process(bunch)

    @abstractmethod
    def read_data(self) -> Bunch:
        raise NotImplementedError

    @abstractmethod
    def process(self, bunch: Bunch) -> Bunch:
        bunch.target = self._attempt_cast_to_int(bunch.target)
        return bunch

    def _read_description(self, path: str = None) -> str:
        descr_path = path or self.local_data_paths['descr']
        with open(descr_path) as descr:
            return descr.read()

    def _make_bunch(self, X, y, target_names,
                    description, images=None) -> Bunch:
        return Bunch(
            data=X, target=y,
            feature_names=self.feature_names,
            target_names=target_names,
            DESCR=description,
            images=images,
            data_filename=self.local_data_paths['X'],
            target_filename=self.local_data_paths['y'],
            filename=self.local_data_paths['X'])  # Backwards compatible (0.21)


class SimpleCSVLocalDatasetLoader(LocalDatasetLoader):
    """Reads a .csv file with:
        The first row containing:
            - the sample size (int)
            - the number of features (int)
            - [Optional] the target names corresponding
              to their 'int' representation in the y column
        The second [Optional] row containing:
            - feature variable names with type 'str'
        X features in columns [:-1] with type 'int' or 'float'
        y target in column [-1] with type 'int' or 'float'

    If you try to read a .csv file containing 'object' or
    'str' type variable values, it will fail!
     """
    def read_data(self) -> Bunch:
        X, y, target_names = self.read_X_y_csv(self.local_data_paths['X'])
        if self.target_names.size > 0:  # class property takes precedence
            target_names = self.target_names
        return self._make_bunch(X, y, target_names, None, None)

    def process(self, bunch: Bunch) -> Bunch:
        return super().process(bunch)

    @staticmethod
    def read_X_y_csv(path: str) -> Tuple[ndarray, ndarray, ndarray]:
        with open(path) as f:
            firstline = f.readline().rstrip().split(',')
            n_features = int(firstline[1])
            col_ixs = tuple(range(n_features + 1))
            target_names = np.array(firstline[2:])

        csv_arr = np.genfromtxt(path, delimiter=',', usecols=col_ixs,
                                skip_header=1, dtype=np.float)
        data, target = csv_arr[:, :n_features], csv_arr[:, n_features]

        mask_row_all_nan = ~np.isnan(data).all(axis=1)
        data, target = data[mask_row_all_nan], target[mask_row_all_nan]
        return data, target, target_names


class Wine(SimpleCSVLocalDatasetLoader):
    X_file = 'wine_data.csv'
    feature_names = np.array([
        'alcohol', 'malic_acid',
        'ash', 'alcalinity_of_ash',
        'magnesium', 'total_phenols',
        'flavanods', 'nonflavanoid_phenols',
        'proanthocyanins', 'color_intensity',
        'hue', 'od280/od315_of_diluted_wines',
        'proline'])
    target_names = np.array([])


class Iris(SimpleCSVLocalDatasetLoader):
    X_file = 'iris.csv'
    feature_names = np.array([
        'sepal length (cm)', 'sepal width (cm)',
        'petal length (cm)', 'petal width (cm)'])
    target_names = np.array([])


class BreastCancer(SimpleCSVLocalDatasetLoader):
    X_file = 'breast_cancer.csv'
    feature_names = np.array([
        'mean radius', 'mean texture',
        'mean perimeter', 'mean area',
        'mean smoothness', 'mean compactness',
        'mean concavity', 'mean concave points',
        'mean symmetry', 'mean fractal dimension',
        'radius error', 'texture error',
        'perimeter error', 'area error',
        'smoothness error', 'compactness error',
        'concavity error', 'concave points error',
        'symmetry error', 'fractal dimension error',
        'worst radius', 'worst texture',
        'worst perimeter', 'worst area',
        'worst smoothness', 'worst compactness',
        'worst concavity', 'worst concave points',
        'worst symmetry', 'worst frctal dimension'])
    target_names = np.array([])


class Digits(LocalDatasetLoader):
    X_file = 'digits.csv.gz'
    feature_names = np.array([list(map(lambda d: (d, x).__repr__(),
                                       np.arange(1, 9, dtype='int')))
                              for x in np.arange(1, 9, dtype='int')])
    target_names = np.arange(10)

    def __init__(self, n_class=10):
        self.n_class = n_class

    def read_data(self):
        X = np.loadtxt(self.local_data_paths['X'], delimiter=',')
        y = X[:, -1].astype('int')
        return self._make_bunch(X, y, self.target_names, None, None)

    def process(self, bunch: Bunch) -> Bunch:
        X, y = bunch.data, bunch.target
        flat_X = X[:, :-1]
        images = flat_X.view()
        images.shape = (-1, 8, 8)
        if self.n_class < 10:
            idx = y < self.n_class
            flat_X, y = flat_X[idx], y[idx]
            images = images[idx]
        bunch.data, bunch.target, bunch.images = flat_X,  y, images
        return bunch


class Diabetes(LocalDatasetLoader):
    X_file = 'diabetes_data.csv.gz'
    y_file = 'diabetes_target.csv.gz'
    descr_file = 'diabetes.rst'
    feature_names = np.array([
        'age', 'sex', 'bmi', 'bp',
        's1', 's2', 's3', 's4', 's5', 's6'])
    target_names = np.array(['progression'])

    def read_data(self):
        X = np.loadtxt(self.local_data_paths['X'], dtype='float')
        y = np.loadtxt(self.local_data_paths['y'], dtype='float')
        return self._make_bunch(X, y, self.target_names, None, None)

    def process(self, bunch: Bunch):
        return super().process(bunch)


class Linnerud(LocalDatasetLoader):
    X_file = 'linnerud_exercise.csv'
    y_file = 'linnerud_physiological.csv'
    descr_file = 'linnerud.rst'
    feature_names = np.array(['Chins', 'Situps', 'Jumps'])
    target_names = np.array(['Weight', 'Waist', 'Pulse'])

    def read_data(self) -> Bunch:
        X = np.loadtxt(self.local_data_paths['X'], skiprows=1, dtype='float')
        y = np.loadtxt(self.local_data_paths['y'], skiprows=1, dtype='float')
        return self._make_bunch(X, y, self.target_names, None, None)

    def process(self, bunch: Bunch) -> Bunch:
        return bunch


class Boston(SimpleCSVLocalDatasetLoader):
    X_file = 'boston_house_prices.csv'
    feature_names = np.array([
        'CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE',
        'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT'])
    target_names = np.array(['MEDV'])


class SampleImages(LocalDatasetLoader):
    X_file = np.array(['china.jpg', 'flower.jpg'])
    feature_names = np.array([])
    target_names = np.array(['china', 'flower'])
    _descr_dir = LocalDatasetLoader._images_dir
    descr_file = 'README.txt'

    def __init__(self, image_name: str = None):
        self.image_name = image_name
        self._check_image_name()

    def _check_image_name(self):
        if self.image_name and self.image_name not in self.X_file:
            msg = 'Cannot find sample image: %s' % self.image_name
            raise AttributeError(msg)

    def read_data(self) -> Bunch:
        from ..externals._pilutil import imread  # import PIL only when needed
        image_files = [join(self._images_dir, file) for file in self.X_file]
        images = [imread(img) for img in image_files]
        return self._make_bunch(None, None, self.target_names, None, images)

    def process(self, bunch: Bunch):
        bunch.filenames = bunch.filename  # Backwards compatible (0.21)
        return bunch


def _pkl_filepath(*args, **kwargs):
    """Ensure different filenames for Python 2 and Python 3 pickles

    An object pickled under Python 3 cannot be loaded under Python 2. An object
    pickled under Python 2 can sometimes not be loaded correctly under Python 3
    because some Python 2 strings are decoded as Python 3 strings which can be
    problematic for objects that use Python 2 strings as byte buffers for
    numerical data instead of "real" strings.

    Therefore, dataset loaders in scikit-learn use different files for pickles
    manages by Python 2 and Python 3 in the same SCIKIT_LEARN_DATA folder so as
    to avoid conflicts.

    args[-1] is expected to be the ".pkl" filename. Under Python 3, a suffix is
    inserted before the extension to s

    _pkl_filepath('/path/to/folder', 'filename.pkl') returns:
      - /path/to/folder/filename.pkl under Python 2
      - /path/to/folder/filename_py3.pkl under Python 3+

    """
    py3_suffix = kwargs.get("py3_suffix", "_py3")
    basename, ext = splitext(args[-1])
    if sys.version_info[0] >= 3:
        basename += py3_suffix
    new_args = args[:-1] + (basename + ext,)
    return join(*new_args)


def _sha256(path):
    """Calculate the sha256 hash of the file at path."""
    sha256hash = hashlib.sha256()
    chunk_size = 8192
    with open(path, "rb") as f:
        while True:
            buffer = f.read(chunk_size)
            if not buffer:
                break
            sha256hash.update(buffer)
    return sha256hash.hexdigest()


def _fetch_remote(remote, dirname=None):
    """Helper function to download a remote dataset into path

    Fetch a dataset pointed by remote's url, save into path using remote's
    filename and ensure its integrity based on the SHA256 Checksum of the
    downloaded file.

    Parameters
    -----------
    remote : RemoteFileMetadata
        Named tuple containing remote dataset meta information: url, filename
        and checksum

    dirname : string
        Directory to save the file to.

    Returns
    -------
    file_path: string
        Full path of the created file.
    """

    file_path = (remote.filename if dirname is None
                 else join(dirname, remote.filename))
    urlretrieve(remote.url, file_path)
    checksum = _sha256(file_path)
    if remote.checksum != checksum:
        raise IOError("{} has an SHA256 checksum ({}) "
                      "differing from expected ({}), "
                      "file may be corrupted.".format(file_path, checksum,
                                                      remote.checksum))
    return file_path
