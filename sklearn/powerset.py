# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from warnings import warn

from .basic_checks import check_is_in_range

__all__ = [
    "Powerset",
]


def _get(dic):
    def _t(val):
        try:
            v = dic[val]
        except KeyError:
            raise ValueError('Error: The key {0} is invalid.'.format(val))
        return v
    return _t


class Powerset:

    """
    This class allows to compress or powerize data following
    several fashions.

    Currently, all values from the data must either be 0 or 1.
    """

    def __init__(self):
        self.dict = {}
        self.tcid = {}
        self.nb_cols = 1
        self.uni = False
        self.power_method = self.__NONE

    def __reset(self, array):
        self.dict = {}
        self.tcid = {}
        self.uni = False
        self.power_method = self.__NONE
        try:
            self.nb_cols = np.asanyarray(array).shape[1]
        except IndexError:
            self.nb_cols = 1
            self.uni = True

    def __reshape(self, array):
        if self.uni:
            return np.reshape(array, (array.shape[0],))
        return array

    __AMPLIFIED = 'amplified'
    __COMPRESS_STRING = 'compress_string'
    __COMPRESS = 'compress'
    __DUMMY = 'dummy'
    __NONE = None

    def amplified_powerize(self, array):
        """
        Powerizes the (binary) data in a simple way.
        If an overflow occurs, it will automatically call the compressed
        version (compressed_powerize).
        Such an overflow can occur if the number of columns of the argument
        is larger than the number of unsigned bits in an integer.

        Parameters
        ----------
        array : numpy.array
            The array to powerize.

        Returns
        -------
        numpy.ndarray of shape = array.shape[0] and dtype = int
            The powerized array.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset().amplified_powerize(numpy.array(array))
        [5 3 0]
        """
        self.__reset(array)
        universal = self.universal_powerize(array)
        convert_to_int = np.vectorize(lambda x: int(x, base=2))
        try:
            converted = convert_to_int(universal)
            converted = converted.astype(np.int, copy=False)
            if np.any(converted < 0):
                raise OverflowError('Overflowed possible at decompress stage.')
        except OverflowError as e:
            warn('{0} -> Forced to use the compressed version.'.format(str(e)))
            return self.compressed_powerize(array)
        self.power_method = self.__AMPLIFIED
        return converted

    def compressed_string_powerize(self, array):
        """
        Powerizes the (binary) data in a compressed way.
        The data is compressed in the sense that the integers
        used are the smallest possible, thus reducing the possibility
        of an overflow.

        Parameters
        ----------
        array : numpy.array
            The array to powerize.

        Returns
        -------
        numpy.ndarray of shape = array.shape[0] and dtype = int
            The powerized array.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset().compressed_string_powerize(numpy.array(array))
        [0 1 2]
        """
        self.__reset(array)
        universal = self.universal_powerize(array)
        _, unique_indexes = np.unique(universal, return_index=True)
        unique_indexes.sort()
        for idx, each in enumerate(unique_indexes):
            self.dict.update({universal[each]: idx})
            self.tcid.update({idx: universal[each]})
        zip_fct = np.vectorize(_get(self.dict))
        self.power_method = self.__COMPRESS_STRING
        return zip_fct(universal)

    def compressed_powerize(self, array):
        """
        Powerizes the (binary) data in a compressed way.
        The data is compressed in the sense that the integers
        used are the smallest possible, thus reducing the possibility
        of an overflow.

        Parameters
        ----------
        array : numpy.array
            The array to powerize.

        Returns
        -------
        numpy.ndarray of shape = array.shape[0] and dtype = int
            The powerized array.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset().compressed_powerize(numpy.array(array))
        [0 1 2]
        """
        self.__reset(array)
        universal = self.universal_powerize(array)
        _, unique_indexes = np.unique(universal, return_index=True)
        unique_indexes.sort()
        for idx, each in enumerate(unique_indexes):
            self.dict.update({long(universal[each], base=2): idx})
            self.tcid.update({idx: long(universal[each], base=2)})
        _tmp = (lambda x: _get(self.dict)(long(x, base=2)))
        zip_fct = np.vectorize(_tmp)
        self.power_method = self.__COMPRESS
        return zip_fct(universal)

    def dummy_powerize(self, array):
        """
        Does not powerize the data.
        It can be used transparently.

        Parameters
        ----------
        array : numpy.array
            The array to powerize.

        Returns
        -------
        numpy.array
            The array argument.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset().dummy_powerize(numpy.array(array))
        [[1, 0, 1], [0, 1, 1], [0, 0, 0]]
        """
        self.__reset(array)
        self.power_method = self.__DUMMY
        return array

    def __amplified_unpowerize(self, array):
        if (self.power_method != self.__AMPLIFIED):
            raise ValueError('Wrong unpowerizing method!')
        convert_to_string = np.vectorize(
            lambda x: (("{:0" + str(self.nb_cols) + "b}").format(x)))
        stringified = convert_to_string(
            np.asanyarray(array).astype(object, copy=False))
        return self.__reshape(self.universal_unpowerize(stringified))

    def __compressed_string_unpowerize(self, array):
        if (self.power_method != self.__COMPRESS_STRING):
            raise ValueError('Wrong unpowerizing method!')
        unzip_fct = np.vectorize(_get(self.tcid))
        return self.__reshape(self.universal_unpowerize(unzip_fct(array)))

    def __compressed_unpowerize(self, array):
        if (self.power_method != self.__COMPRESS):
            raise ValueError('Wrong unpowerizing method!')

        def x(y):
            r = _get(self.tcid)(y)
            r = ("{:0" + str(self.nb_cols) + "b}").format(r)
            return r
        unzip_fct = np.vectorize(x)
        return self.__reshape(self.universal_unpowerize(unzip_fct(array)))

    def __dummy_unpowerize(self, array):
        if (self.power_method != self.__DUMMY):
            raise ValueError('Wrong unpowerizing method!')
        return self.__reshape(array)

    @staticmethod
    def _probas_aux(uncompressed, cls_probas_array, n_odds=2):
        """
        uncompressed = mapping from idx to line
        (1 -> [010010] for example 6 labels) must have key 0

        cls_proba_array = 1 line of result of predict_proba
        ([0.2 0.1 0.8 0.6] if 4 classes)

        Returns 1 line with correct probas

        Example:
        uncompressed = { 0 : [0,1,0,0], 1 : [0,0,1,1], 2 : [1,1,0,0]}
        uncompressed = { k : [wxyz]}
        cls_proba_array = [0.4, 0.5, 0.1]
        cls_proba_array = [P(0), P(1), P(2)] = [P(k)]
        >>> list([0.4 + 0.5, 0.1],
                 [0.5, 0.4 + 0.1],
                 [0.4 + 0.1, 0.5],
                 [0.4 + 0.1, 0.5])
        >>> # i.e. np.array([[P(dict(0)[0] = 0) * P(0), [], [], []])
        sum(P(dict(i)[0] = 0) * P(i))
        """
        n_labels = len(uncompressed[0])
        retv = [np.zeros((n_odds,), dtype=float) for _ in xrange(n_labels)]
        for label, odd in enumerate(retv):
            for cls, cls_odd in enumerate(cls_probas_array):
                odd[uncompressed[cls][label]] += cls_odd
        return retv

    def probas_unpowerize(self, probas_array):
        """
        Takes the output of predict_proba and gives the probabilities
        for each original class.

        It decompresses the classes, thus a compression
        method should first have been selected and used to
        fit the estimator.

        Parameters
        ----------
        probas_array : array-like
            Output of predict_proba of a classifier.

            The array must be a list of size n_labels containing
            numpy.ndarray of floats with a size of (n_samples, n_classes).

        Returns
        -------
        list of numpy.ndarray
            This list has the same format as `probas_arrays`
            but contains the decompressed probabilities.

        Example
        -------
        >>> powerset = Powerset()
        >>> p = powerset.amplified_powerize(np.array(
                [[0, 0, 0, ], [0, 0, 1, ], [0, 1, 0, ], [0, 1, 1, ],
                 [1, 0, 0, ], [1, 0, 1, ], [1, 1, 0, ], [1, 1, 1, ]]))
        [0, 1, 2, 3, 4, 5, 6, 7]
        >>> c = np.array([[0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16]])
        >>> powerset.probas_unpowerize(c)
        [array([[ 0.42,  0.58]]),
         array([[ 0.46,  0.54]]),
         array([[ 0.48,  0.52]])]
        """
        methods = {self.__AMPLIFIED: None,
                   self.__COMPRESS_STRING: None,
                   self.__COMPRESS: None,
                   self.__DUMMY: True,
                   self.__NONE: False
                   }
        tmp = methods.get(self.power_method, False)
        if tmp is None:
            pass
        elif tmp:
            return probas_array
        else:
            raise ValueError('No powerization already.')
        try:
            n_samples, n_classes = probas_array.shape
        except (AttributeError, ValueError):
            n_samples, n_classes = probas_array[0].shape
            used_array = probas_array
        else:
            used_array = [probas_array]

        uncompressed = {}
        for i in xrange(n_classes):
            tmp = self.unpowerize(np.array([i]))[0]
            if self.uni:
                tmp = [tmp]
            uncompressed.update({i: tmp})

        n_odds = 2

        retv = [np.empty((n_samples, n_odds), dtype=float)
                for _ in xrange(self.nb_cols)]
        for idx_label, label in enumerate(used_array):
            for idx_sample, sample in enumerate(label):
                tmp = self._probas_aux(uncompressed, sample, n_odds=n_odds)
                for idx_int_label, int_label in enumerate(tmp):
                    retv[idx_int_label][idx_sample] = int_label

        if self.uni:
            return retv[0]
        return retv

    def majority_unpowerize(self, probas_array):
        """
        Applies the reverse function of the powerization but only
        gives back the highest probable class from a probability array.

        Parameters
        ----------
        probas_array : array-like
            Output of predict_proba of a classifier.

            The array must be a list of size n_labels containing
            numpy.ndarray of floats with a size of (n_samples, n_classes).

        Returns
        -------
        list of numpy.array
            This list has a length of the number of labels, and each
            array a length of the number of samples.

            Each value corresponds to the class of the sample for the
            corresponding label.

        Example
        -------
        >>> powerset = Powerset()
        >>> p = powerset.amplified_powerize(np.array(
                [[0, 0, 0, ], [0, 0, 1, ], [0, 1, 0, ], [0, 1, 1, ],
                 [1, 0, 0, ], [1, 0, 1, ], [1, 1, 0, ], [1, 1, 1, ]]))
        [0, 1, 2, 3, 4, 5, 6, 7]
        >>> c = np.array([[0.09, 0.10, 0.11, 0.12, 0.13, 0.16, 0.15, 0.14],
                          [0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16]])
        >>> powerset.majority_unpowerize(c)
        [array([1, 1]), array([0, 1]), array([1, 1])]
        """
        methods = {self.__AMPLIFIED: None,
                   self.__COMPRESS_STRING: None,
                   self.__COMPRESS: None,
                   self.__DUMMY: True,
                   self.__NONE: False
                   }
        tmp = methods.get(self.power_method, False)
        if tmp is None:
            pass
        elif not tmp:
            raise ValueError('No powerization already.')
        try:
            n_samples, n_classes = probas_array.shape
        except (AttributeError, ValueError):
            n_samples, n_classes = probas_array[0].shape
            used_array = probas_array
        else:
            used_array = [probas_array]

        if len(used_array) != 1:
            raise ValueError('Currently only supports single label.')
        uncompressed = {}
        for i in xrange(n_classes):
            tmp = self.unpowerize(np.array([i]))[0]
            if self.uni:
                tmp = [tmp]
            uncompressed.update({i: tmp})

        retv = [np.array([0 for _ in xrange(n_samples)])
                for i in xrange(self.nb_cols)]
        for idx, b in enumerate(used_array):
            r = np.argmax(b, axis=1)
            for idx_sample, sample in enumerate(r):
                l = uncompressed[sample]
                for x in xrange(len(l)):
                    retv[x][idx_sample] = l[x]

        if self.uni:
            return retv[0]
        return retv

    def unpowerize(self, array):
        """
        Applies the reverse function of the powerization previously used.

        Parameters
        ----------
        array : numpy.array
            The array to unpowerize.

        Returns
        -------
        numpy.ndarray of shape = (array.shape[0], powerized_array.shape[1])
        and dtype = int
            The unpowerized array.

        Raises
        ------
        ValueError
            If the Powerset object did not powerize first.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> p = Powerset()
        >>> p.unpowerize(p.compressed_powerize(numpy.array(array)))
        [[1, 0, 1], [0, 1, 1], [0, 0, 0]]
        """
        methods = {self.__AMPLIFIED: self.__amplified_unpowerize,
                   self.__COMPRESS_STRING: self.__compressed_string_unpowerize,
                   self.__COMPRESS: self.__compressed_unpowerize,
                   self.__DUMMY: self.__dummy_unpowerize,
                   self.__NONE: None
                   }
        tmp = methods.get(self.power_method)
        if (tmp is None):
            raise ValueError('No powerization already.')
        return tmp(array)

    @staticmethod
    def universal_powerize(array):
        """
        Powerizes data into string.

        Parameters
        ----------
        array : numpy.array
            The array to powerize.

        Returns
        -------
        numpy.ndarray of shape = array.shape[0] and dtype = object
            The powerized array.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset.universal_powerize(numpy.array(array))
        ["101", "011", "000"]
        """
        array_ = np.asanyarray(array)
        check_is_in_range(len(array_.shape), lower_bound=1, higher_bound=2,
                          low_exclusive=False, high_exclusive=False,
                          launch_exception=True,
                          msg=("Only uni-bidimensional arrays are supported ("
                               "-D passed).").format(str(len(array_.shape))),
                          exception_type=ValueError)
        result = np.empty((array_.shape[0], ), dtype=object)
        try:
            view = np.reshape(array_, (array_.shape[0], array_.shape[1]))
        except IndexError:
            view = np.reshape(array_, (array_.shape[0], 1))
        uniques = np.unique(view)

        check_is_in_range(elem=len(uniques), lower_bound=1, higher_bound=2,
                          low_exclusive=False, high_exclusive=False,
                          launch_exception=True,
                          msg="The array does not contain 2 classes or less.",
                          exception_type=ValueError)
        if not (((len(uniques) == 2) and (0 in uniques) and (1 in uniques)) or
                ((len(uniques) == 1) and ((0 in uniques) or (1 in uniques)))):
            raise ValueError("Only binary arrays (0 and 1) are supported.")

        for line in xrange(0, view.shape[0]):
            l = ''
            for elem in view[line, :]:
                l += str(elem)
            result[line] = l

        return result

    @staticmethod
    def universal_unpowerize(array):
        """
        Applies the reverse function of universal_powerize.

        Parameters
        ----------
        array : numpy.array
            The array to unpowerize.

        Returns
        -------
        numpy.ndarray of shape = (array.shape[0], powerized_array.shape[1])
        and dtype = int
            The unpowerized array.

        Raises
        ------
        ValueError
            If the array is not unidimensional.

        Example
        -------
        >>> array = [[1, 0, 1],
                     [0, 1, 1],
                     [0, 0, 0]]
        >>> Powerset.universal_unpowerize(
        Powerset.universal_powerize(numpy.array(array)))
        [[1, 0, 1], [0, 1, 1], [0, 0, 0]]
        """
        array_ = np.asanyarray(array)
        check_is_in_range(len(array_.shape), higher_bound=2,
                          launch_exception=True,
                          msg=("Only unidimensional arrays are supported "
                               "({0}-D passed).").format(
                                   str(len(array_.shape))),
                          exception_type=ValueError)
        view = np.reshape(array_, (array_.shape[0],))
        result = np.zeros((view.shape[0], len(view[0])), dtype=int)

        for line in xrange(result.shape[0]):
            for col in xrange(result.shape[1]):
                result[line, col] = int(view[line][col], base=2)

        return result
