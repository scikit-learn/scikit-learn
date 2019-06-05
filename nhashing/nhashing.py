import pandas as pd
import multiprocessing
from category_encoders.hashing import HashingEncoder

__author__ = 'LiuShulun'


class NHashingEncoder:
    """A extended multivariate hashing implementation using multi-process

    !!! Significant !!! improve the encoding speed of hashing encoder on
    multi-kernel device. This NHashingEncoder allows setting
    max_process & max_sample to use more CPU resource when hashing encoding.
    On a 4C8T CPU, NHashingEncoder with max_process=4 makes
    encoding 3+ times faster than hashingEncoder.

    Parameters
    ----------

    verbose: int
        integer indicating verbosity of the output. 0 for none.
    cols: list
        a list of columns to encode, if None all string columns'll be encoded.
    drop_invariant: bool
        boolean for whether or not to drop columns with 0 variance.
    return_df: bool
        boolean for whether to return a pandas DataFrame from transform
        (otherwise it will be a numpy array).
    hash_method: str
        which hashing method to use. Any method from hashlib works.
    max_process: int
        how many PROCESS to use, NOT Thread. limited in range(1, 64)
    max_sample: int
        how many samples will be encode by each process at each time

    Example
    -------
    >>> from category_encoders.hashing import HashingEncoder
    >>> from nhashing import NHashingEncoder
    >>> import time
    >>> import pandas as pd
    >>> from sklearn.datasets import load_boston
    >>> bunch = load_boston()
    >>> X = pd.DataFrame(bunch.data, columns=bunch.feature_names)
    >>> DL = []
    >>> for i in range(1000): DL.append(X)
    >>> DF = pd.concat(DL, ignore_index=True).reset_index(drop=True)
    >>> he = HashingEncoder(cols=['CHAS', 'RAD'])
    >>> start = time.time()
    >>> he.fit_transform(DF)
    >>> he_time = time.time() - start
    >>> nhe = NHashingEncoder(cols=['CHAS', 'RAD'])
    >>> start = time.time()
    >>> nhe.fit_transform(DF)
    >>> nhe_time = time.time() - start
    >>> print("500000samples HashingEncoder Time:", he_time, "NHashingEncoder Time:", nhe_time)

    500000samples HashingEncoder Time: 358.14579820632935 NHashingEncoder Time: 110.51188397407532

    References
    ----------
    .. [1] Feature Hashing for Large Scale Multitask Learning, from
    https://alex.smola.org/papers/2009/Weinbergeretal09.pdf

    """
    def __init__(self, verbose=0, n_components=8, cols=None, drop_invariant=False,
                 return_df=True, hash_method='md5', max_process=1, max_sample=0):
        self.verbose = verbose
        self.n_components = n_components
        self.cols = cols
        self.drop_invariant = drop_invariant
        self.return_df = return_df
        self.hash_method = hash_method
        self.max_process = 1 if max_process not in range(1, 64) else max_process
        self.max_sample = max_sample
        self.data_lock = multiprocessing.Lock()
        self.start_state = multiprocessing.Manager().Queue()
        self.start_state.put(-1)
        self.origin_parts = multiprocessing.Manager().Queue()
        self.hashing_parts = multiprocessing.Manager().Queue()
        self.n_process = []
        self.data_lines = 0
        self.data = None

    def __hash_encode(self, data_part, cols):
        he = HashingEncoder(verbose=self.verbose, n_components=self.n_components,
                            cols=cols, drop_invariant=self.drop_invariant,
                            return_df=self.return_df, hash_method=self.hash_method)
        data_part = he.fit_transform(data_part)

        return data_part

    def __require_data(self, cols, process_index):
        if self.data_lock.acquire():
            if not self.start_state.empty():
                done_index = 0
                while not self.start_state.empty():
                    self.start_state.get()
            else:
                if self.origin_parts.empty():
                    done_index = self.data_lines
                else:
                    done_index = self.origin_parts.get()

            if all([self.data_lines > 0, done_index < self.data_lines]):
                start_index = done_index
                is_last_part = (self.data_lines - done_index) <= self.max_sample
                if is_last_part:
                    done_index = self.data_lines
                else:
                    done_index += self.max_sample
                self.origin_parts.put(done_index)
                self.data_lock.release()
                data_part = self.data.iloc[start_index: done_index]
                self.hashing_parts.put(self.__hash_encode(data_part, cols))
                if done_index < self.data_lines:
                    self.__require_data(cols, process_index)
            else:
                self.data_lock.release()
        else:
            self.data_lock.release()

    def fit_transform(self, data=None):
        """Fit encoder according to X and y.

                Parameters
                ----------

                data : array-like, shape = [n_samples, n_features]
                    Training vectors, where n_samples is the number of samples
                    and n_features is the number of features.

                Returns
                -------

                data : data after hashing encoding, concat, reset_index

                """
        if data is None or self.cols is None:
            raise AttributeError("None data or feature input")

        self.data = data
        self.data_lines = len(data)
        if self.max_sample == 0 and self.max_process == 1:
            self.max_sample = self.data_lines

            self.__require_data(self.cols, 1)
        else:
            if self.max_sample == 0:
                self.max_sample = int(self.data_lines / self.max_process)
            for thread_index in range(self.max_process):
                process = multiprocessing.Process(target=self.__require_data,
                                                  args=(self.cols, thread_index + 1,))
                process.daemon = True
                self.n_process.append(process)
            for process in self.n_process:
                process.start()
            for process in self.n_process:
                process.join()
        if self.max_sample == 0 or self.max_sample == self.data_lines:
            return None if self.hashing_parts.empty() else self.hashing_parts.get()
        else:
            concat_data = []
            while not self.hashing_parts.empty():
                concat_data.append(self.hashing_parts.get())
            data = pd.concat(concat_data, ignore_index=True).reset_index(drop=True)
            return data
