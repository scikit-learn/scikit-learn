# Authors: Nishant Baheti <nishantbaheti.it19@gmail.com>

from typing import Optional, Union
import numpy as np


class PSDDenoiser:
    """PSD (Power Spectral Density) Based Denoiser

    This method takes the FFT transform of the signal to calculate PSD
    based on the PSD results and cutoff threshold the signal is filtered and
    a FFT inverse is applied to regenerate denoised signal.

    Parameters
    ----------
    threshold : Optional[Union[int, float, str]], optional
        threshold to create cutoff mask, but any threshold can be applied,
        if it is precalculated by any method chosen by the process, by default auto-mean
        { auto-mean, auto-max }

    Examples
    --------
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.preprocessing import PSDDenoiser
        >>> rng = np.random.default_rng()
        >>> fs = 10e3
        >>> N = 100
        >>> amp = 2 * np.sqrt(2)
        >>> freq = 1234.0
        >>> noise_power = 0.001 * fs / 2
        >>> time = np.arange(N) / fs
        >>> X = amp * np.sin(2 * np.pi * freq * time)
        >>> X += rng.normal(scale=np.sqrt(noise_power), size=time.shape)

        >>> denoiser = PSDDenoiser()
        >>> cleaned_signal = denoiser.transform(X)
        >>> plt.plot(X, label="noisy")
        >>> plt.plot(cleaned_signal, label="cleaned")
        >>> plt.title(f"Threshold : {denoiser.threshold}")
        >>> plt.legend(loc="best")
        >>> plt.show()

        >>> denoiser = PSDDenoiser(10)
        >>> cleaned_signal = denoiser.transform(X)
        >>> plt.plot(X, label='noisy')
        >>> plt.plot(cleaned_signal, label='cleaned')
        >>> plt.title(f"Threshold : {denoiser.threshold}")
        >>> plt.legend(loc='best')
        >>> plt.show()
    """

    def __init__(
        self, threshold: Optional[Union[int, float, str]] = "auto-mean"
    ) -> None:
        self._threshold = self.__init_threshold(threshold)
        self._f_hat = None
        self._pxx = None
        self._cutoff_mask = None
        self._filtered_f_hat = None
        self._denoised_X = None

    def __init_threshold(self, threshold):
        if isinstance(threshold, (str,)):
            threshold = threshold.lower()
            assert threshold in (
                "auto-mean",
                "auto-max",
            ), "available auto threshold methods { auto-mean, auto-max, auto-min }"
        return threshold

    def _reshape_X(self, X: np.ndarray):
        if len(X.shape) == 1:
            X = X.reshape(-1, 1)
        return X

    def psd(self, f_hat: np.ndarray, tau: int) -> np.ndarray:
        """Power Spectral Density

        Parameters
        ----------
        f_hat : np.ndarray
            Signal in Frequency Domain
        tau : int
            Interval

        Returns
        -------
        np.ndarray
            Power spectrum
        """
        return ((f_hat * np.conjugate(f_hat)) / tau).real

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Apply PSD

        Parameters
        ----------
        X : np.ndarray
            Input matrix, signal in IOT Terms.

        Returns
        -------
        np.ndarray
            Denoised Signal
        """
        X = self._reshape_X(X)
        tau = len(X)
        self._f_hat = np.fft.fft(X, tau, axis=0)
        self._pxx = self.psd(self._f_hat, tau)

        if isinstance(self._threshold, str):
            agg_func = getattr(np, self._threshold.split("-")[1])
            self._threshold = agg_func(self._pxx[np.int16(np.floor(tau / 2)) :])
        self._cutoff_mask = self._pxx > self._threshold
        self._filtered_f_hat = self._f_hat * self._cutoff_mask
        self._denoised_X = np.fft.ifft(self._filtered_f_hat, axis=0).real
        return self._denoised_X

    @property
    def threshold(self) -> Union[str, float, int]:
        """Threshold calculated by the process

        In Power spectrum after the half length it takes the aggregation of the values
        and use that as a threshold to cutoff frequencies that are insignificant.

        Returns
        -------
        Union[str, float, int]
            cutoff threshold value
        """
        return self._threshold

    @property
    def f_hat(self) -> Optional[np.ndarray]:
        """FFT of input signal"""
        return self._f_hat

    @property
    def filtered_f_hat(self) -> Optional[np.ndarray]:
        """filtered FFT of input signal"""
        return self._filtered_f_hat
