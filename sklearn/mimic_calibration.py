import numpy as np
from .base import BaseEstimator, RegressorMixin
from .utils import indexable, column_or_1d
import sys


class _MimicCalibration(BaseEstimator, RegressorMixin):
    """ mimic calibration
    reference:
    https://www.youtube.com/watch?v=Cg--SC76I1I
    http://tech.magnetic.com/2015/06/click-prediction-with-vowpal-wabbit.html
    """
    def __init__(self, threshold_pos=5, boundary_choice=2, record_history=False):
        self.threshold_pos = threshold_pos
        self.boundary_choice = boundary_choice
        self.record_history = record_history
        self.history_record_table = []

    def get_bin_boundary(self, current_binning, boundary_choice):
        """
        current_binning:
        [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, PosRate_temp]]

        boundary_choice:
        0: choose socre_min, ie left boundary of bin
        1: choose socre_max, ie right boundary of bin
        2: choose socre_mean, ie mean score of bin
        """
        num_rows = len(current_binning)
        boundary_table_temp = []
        k = None
        if (boundary_choice == 0):
            k = 1
        elif (boundary_choice == 1):
            k = 2
        elif (boundary_choice == 2):
            k = 3
        else:
            raise Exception("Un-identified boundary choice: {x}".format(x=boundary_choice))
        for i in range(num_rows):
            boundary_table_temp += [current_binning[i][k]]
        return boundary_table_temp

    def construct_initial_bin(self,
                              sorted_score,
                              sorted_target,
                              threshold_pos):
        """make each bin having the number of positives equal to threshold_pos.
        the default number of threshold_pos = 5.

        Parameters
        ----------
        sorted_score: the sorted probability from the model, ie pre-calibrated score.
        sorted target: the target in the order of increasing score. the number of target = 2.
        threshold_pos: number of positive in each bin, default=5

        Returns
        ----------
        it return the info of each bin and every bin has nPos = threshold_pos.
        bl_index: corresponding left index in the sorted score/target.
        score_min, score_max, score_mean: score/probability info in each bin.
        nPos: the number of positive in the bin which should be threshold_pos except the last bin.
        nPosRate: the positive rate in the bin.
        bin_info: [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, nPosRate_temp]]

        """
        bin_right_index_array = []
        last_index = len(sorted_target)-1
        count = 0
        # make each bin having number of positive = threshold positive.
        # bin_right_index_array: its element is right-boundary index of each bin.
        for i in range(len(sorted_target)):
            y = sorted_target[i]
            if y > 0:
                count += 1
            if (count == threshold_pos):
                bin_right_index_array += [i]
                count = 0
        if (len(sorted_target)-1 not in bin_right_index_array):
            bin_right_index_array += [last_index]

        # bl_index: left boundary index of each bin.
        bl_index = 0
        bin_info = []
        total_number_pos = 0
        for br_index in bin_right_index_array:
            # score stats
            score_temp = sorted_score[bl_index: br_index + 1]
            score_min = min(score_temp)
            score_max = max(score_temp)
            score_mean = np.mean(score_temp)
            # target
            target_row = sorted_target[bl_index: br_index + 1]
            nPos_temp = np.sum(target_row)
            if (br_index != last_index):
                assert (nPos_temp == threshold_pos), "The sum of positive must be equal to threshold pos except the last index."
            total_number_per_bin = len(target_row)
            nPosRate_temp = 1.0*nPos_temp/total_number_per_bin
            bin_info += [[bl_index, score_min, score_max, score_mean, nPos_temp, total_number_per_bin, nPosRate_temp]]
            total_number_pos += nPos_temp
            bl_index = br_index + 1
        return bin_info, total_number_pos

    def merge_bins(self, binning_input, increasing_flag):
        # binning_input
        # [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, PosRate_temp]]
        nbins = len(binning_input)
        result = []
        for i in range(1, nbins):
            # current_bin: latest new bin in the result
            if (i == 1):
                result += [binning_input[0]]
            current_bin = result[-1]
            current_bin_PosRate = current_bin[-1]
            next_bin = binning_input[i]
            next_bin_PosRate = next_bin[-1]
            if(current_bin_PosRate > next_bin_PosRate):
                increasing_flag = False
                # merge two bins:
                # [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, PosRate_temp]]
                new_bin_index_temp = min(current_bin[0], next_bin[0])
                new_score_min_temp = min(current_bin[1], next_bin[1])
                new_score_max_temp = max(current_bin[2], next_bin[2])
                new_score_mean_temp = (current_bin[3] + next_bin[3])/2.0
                new_pos_temp = current_bin[4] + next_bin[4]
                new_total_temp = current_bin[5] + next_bin[5]
                new_PosRate_temp = 1.0*new_pos_temp/new_total_temp
                # update the latest bin info in the latest result
                result[-1] = [new_bin_index_temp, new_score_min_temp, new_score_max_temp,
                              new_score_mean_temp, new_pos_temp, new_total_temp, new_PosRate_temp]
            else:
                result += [next_bin]
        return result, increasing_flag

    def run_merge_function(self, current_binning, record_history=False):
        # current_binning
        # [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, PosRate_temp]]
        self.history_record_table = []
        if (record_history):
            self.history_record_table += [current_binning]

        keep_merge = True
        while(keep_merge):
            new_bin_temp, increasing_flag = self.merge_bins(current_binning, True)
            if (record_history):
                self.history_record_table += [new_bin_temp]

            # update the current_binning
            current_binning = new_bin_temp
            # if it increasing monotonically, we stop merge
            keep_merge = not increasing_flag
        # if (record_history):
        #     return self.history_record_table
        return [new_bin_temp]

    def _mimic_calibration(self,
                           y_score,
                           y_target,
                           number_positive_within_bin=5):
        """ y_score: the prediction from the model, probability
        y_target: 0 or 1
        """
        assert ((y_score.min() >= 0) & (y_score.max() <= 1.0)), "y_score is a probability which is between 0 and 1."
        assert (len(np.unique(y_score)) > 2), "y_score should be at least 3 different probability."
        assert np.array_equal(np.unique(y_target), np.array([0, 1])), "y_traget must be 0 and 1."
        y_score = column_or_1d(y_score)
        y_target = column_or_1d(y_target)
        # sort y_score
        sorted_index = y_score.argsort()
        y_score = y_score[sorted_index]
        y_target = y_target[sorted_index]
        threshold_pos = number_positive_within_bin
        initial_binning, total_number_pos = self.construct_initial_bin(y_score, y_target, threshold_pos)
        final_binning = self.run_merge_function(initial_binning, self.record_history)
        latest_bin_temp = final_binning[-1]
        boundary_table = self.get_bin_boundary(latest_bin_temp, self.boundary_choice)
        return boundary_table, latest_bin_temp

    def fit(self, X, y, sample_weight=None):
        X = column_or_1d(X)
        y = column_or_1d(y)
        X, y = indexable(X, y)
        self.boundary_table, self.calibrated_model = self._mimic_calibration(X, y, self.threshold_pos)
        return self

    def predict_per_element(self, x, num_boundary):
        # linear interpolation
        which_bin = np.digitize([x], self.boundary_table, right=True)[0]
        if ((which_bin == 0) or (which_bin == len(self.boundary_table)-1)):
            y = self.calibrated_model[which_bin][6]
        else:
            if (which_bin >= num_boundary):
                # outside bin boundary
                delta_y = 1.0 - self.calibrated_model[num_boundary-1][6]
                delta_x = x - self.boundary_table[num_boundary-1]
                if (delta_x != 0):
                    y = self.calibrated_model[num_boundary-1][6] + \
                        (1.0*delta_y/delta_x) * (x - self.boundary_table[num_boundary-1])
                else:
                    y = self.calibrated_model[num_boundary-1][6]
            else:
                delta_y = self.calibrated_model[which_bin][6] - self.calibrated_model[which_bin-1][6]
                delta_x = self.boundary_table[which_bin] - self.boundary_table[which_bin-1]
                y = self.calibrated_model[which_bin-1][6] + \
                    (1.0*delta_y/delta_x) * (x - self.boundary_table[which_bin-1])
        return y

    def predict(self, pre_calib_prob):
        """x_array: a raw probability, ie pre-calibrated probability.
        calibrated_model:
        [[bl_index, score_min, score_max, score_mean, nPos_temp, total_temp, PosRate_temp]]
        Use calibrated_model to predict the calibrated probability.
        Perform linear interpolation to find the calibrated probability.
        """
        pre_calib_prob = column_or_1d(pre_calib_prob)
        if(self.calibrated_model is None):
            sys.exit("Please calibrate model first by calling fit function.")
        else:
            calib_prob = []
            num_boundary = len(self.boundary_table)
            for x in pre_calib_prob:
                y = self.predict_per_element(x, num_boundary)
                calib_prob += [y]
            calib_prob = np.array(calib_prob)
        return calib_prob

    def get_one_history(self, one_history):
        score_array = []
        nP_array = []
        for row in one_history:
            # the mean of score at each bin
            score = row[3]
            # the nPos rate at each bin
            nP = row[6]
            score_array += [score]
            nP_array += [nP]
        return score_array, nP_array

    def plot_merge_result(self):
        import matplotlib.pyplot as plt
        # self.calibrated_model
        data = None
        if (self.record_history):
            data = self.history_record_table
        else:
            data = self.calibrated_model

        number_of_history = len(data)
        print("plot history size: {x}".format(x=number_of_history))
        for i in range(number_of_history):
            one_history = data[i]
            score_array, nPosRate_array = self.get_one_history(one_history)
            plt.plot(score_array, nPosRate_array, label=str(i))
        plt.xlabel("score")
        plt.ylabel("nPos Rate")
        plt.legend()
        plt.show()
