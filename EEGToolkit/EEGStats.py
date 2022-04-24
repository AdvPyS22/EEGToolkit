from itertools import groupby
from operator import itemgetter
import numpy as np
import scipy.stats as sp
from statsmodels.stats.weightstats import CompareMeans
import matplotlib.pyplot as plt


def mean(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    """
    return np.apply_along_axis(np.mean, 0, extracted_EEGData)

def sem(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    """
    return np.apply_along_axis(sp.sem, 0, extracted_EEGData)

def plot(extracted_EEGData:np.ndarray,
         sampling_frequency:float,
         start_sec:float,
         stop_sec:float,
         x_scale=10**3,
         y_scale=10**3) -> None:
    """
    """
    
    _mean = mean(extracted_EEGData)
    _sem = sem(extracted_EEGData)
    
    x_scale = np.array(range(int(start_sec*sampling_frequency),
                             int(stop_sec*sampling_frequency)))/sampling_frequency*x_scale
    
    plt.plot(x_scale, _mean*y_scale)
    plt.fill_between(x_scale, (_mean - _sem)*y_scale, (_mean + _sem)*y_scale, alpha=0.2)
    plt.axvline(x=0)
    plt.title('Average EEG Signal (Shaded area SEM)')
    plt.ylabel("Signal amplitude")
    plt.xlabel("Time relative to event")
    plt.show()

def difference_plot(extracted_EEGData_1:np.ndarray,
                    extracted_EEGData_2:np.ndarray,
                    significance_level:float,
                    sampling_frequency:float,
                    start_sec:float,
                    stop_sec:float,
                    x_scale=10**3,
                    y_scale=10**3) -> None:
    """
    """

    diff = CompareMeans.from_data(extracted_EEGData_1, extracted_EEGData_2)
    pvalue = diff.ttest_ind()[1]

    mean_1 = mean(extracted_EEGData_1)
    mean_2 = mean(extracted_EEGData_2)

    max_y = np.max([np.max(mean_1), np.max(mean_2)])*y_scale
    min_y = np.min([np.min(mean_1), np.min(mean_2)])*y_scale

    x_scale = np.array(range(int(start_sec*sampling_frequency),
                             int(stop_sec*sampling_frequency)))/sampling_frequency*x_scale

    plt.plot(x_scale, mean_1*y_scale)
    plt.plot(x_scale, mean_2*y_scale)
    plt.fill_between(x_scale, min_y, max_y, where=pvalue < significance_level,
                facecolor='red', alpha=0.2)
    plt.title('Average EEG Signal (Shaded area significant regions)')
    plt.ylabel("Signal amplitude")
    plt.xlabel("Time relative to event")
    plt.show()