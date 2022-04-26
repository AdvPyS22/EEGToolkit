"""
TODO
"""

import numpy as np
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
import matplotlib.pyplot as plt


def mean(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    TODO
    """
    return DescrStatsW(extracted_EEGData).mean

def sem(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    """
    return DescrStatsW(extracted_EEGData).std / (len(extracted_EEGData)**(1/2))

def _plot_(ax, 
         extracted_EEGData:np.ndarray,
         sampling_frequency:float,
         start_sec:float,
         stop_sec:float,
         x_scale=10**3,
         y_scale=10**3) -> None:
    
    _mean = mean(extracted_EEGData)
    _sem = sem(extracted_EEGData)

    x_scale = np.array(range(int(start_sec*sampling_frequency),
                             int(stop_sec*sampling_frequency)))/sampling_frequency*x_scale

    ax.plot(x_scale, _mean*y_scale)
    ax.fill_between(x_scale,
                     (_mean - _sem)*y_scale,
                     (_mean + _sem)*y_scale,
                     alpha=0.2)
    ax.axvline(x=0)
    ax.set_title('Average EEG Signal (Shaded area SEM)')
    ax.set_ylabel("Signal amplitude")
    ax.set_xlabel("Time relative to event")

def plot(extracted_EEGData:np.ndarray,
         sampling_frequency:float,
         start_sec:float,
         stop_sec:float,
         x_scale=10**3,
         y_scale=10**3) -> None:
    """
    TODO
    """
    fig, ax = plt.subplots()
    
    _plot_(ax,
           extracted_EEGData,
           sampling_frequency,
           start_sec,
           stop_sec,
           x_scale,
           y_scale)
    
    plt.show()

def _difference_plot_(ax,
                      extracted_EEGData_1:np.ndarray,
                      extracted_EEGData_2:np.ndarray,
                      significance_level:float,
                      sampling_frequency:float,
                      start_sec:float,
                      stop_sec:float,
                      x_scale=10**3,
                      y_scale=10**3) -> None: 

    """
    TODO
    """

    diff = CompareMeans.from_data(extracted_EEGData_1, extracted_EEGData_2)
    pvalue = diff.ttest_ind()[1]

    mean_1 = mean(extracted_EEGData_1)
    mean_2 = mean(extracted_EEGData_2)

    max_y = np.max([np.max(mean_1), np.max(mean_2)])*y_scale
    min_y = np.min([np.min(mean_1), np.min(mean_2)])*y_scale

    x_scale = np.arange(int(start_sec*sampling_frequency),
                        int(stop_sec*sampling_frequency), dtype=int)/sampling_frequency*x_scale

    ax.plot(x_scale, mean_1*y_scale)
    ax.plot(x_scale, mean_2*y_scale)
    ax.axvline(x=0)
    ax.fill_between(x_scale,
                     min_y,
                     max_y,
                     where=pvalue < significance_level,
                     facecolor='red',
                     alpha=0.2)
    ax.set_title('Average EEG Signal (Shaded area significant regions)')
    ax.set_ylabel("Signal amplitude")
    ax.set_xlabel("Time relative to event")

def difference_plot(extracted_EEGData_1:np.ndarray,
                    extracted_EEGData_2:np.ndarray,
                    significance_level:float,
                    sampling_frequency:float,
                    start_sec:float,
                    stop_sec:float,
                    x_scale=10**3,
                    y_scale=10**3) -> None:
    """
    TODO
    """
    
    fig, ax = plt.subplots()
    _difference_plot_(ax,
                      extracted_EEGData_1,
                      extracted_EEGData_2,
                      significance_level,
                      sampling_frequency,
                      start_sec,
                      stop_sec,
                      x_scale,
                      y_scale)
    plt.show()