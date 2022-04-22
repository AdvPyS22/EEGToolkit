import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt

"""
"""
def mean(extracted_EEGData:np.ndarray) -> np.ndarray:
    return np.apply_along_axis(np.mean, 0, extracted_EEGData)

"""
"""
def sem(extracted_EEGData:np.ndarray) -> np.ndarray:
    return np.apply_along_axis(sp.sem, 0, extracted_EEGData)

"""
"""
def plot(extracted_EEGData:np.ndarray, \
         sampling_frequency:float, \
         start_sec:float, \
         stop_sec:float, \
         x_scale=10**3, \
         y_scale=10**3) -> None:
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