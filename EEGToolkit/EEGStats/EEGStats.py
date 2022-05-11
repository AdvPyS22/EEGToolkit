"""
Computes time-point-wise the mean and SEM values of EEG data stored in numpy ndarrays.
"""

import numpy as np
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW
import matplotlib.pyplot as plt

# for custom legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ----------------------------------------------------------------
#                        Data Statistics 
# ----------------------------------------------------------------


def mean(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    Computes the element-wise mean of an `m x n numpy ndarray` 
    with `m` repeated datasets and `n` entries per set along axis 1
    (between repeated sets).

    Parameters
    ----------
    extracted_EEGData : np.ndarray
        An m x n numpy ndarray containing EEG data

    Returns
    -------
    means : np.ndarray
        An 1 x n numpy ndarray containing the element-wise means between all m repeated sets.
    """
    means = DescrStatsW( extracted_EEGData ).mean
    return means

def sem(extracted_EEGData:np.ndarray) -> np.ndarray:
    """
    Compute the standard error of the mean (SEM)
    of an `m x n numpy ndarray` with `m` repeated 
    datasets and `n` entries per set along axis 1
    (between repeated sets).

    Parameters
    ----------
    extracted_EEGData : np.ndarray
        An m x n numpy ndarray containing EEG data

    Returns
    -------
    sems : np.ndarray
        An 1 x n numpy ndarray containing the SEM between all m repeated sets.

    """
    sems = DescrStatsW( extracted_EEGData ) 
    length = len( extracted_EEGData )
    sems = sems.std / np.sqrt( length )
    return sems



def compare_signals(
                        extracted_EEGData_1 : np.ndarray, 
                        extracted_EEGData_2 : np.ndarray
                    ) -> np.ndarray:
    """
    Compares two EEG signal datasets stored as identical `m x n` np ndarrays
    with `m` repeated datasets and `n` entries per set, element-wise (along axis 1).
    It performs T-Tests to check for siginificant differences between the two datasets
    and returns the corresoponding p-values as a new array.

    Parameters
    ----------
    extracted_EEGData_1 : np.ndarray
        The first EEG signal dataset.
        
    extracted_EEGData_2 : np.ndarray
        The second EEG signal dataset.
    
    Returns
    -----
    pvalues : np.ndarray
        A 1 x n numpy ndarray of p-values for the 
        difference between Signals 1 and 2 at each position.
    """

    # compare the signal means
    diff = CompareMeans.from_data(
                                    extracted_EEGData_1, 
                                    extracted_EEGData_2
                                )

    # perform t-test comparison
    diff = diff.ttest_ind()

    # get the p-values 
    pvalues = diff[1]
    return pvalues

# ----------------------------------------------------------------
#                        Data Visualisation 
# ----------------------------------------------------------------

# set up some default style settings 

# colorscheme
signal_color = "gray"
signal1_color = "xkcd:dark royal blue"
signal2_color = "xkcd:watermelon"
signif_shade_color = "xkcd:yellow tan"

# opacities
signal_alpha = 0.8
sem_alpha = 0.2
signif_shade_alpha = 0.5

def plot_signal(
                    extracted_EEGData : np.ndarray,
                    sampling_frequency : float,
                    start_sec : float,
                    stop_sec : float,
                    x_scale = 10**3,
                    y_scale = 10**3,
                    baseline : np.ndarray = None,
                    significance_level:float = 0.05,
                    make_legend = False,
                    ax = None ) -> None:
    """
    Visualises a single EGG signal dataset from an `m x n numpy ndarray`
    with `m` repeated datasets and `n` entries per set. It generates a solid
    line for the time-point-wise mean and a shaded area of the corresponding SEM. 

    Parameters
    ----------

    extracted_EEGData : np.ndarray
        An m x n numpy ndarray of repeated EEG data.
    
    sampling_frequency : float
        The frequency in which the EEG data was recorded in `Hertz`.
    
    start_sec : float  
        The initial / low bound of the time-scale in seconds.
    
    stop_sec : float
        The final / upper bound of the time-scale in seconds.

    x_scale : float
        A scaling factor to adjust the data's x-value range. 
        E.g. `x_scale = 1000` to adjust the time-scale to milliseconds.
    
    y_scale : float
        A scaling factor for the data's y-value range.
        E.g. `y_scale = 1000` to adjust the signal-scale to millivolts.
    
    baseline : np.ndarray
        An array of baseline data to compare the signal to using a pair-wise t-test.
        This will introduce shaded boxes for regions that show significant differences
        between the signal and the baseline.
    
    significance_level : float
        The significance threshold for which to accept a signal difference
        as significant. Default is `0.05`.

    make_legend : bool
        Generates a legend for the plot.

    Returns
    -----
    pvalues : np.ndarray    
        The p-value from pair-wise T-Test comparison between the 
        signal and the baseline (if provided) at each signal position (timepoint).
    """

    # generate a new figure if no ax is specified
    if ax is None:
        fig, ax = plt.subplots()
    else: 
        fig = None

    pvalues = _plot_(ax,
                    extracted_EEGData,
                    sampling_frequency,
                    start_sec,
                    stop_sec,
                    x_scale,
                    y_scale, 
                    baseline,
                    make_legend,
                    significance_level = significance_level 
                )

    # we show the figure only if no ax was provided
    # and thus no "bigger" figure is assembled elsewhere...
    if fig is not None:
        plt.show()
    
    return pvalues

def difference_plot(extracted_EEGData_1:np.ndarray,
                    extracted_EEGData_2:np.ndarray,
                    sampling_frequency:float,
                    start_sec:float,
                    stop_sec:float,
                    significance_level:float = 0.05,
                    x_scale=10**3,
                    y_scale=10**3,
                    make_legend = False, 
                    ax = None ) -> None:

    """
    Visualises the difference between two EEG signals and tests
    time-point-wise the difference using T-Tests. Individual EEG
    Signals are shown as mean-line with shaded SEM, and time-points
    with significant differences are highlighted with overall-shading.

    Parameters
    ----------
    extracted_EEGData_1 : np.ndarray
        The first EEG signal dataset.
        
    extracted_EEGData_2 : np.ndarray
        The second EEG signal dataset.
    
    sampling_frequency : float
        The frequency in which the EEG data was recorded in `Hertz`.
    
    significance_level : float
        The significance threshold for which to accept a signal difference
        as significant. Default is `0.05`.

    start_sec : float  
        The initial / low bound of the time-scale in seconds.

    stop_sec : float
        The final / upper bound of the time-scale in seconds.

    x_scale : float
        A scaling factor to adjust the data's x-value range. 
        E.g. `x_scale = 1000` to adjust the time-scale to milliseconds.
    
    y_scale : float
        A scaling factor for the data's y-value range.
        E.g. `y_scale = 1000` to adjust the signal-scale to millivolts.
    
    make_legend : bool
        Generates a legend for the plot.
    
    Returns
    -----
    pvalues : np.ndarray    
        The p-value from pair-wise T-Test comparison between the 
        two signals at each signal position (timepoint).
    """

    # generate a new figure if no ax is given
    if ax is None:
        fig, ax = plt.subplots()
    else: 
        fig = None

    pvalues = _difference_plot_(ax,
                        extracted_EEGData_1,
                        extracted_EEGData_2,
                        significance_level,
                        sampling_frequency,
                        start_sec,
                        stop_sec,
                        x_scale,
                        y_scale,
                        make_legend 
                        )

    # we show the figure only if no ax was provided
    # and thus no "bigger" figure is assembled elsewhere...
    if fig is not None: 
        plt.show()

    return pvalues

def _plot_(
            ax, 
            extracted_EEGData:np.ndarray,
            sampling_frequency:float,
            start_sec:float,
            stop_sec:float,
            x_scale=10**3,
            y_scale=10**3, 
            baseline : np.ndarray = None, 
            make_legend = False,
            **kwargs ) -> None:
    """
    Generates a mean signal line with shaded 
    SEM area from an m x n numpy ndarray.

    Note
    -----
    This is the core for `plot_signal` 
    """
    # compute mean and SEM for the signal
    _mean = mean(extracted_EEGData)
    _sem = sem(extracted_EEGData)

    # now generate scaled xvalues
    x_values = _scale_xboundries(sampling_frequency, start_sec, stop_sec, x_scale)

    # now plot the signal's scaled mean line
    signal = _mean * y_scale 
    ax.plot(x_values, signal, color = signal_color )

    # and add the shading for SEM 
    lower = (_mean - _sem) 
    upper = (_mean + _sem)
    lower *= y_scale
    upper *= y_scale

    ax.fill_between(
                        x_values,
                        lower,
                        upper,
                        color = signal_color,
                        edgecolor = None,
                        linewidth = 0,
                        alpha = sem_alpha
                        
                )

    # if we have baseline data, we compare also to the baseline
    # and shade siginificantly different regions...
    pvalues = None
    if baseline is not None:

        pvalues = compare_signals( extracted_EEGData, baseline )

        # generate the y-value boundries for the plot
        # and scale the y-values to some user-defined range
        # plus add a little padding
        max_y = np.max( _mean )
        min_y = np.min( _mean )
        yvalues = _scale_yboundries( y_scale, max_y, min_y, pad = 1.2 )

        # add the shaded fillings for significantly different
        # timepoint areas
        signif_level = kwargs.pop( "significance_level", 0.05 )
        _shade_singificant_regions(
                                    ax, 
                                    significance_level = signif_level, 
                                    pvalues = pvalues, 
                                    xvalues = x_values,
                                    ylim = yvalues, 
                                )



    # and add a line for the start of the signal
    ax.axvline( x=0, linewidth = 2, color = "black" )

    # and some axes formatting...
    ax.set_title("Average EEG Signal (Shaded area SEM)")
    ax.set_ylabel("Signal\namplitude")
    ax.set_xlabel("Time relative to event")
    
    if make_legend:
        handles = [
                    Line2D(     [0], [0], 
                                color = signal_color, 
                                label = "Mean Signal" 
                        ),
                    Patch( 
                                facecolor = signal_color, 
                                edgecolor = None, linewidth = 0,
                                alpha = sem_alpha,
                                label = "SEM"
                        )
                ]
        if baseline is not None: 
            handles.append(
                            Patch( 
                                    facecolor = signif_shade_color, 
                                    edgecolor = None, linewidth = 0,
                                    alpha = signif_shade_alpha,
                                    label = f"pvalue < {signif_level}"
                            )
                        )
        # now add a custom legend
        ax.legend( handles = handles,
                            bbox_to_anchor = (1, -0.5),
                            frameon = False
                        )
    return pvalues

def _difference_plot_(ax,
                      extracted_EEGData_1:np.ndarray,
                      extracted_EEGData_2:np.ndarray,
                      significance_level:float,
                      sampling_frequency:float,
                      start_sec:float,
                      stop_sec:float,
                      x_scale=10**3,
                      y_scale=10**3,
                      make_legend = False ) -> None: 

    """
    Generates a plot of two EEG signals and their time-point-wise 
    difference using a T-Test. 

    Note
    -----
    This is the core of difference_plot
    """

    # compare the two EEG signals and generate a p-value 
    # from t-test position-wise comparisons...
    pvalues = compare_signals(extracted_EEGData_1, extracted_EEGData_2)

    # generate the mean lines for both signals
    mean_1 = mean(extracted_EEGData_1)
    mean_2 = mean(extracted_EEGData_2)

    # generate the y-value boundries for the plot
    max_y = np.max( [np.max(mean_1), np.max(mean_2)] )
    min_y = np.min( [np.min(mean_1), np.min(mean_2)] )

    # and scale the y-values to some user-defined range
    # plus add a little padding
    yvalues = _scale_yboundries( y_scale, max_y, min_y, pad = 1.2 )

    # generate correspondingly scaled x-values
    x_values = _scale_xboundries(sampling_frequency, start_sec, stop_sec, x_scale)

    # add the shaded fillings for significantly different
    # timepoint areas
    _shade_singificant_regions(
                                ax, 
                                significance_level = significance_level, 
                                pvalues = pvalues, 
                                xvalues = x_values,
                                ylim = yvalues, 
                            )

    # plot scaled signals 1 and 2
    signal_1 = mean_1 * y_scale
    signal_2 = mean_2 * y_scale

    ax.plot(x_values, signal_1, color = signal1_color, alpha = signal_alpha )
    ax.plot(x_values, signal_2, color = signal2_color, alpha = signal_alpha )

    # plot a vertical line at the signal start
    ax.axvline(x = 0, color = "black", linewidth = 2 )

    # and some axes formatting...
    ax.set_title("Average EEG Signal (Shaded area significant regions)")
    ax.set_ylabel("Signal\namplitude")
    ax.set_xlabel("Time relative to event")

    if make_legend: 
        # now add a custom legend
        ax.legend( handles = [
                                Line2D(     [0], [0], 
                                            color = signal1_color, 
                                            label = "Mean horizontal Signal" 
                                    ),
                                Line2D(     [0], [0], 
                                            color = signal2_color, 
                                            label = "Mean vertical Signal" 
                                    ),
                                Patch( 
                                            facecolor = signif_shade_color, 
                                            edgecolor = None, linewidth = 0,
                                            alpha = signif_shade_alpha,
                                            label = f"pvalue < {significance_level}"
                                    )

                            ],
                            bbox_to_anchor = (1, -0.5),
                            frameon = False
                        )
    return pvalues

def _shade_singificant_regions(ax, significance_level : float , pvalues : np.ndarray , xvalues : np.ndarray , ylim : tuple ):
    """
    Shades the background of regions (x-value ranges) where corresponding p-values
    are below a given significance level. 

    Parameters
    ----------
    ax : plt.axes.Axes 
        An Axes object to plot to.
    significance_level : float
        The significance level to use.
    pvalues : np.ndarray
        An array of pvalues for each x-value.
    xvalues : np.ndarray
        An array of x-values.
    ylim : tuple
        A tuple of minimal and maximal y-values to shade.
    """
    ax.fill_between(
                     xvalues,
                     ylim[0],
                     ylim[1],
                     # now fill areas where the pvalues are 
                     # below our significance_level
                     where = pvalues < significance_level, 
                     facecolor = signif_shade_color,
                     alpha = signif_shade_alpha
                )

def _scale_xboundries(sampling_frequency, start_sec, stop_sec, x_scale):
    """
    Scales minimal and maximal x values from starting and final timestep values
    given some scale and sampling frequency.
    """
    x_values = np.arange(
                            int( start_sec * sampling_frequency ),
                            int( stop_sec * sampling_frequency ), 
                            dtype=int
                        )
    x_values = x_values / sampling_frequency * x_scale
    return x_values

def _scale_yboundries(y_scale, max_y, min_y, pad = 1 ):
    """
    Scales minimal and maximal y values from some data
    by a given scale and adds additional padding as a scalar factor
    (pad = 1 means no padding).
    """
    max_y *= y_scale * pad
    min_y *= y_scale * pad
    return min_y, max_y