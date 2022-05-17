"""

This module provides a data class `EEGData` to work with EEG signal data for event-reaction-time delay experiments.
It works with two separate input datafiles, one storing the EEG signal itself as a 1D array, and one 
describing event metadata as a 2D array, describing both the timepoints and the type of event in two columns.

"""

import subprocess
import inspect
import os 
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from EEGToolkit.EEGStats import plot_signal, difference_plot

supported_filetypes = [ "npy", "tsv", "csv", "txt" ]
class EEGData():
    """
    Handles EEG data stored from two separate datafiles, one for the EEG signal
    and one for the corresponding events. 

    Input Datafiles
    --------------
    The EEG signal datafile must be a 1D array of values. 
    The events datafile must be a 2D array of timepoints at which an event occurs,
    as well as a categorical (numerically encoded) label of the kind of event that occured. Note, that the
    datafiles must **not** contain any headers!

    Supported file types are:
    - `npy`
    - `txt`     ( space-separated for events datafiles )
    - `tsv`
    - `csv`     (both `,` and `;` separated )   

    Parameters
    ----------
    signal_path : str
        A filepath to a valid datafile containing EEG signal data.
    
    event_path : str
        A filepath to a valid datafile containing corresponding 
        event information for the signal datafile.

    sampling_frequency : float
        The frequency at which the EEG was recorded in `Hertz`.

    """

    def __init__(self,
                 signal_path:str,
                 event_path:str,
                sampling_frequency:float) -> None:

        # store the used datafiles...
        self._signal_src = signal_path
        self._event_src = event_path

        # first check for valid input data
        self._check_sanity(signal_path, event_path, sampling_frequency)

        # read the datafiles...
        self.read( signal_path = signal_path, event_path = event_path )    

        # now setup the frames for the events 
        self._n_frames = len(self.signal)

        # this will setup self._events which is a
        # dictionary of event identifiers : number of repeated measurements
        self._set_n_events()

        # setup a _data argument for the 
        # extracted event datasets
        self._data = None 

        # store the frequency
        self.sampling_frequency = sampling_frequency
        
        # setup default parameters for the window around
        # which each event signal should be extracted
        # in seconds
        self._start_sec = -0.5
        self._stop_sec = 1

        # and setup the extracted events in case 
        # only a subset are being extacted
        self._extracted_events = None

        # save a baseline for each event type which will be an 
        # np.ndarray to store the timepoints (or a subset thereof)
        # before the signal onset. The time points will be sampled
        # from the extracted timepoints...
        self._baseline = None

        # setup a dictionary to store the p-values of pair-wise comparison
        # between either two signals or a signal with it's baseline.
        # keys will be tuples of signal1, signal2, for baseline comparison
        # signal1 = signal2...
        self._pvalues = {}

    def read( self, signal_path : str = None , event_path : str = None ) -> None: 
        """
        Read the provided data files and stores the
        data into numpy ndarrays.

        Note
        ----
        This method is automatically called at initiation. However, new data
        can be loaded using this method manually.

        Input Datafiles
        --------------
        The EEG signal datafile must be a 1D array of values. 
        The events datafile must be a 2D array of timepoints at which an event occurs,
        as well as a categorical (numerically encoded) label of the kind of event that occured. Note, that the
        datafiles must **not** contain any headers!

        Supported file types are:
        - `npy`
        - `txt`     ( space-separated for events datafiles )
        - `tsv`
        - `csv`     (both `,` and `;` separated )   

        Parameters
        ----------
        signal_path : str
            A filepath to a valid datafile containing EEG signal data.
        
        event_path : str
            A filepath to a valid datafile containing corresponding 
            event information for the signal datafile.

        """
        
        # first read the signal data file
        if signal_path is not None:
            suffix = self._filesuffix( signal_path )
            if suffix == "npy": 
                signal = self._read_npy( signal_path )
            else: 
                signal = self._read_datafile( signal_path )

            # now save
            self.signal = signal

        # now the same for the events data file 
        if event_path is not None: 
            suffix = self._filesuffix( event_path )
            if suffix == "npy": 
                events = self._read_npy( event_path )
            else: 
                events = self._read_datafile( event_path )

            # now save
            self._events_data = events

    def extract(self,
                start_sec:float,
                stop_sec:float,
                event_type : ( int or tuple or list or np.ndarray ) = None, 
                **kwargs ) -> np.ndarray:
        """
        Extracts data for a specific (set of) event(s) from the loaded data. 
        And returns the data as numpy ndarrays (or a list thereof, in case of 
        multiple events).

        Note
        ----
        This method is automatically called by `EGGData.summary` so it is not necessary
        to manually extract data, unless only a specific subset of event types should
        be extracted!

        Parameters
        ----------
        start_sec : float
            The timepoint relative to the provided 
            events data in seconds at which to begin extraction. 
            E.g. `-0.05` would correspond to `0.05 seconds` before
            the actual onset of the recorded event.

        stop_sec : float
            The timepoint relative to the provided 
            events data in seconds at which to end extraction. 
            E.g. `0.75` would correspond to `0.75 seconds` after
            the actual onset of the recorded event.
        
        event_type : int or tuple or list or np.ndarray
            Either a single event type or an iterable of multiple event types.
            If `event_type = None` (default) data for **all** event types will be extracted!

        Returns
        -------
        extracted_data : np.ndarray or list
            The data of the provided events as an ndarray or a 
            list of ndarrays in the same order as the provided event-type labels.
        """

        # check if we should extract the data for all events
        if event_type is None:

            # get all events
            events_to_extract = self._events.keys()

            # extract each type from the loaded data
            data = [ 
                        self.extract(start_sec, stop_sec, etype) 
                        for etype in events_to_extract 
                ]

            self._data = data
            self._extracted_events = events_to_extract
            return data

        # check if there is a provided subset of events to extract
        if isinstance(event_type, (tuple, np.ndarray, list)):

            events_to_extract = event_type

            # extract provided type from the loaded data
            data = [ 
                        self.extract(start_sec, stop_sec, etype) 
                        for etype in events_to_extract 
                ]

            self._data = data
            self._extracted_events = events_to_extract
            return data

        # now the part for extracting only a 
        # single event type data
        data = self._extract_window(start_sec, stop_sec, event_type)
        self._data = data

        self._extracted_events = event_type

        # store start and stop sec values
        # for later use in summary()
        self._start_sec = start_sec
        self._stop_sec = stop_sec

        return data

    def baseline( self, size : int or float = None ):
        """
        Generates a baseline distribution for EEG Signals,
        using random sampling from pre-signal timepoints accross 
        replicates and events.

        Note
        ----
        This requires that events have already been extacted!

        Parameters
        ----------
        size : int or float
            The number of random samples to draw. If `None` are provided (default)
            the entire available pre-signal data is used. If an `integer` is provided
            then the final baseline data contains exactly the given number of datapoints.
            Alternatively, a `float` `0 < size <= 1` may be provided to specify a fraction
            of data to sample from. E.g. `size = 0.2` would incorporate 20% of the available
            datapoints into the baseline.
        
        Returns
        -------
        baseline : np.ndarray
            An np.ndarray of the randomly drawn samples.
        """
        start_sec, stop_sec = self._start_sec, self._stop_sec

        # first get the time period before t=0, beginning at the starting time...
        if isinstance( self._data, list ):
            random_samples = [ self._extract_window( start_sec, 0, e ) for e in self._extracted_events ]
        elif isinstance( self._data, np.ndarray ):
            random_samples = [ self._extract_window( start_sec, 0, self._extracted_events ) ]        
        elif self._data is None:
            raise Exception( f"No events data has been extracted yet! Make sure to run extract() before computing a baseline." )

        # collapse the repeats into a single dataset
        random_samples = [ np.reshape( i, i.size ) for i in random_samples ] 

        # now if there is a provided size we subsample
        if size is not None: 
            if size <= 1:
                random_samples = [ np.random.choice( i, size = size * i.size ) for i in random_samples ]
            elif size > 1:
                random_samples = [ np.random.choice( i, size = size ) for i in random_samples ]
            else:
                raise ValueError( f"size needs to be a fraction in [0,1] or an integer > 1 (got size = {size})" )

        self._baseline = random_samples

        # Alright, we currently have the entire sets of pre-timeframes for the baseline and we
        # will use them as they are completely to use for the baseline comparison. 
        # With the code below we compare a sub-sampled versions thereof. Long story short,
        # it works also pretty well with sub-sampled versions as well...
        # import statsmodels.api as sm
        # from matplotlib import colors

        # fig, ax = plt.subplots( 2,3 ) 
        # for r in random_samples:
        #     ax[0,0].plot( r )
        #     r1 = r.reshape( r.size )
        #     ax[0,1].hist( r1, bins = 50 )

        #     sm.qqplot( r1, ax = ax[0,2], alpha = 0.3, line = "s", color = list(colors.cnames.values())[ int(np.random.randint(low = 0, high = 10, size = 1))]  )
        # random_samples = [ np.random.choice( r.reshape(r.size), size = size, replace = False ) for r in random_samples ]
        
        # for r in random_samples:
        #     ax[1,0].plot( r )
        #     r1 = np.reshape( r, r.size )
        #     ax[1,1].hist( r1, bins = 50 )
            
        #     sm.qqplot( r1, ax = ax[1,2], alpha = 0.3, line = "s", color =  list(colors.cnames.values())[ int(np.random.randint(low = 0, high = 10, size = 1))] )
        # # ax[1].hist( np.reshape( random_samples, random_samples.size)  )
        # plt.tight_layout()
        # plt.show()


    def pvalues( self, event1 : int, event2 : int = None ):
        """
        Gets the p-value np.ndarray for each signal timepoint from a comparison of 
        either two separate event types or one event with its baseline. 

        Parameters
        ----------
        event1 : int
            The numeric event identifier of the (first) signal to get.
            If `None` is provided, the entire dictionary of pvalues is returned.

        event2 : int
            The numeric event identifier of the (second) signal from the comparison to get.
            If `None` is provided then the first signals comparison to it's baseline will be 
            returned (if baseline comparison was performed).
        
        Returns
        -------
        pvalues : np.ndarray or dict
            An np.ndarray of p-values from a given comparison.
        """

        if event1 is None: 
            return self._pvalues

        if event2 is None:
            key = (event1, event1)
        else: 
            key = (event1, event2)
        pvalues = self._pvalues.get( key, None )
        return pvalues 

    @property
    def events( self ):
        """
        Returns 
        -------
        list
            A list of all different event types from from the loaded metadata.
        """
        return list( self._events.keys() )

    @property
    def timeframe( self ):
        """
        Returns
        -------
        tuple   
            The used timeframe for event data extraction.
            This consists of the pre-trigger and post-trigger
            time offsets in seconds.
        """
        return ( self._start_sec, self._stop_sec )



    def summary(self,
                x_scale:float,
                y_scale:float,
                significance_level:float = 0.05,
                output:str = None,
                **kwargs ) -> None:
        """
        Performs pair-wise T-Tests to compare extracted event data 
        (automatically extracts data for all events if no events were extracted yet). 
        Results are summarised in a figure. Individual signals are plotted
        on the diagonal by their mean signal accross replicates with indicated SEM.
        On non-diagonal plots pair-wise comparisons between two signals (one "horizontal"
        and one "vertical" ) are shown. Regions of significant differences are shaded.

        Parameters
        ----------

        x_scale : float
            A scaling factor to adjust the data's x-value range. 
            E.g. `x_scale = 1000` to adjust the time-scale to milliseconds.
        
        y_scale : float
            A scaling factor for the data's y-value range.
            E.g. `y_scale = 1000` to adjust the signal-scale to millivolts.

        significance_level : float
            The threshold for accepting a signal difference as significant.
            Default is `0.05`.
        
        output : str
            The output filename to save the summary figure into.

        **kwargs
            Any additional keyword arguments to pass to `EEGData.extract` in case
            no event data has been extracted yet.
        """

        # extract the event data if not yet done already
        start_sec = kwargs.pop( "start_sec", self._start_sec )
        stop_sec = kwargs.pop( "stop_sec", self._stop_sec )
        if self._data is None: 

            self.extract( start_sec = start_sec, stop_sec = stop_sec, **kwargs )
            self.baseline() 

        data = list( self._data ) 
        signals = list(self._events.keys())
        n = len(data)

        # generate a new figure
        figsize = kwargs.pop( "figsize", ( 3*n,2*n ) )
        fig, ax = plt.subplots(n,n, figsize = figsize )

        # setup a baseline reference, either with the computed
        # baselines or None ...
        baseline = self._baseline if self._baseline is not None else [ None for i in range(n) ]
    
        # now first plot the individual signals
        # on their own on diagonal plots
        for i in range(n):

            # only the last subplot should make a legend
            make_legend = i == n-1 
            p = plot_signal(
                    data[i], 
                    self.sampling_frequency, 
                    start_sec, stop_sec, 
                    x_scale, y_scale,
                    baseline = baseline[i],
                    make_legend = make_legend,
                    significance_level = significance_level,
                    ax = ax[i,i] )
                
            ax[i,i].set_title(f"Signal {signals[i]}")

            # if we got a baseline to compare to we also want to 
            # store the resulting p-values
            if p is not None: 
                self._pvalues[ (i,i) ] = p 

            # hide all "left-over" subplots from the layout
            # i.e. hide the upper-right half of the figure...
            for a in ax[ i, i+1: ]:
                a.axis("off")

        # now make pair-wise comparisons between two signals
        # plotting the results on the lower-left half of the 
        # figure...
        for i,j in [(i,j) for i in range(n) for j in range(i)]:

            # only the last plot shall make a legend
            make_legend = i == n-1 and j == i-1 

            p = difference_plot( 
                                data[i], 
                                data[j], 
                                self.sampling_frequency, 
                                start_sec, stop_sec, 
                                significance_level, 
                                x_scale, y_scale,
                                make_legend = make_legend,
                                ax = ax[i,j]
                            )
            ax[i,j].set_title(f"Signals: {signals[j]} vs {signals[i]}")

            # we also want to store the resulting p-values of the 
            # signal comparison
            self._pvalues[ ( signals[j],signals[i] ) ] = p

        fig.tight_layout()
        
        if output is None:
            plt.show()
            return fig
        plt.savefig(output, bbox_inches = "tight" )


    def _extract_window(self, start_sec, stop_sec, event_type):
        """
        Extracts a set of time-frame windows from the data 
        and returns them as a numpy ndarray.
        """

        # first adjust the start and end to 
        # match the sampling frequency
        start_frame, stop_frame = self._adjust_timesteps(start_sec, stop_sec)

        # next generate a set of slices for the EEG data around the timepoints for
        # the events
        firing_slices = [
                            slice( event[0]+start_frame, event[0]+stop_frame ) 
                            for event in self._events_data 
                            if event[1] == event_type
                    ]

        # now get the actual data of the event
        data = [ self.signal[ slice ] for slice in firing_slices ]
        data = np.array( data )
        return data

    def _adjust_timesteps(self, start_sec, stop_sec):
        """
        Adjusts time steps / time points with the used recording frequency,
        to match the indices within the data.
        """
        start_frame = int( start_sec * self.sampling_frequency )
        stop_frame = int( stop_sec * self.sampling_frequency )
        return start_frame,stop_frame

    def _set_n_events(self) -> None:
        """
        Sets up a dictionary of the different event types
        found in the events data.
        """

        event_types = {event[1] for event in self._events_data}
        self._events = {event_type: len([event for event in self._events_data if event[1] == event_type]) for event_type in event_types}            

    def _check_sanity(self, signal_path, event_path, sampling_frequency):
        """
        Checks if valid data inputs were provided
        """

        # check if input files exist
        if not os.path.isfile(signal_path):
            raise FileNotFoundError( f"The signal datafile could not be found at '{signal_path}'!" )
        
        if not os.path.isfile(event_path):
            raise FileNotFoundError( f"The event datafile could not be found at '{event_path}!" )
        
        # check if the datafiles conform to suppored filetypes
        fname = os.path.basename(signal_path)
        if not any( [ fname.endswith(suffix) for suffix in supported_filetypes ] ):
            suffix = fname.split(".")[-1]
            raise TypeError( f"The signal datafile could not be interpreted ('.{suffix}'), only {supported_filetypes} files are supported!" )

        fname = os.path.basename(event_path)
        if not any( [ fname.endswith(suffix) for suffix in supported_filetypes ] ):
            suffix = fname.split(".")[-1]
            raise TypeError( f"The event datafile could not be interpreted ('.{suffix}'), only {supported_filetypes} files are supported!" )
        
        # check if the frequency is a positive number
        if not sampling_frequency > 0:
            raise ValueError( f"Sampling frequency {sampling_frequency} must be a positive value" )

    def _read_npy( self, filepath ) -> np.ndarray:
        """
        Reads data from npy files
        """
        data = np.load(filepath)
        return data 

    def _read_datafile( self, filepath) -> np.ndarray: 
        """
        Reads data from tsv, csv, and txt files
        """

        # first get the filetype 
        suffix = self._filesuffix(filepath)

        # now get the corresponding delimiter
        delimiters = {
                        "csv" : self._csv_delimiter( filepath ),
                        "tsv" : "\t",
                        "txt" : " "
                    }
        delimiter = delimiters[ suffix ]
        
        # now read the file
        data = pd.read_csv( filepath, header = None, sep = delimiter )

        # convert to numpy ndarray
        data = data.to_numpy()
        data = np.squeeze( data )
        return data 

    def _filesuffix(self, filepath):
        """
        Returns the suffix from a filepath
        """
        suffix = os.path.basename( filepath )
        suffix = suffix.split(".")[-1]
        return suffix

    def _csv_delimiter(self, filepath):
        """
        Checks if a csv file is , or ; delimited and returns the 
        correct delmiter to use...
        """
        # open the file and read 
        with open( filepath, "r" ) as f:
            content = f.read()

        # check if a semicolon is present
        # if so, we delimit at ; 
        has_semicolon = ";" in content
        delimiter = ";" if has_semicolon else ","

        return delimiter

def main():
    """
    The main function called through the CLI

    Example Usage
    -------------
    python3 ./EEGToolkit/EEGData.py \
        --eeg_path "data/eeg.npy" \
        --event_path "data/events.npy" \
        --sampling_frequency 500 \
        --p_value 0.05 \
        --start_sec -0.3 \
        --stop_sec 1.0 \
        --x_scale 1000 \
        --y_scale 1000 \
        --output "./test.png"
    """

    descr1 = """

-----------------------------------------------------
▒█▀▀▀ ▒█▀▀▀ ▒█▀▀█ ▀▀█▀▀ █▀▀█ █▀▀█ █░░ ▒█░▄▀ ░▀░ ▀▀█▀▀
▒█▀▀▀ ▒█▀▀▀ ▒█░▄▄ ░▒█░░ █░░█ █░░█ █░░ ▒█▀▄░ ▀█▀ ░░█░░
▒█▄▄▄ ▒█▄▄▄ ▒█▄▄█ ░▒█░░ ▀▀▀▀ ▀▀▀▀ ▀▀▀ ▒█░▒█ ▀▀▀ ░░▀░░
-----------------------------------------------------

This script takes in two data files of EEG signal data and accompanying event-trigger metadata. It performs intra- and inter-signal type comparisons using pair-wise T-Tests over the time-series, highlighting significantly different stretches and producing a summary figure. 
    """
    descr2 = f"""

Input Data
----------
Accepted input file types are {supported_filetypes}. The EEG-signal datafile must specify a 1D array of measurements, while the trigger metadata file must specify a 2D array (2 columns) of trigger time points and event classifier labels (numerically encoded). 
    """
    
    parser = argparse.ArgumentParser( prefix_chars = "-", 
    formatter_class=argparse.RawDescriptionHelpFormatter,description = descr1, epilog = descr2 )
    parser.add_argument(
                            "--eeg_path", "--eeg", 
                            type=str,
                            help = f"A file containing EEG signal data. Supported filetypes are {supported_filetypes}" 
                    )
    parser.add_argument(
                            "--event_path", "--event", 
                            type=str,
                            help = f"A file containing event metadata for the signal file. Supported filetypes are {supported_filetypes}"
                    )
    parser.add_argument(
                            "--output", "-o", 
                            type=str, default = None,
                            help = "An output file into which the output summary figure should be saved. If none is provided (default) the Figure will simply be shown."
                    )
    parser.add_argument(
                            "--sampling_frequency", "--freq", "-f", 
                            type=float,
                            help = "The frequency at which the EEG signal data was recorded (in Hertz)."
                    )
    parser.add_argument(
                            "--p_value", "-p", 
                            type=float, default = 0.05, 
                            help = "The significance threshold at which to accept two signals being significantly different at a T-Test comparison. Default is 0.05."
                    )
    parser.add_argument(
                            "--start_sec", "--start", "-s", 
                            type=float,
                            help = "The upstream time-padding for event extraction (in seconds)."
                    )
    parser.add_argument(
                            "--stop_sec", "--stop", "-e", 
                            type=float,
                            help = "The downstream time-padding for event extraction (in seconds)."
                    )
    
    parser.add_argument(
                            "--baseline", "-b", 
                            type=bool, default = True,
                            help = "Perform baseline comparison for each event type using the same significance threshold as used for inter-signal comparisons. Will be performed by default."
                    )
    parser.add_argument(
                            "--x_scale", "-x", 
                            type=float, default = 1000,
                            help = "A scaling factor for the time-scale (x-values) from seconds to some other unit. Default is 1000 (= milliseconds)."
                    )
    parser.add_argument(
                            "--y_scale", "-y", 
                            type=float, default = 1000,
                            help = "A scaling factor for the signal-scale (y-values) from volts to some other unit. Default is 1000 (= millivolts)."
                    )
    
    parser.add_argument(
                            "--viewer", "-i", 
                            action="store_true",
                            default = False,
                            help = "Open the EEGToolKit Viewer GUI in a web browser."
                    )

    args = parser.parse_args()
    
    # if the viewer is being called then we want to just open the 
    # viewer and nothing else
    if args.viewer:
        # first we need to get the relative location of the main.
        # py file within the package. 
        directory = os.path.dirname( 
                                            inspect.getfile( plot_signal ) 
                                        )
        directory = os.path.dirname( directory )
        main_file = f"{directory}/main.py"

        # then we call the web interface
        print( "Starting the \033[94mEEGToolKit \033[96mViewer" )
        subprocess.run( f"streamlit run {main_file}", shell = True )
    
    else: 

        # the main program (reading datafiles, extracting, and summarizing)
        try:
            data = EEGData(args.eeg_path, args.event_path, args.sampling_frequency)
            data.extract( args.start_sec, args.stop_sec )
            if args.baseline:
                data.baseline()
            data.summary(
                            significance_level = args.p_value,
                            x_scale = args.x_scale,
                            y_scale = args.y_scale,
                            output = args.output
                        )

            if args.output is not None: 
                print( f"Output saved successfully to: '{args.output}'" )
        except FileNotFoundError as e:
            print(e)
            return
        except TypeError as e:
            print(e)
            return
        except ValueError as e:
            print(e)
            return

if __name__ == "__main__":

    test_mode = False
    if not test_mode:
        main()
    else: 
        print( "Running in Test Mode" )

        eeg = "./data/eeg.npy"
        events = "./data/events.npy"

        e = EEGData( eeg, events, 500 )
        e.extract( -0.3, 1 )
        e.baseline()
        e.summary( 1000, 1000, output = "./test.pdf" )
        plt.show()