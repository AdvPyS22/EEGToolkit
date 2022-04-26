"""

This module provides a data class `EEGData` to work with EEG signal data for event-reaction-time delay experiments.
It works with two separate input datafiles, one storing the EEG signal itself as a 1D array, and one 
describing event metadata as a 2D array, describing both the timepoints and the type of event in two columns.

Supported file types are:
- `npy`
- `txt`     ( space-separated for events datafiles )
- `tsv`
- `csv`     (both `,` and `;` separated )   

### Example Usage
To use this module for data analysis, only three steps are necessary,
(1st) setup of the `EEGData` object, (2nd) event data extraction, and (3rd)
data summary (which performs signal comparison).

```
# setting up the EEGData with some datafiles
eeg = EEGData( eeg_path = "data/eeg.npy", event_path = "data/events.npy", sampling_frequency = 500 )

# extract the events data
data.extract( start_sec = -0.3 , stop_sec = 1 )

# summarize and pair-wise compare event-signal types.
data.summary(
                significance_level = 0.05,
                x_scale = 1000,
                y_scale = 10000,
            )
```

### CLI 
This module additionally offers a CLI to directly call the full analysis procedure from the terminal.

```
python3 EEGData.py \
                    --eeg_path "./data/eeg.npy" \
                    --event_path "./data/events.npy" \
                    --sampling_frequency 500 \
                    --p_value 0.05 \
                    --start_sec -0.3 \
                    --stop_sec 1.0 \
                    --x_scale 1000 \
                    --y_scale 10000 \
                    --output "./test_output.pdf"
```

"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from EEGStats import plot_signal, difference_plot

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
        self.n_frames = len(self.signal)
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
            self.events = events


    def extract(self,
                start_sec:float,
                stop_sec:float,
                event_type : ( int or tuple or list or np.ndarray ) = None) -> np.ndarray:
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
            events_to_extract = self.n_events.keys()

            # extract each type from the loaded data
            data = [ 
                        self.extract(start_sec, stop_sec, etype) 
                        for etype in events_to_extract 
                ]

            self._data = data
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
            return data

        # now the part for extracting only a 
        # single event type data
        
        # first adjust the start and end to 
        # match the sampling frequency 
        start_frame = int( start_sec * self.sampling_frequency )
        stop_frame = int( stop_sec * self.sampling_frequency )

        # next generate a set of slices for the EEG data around the timepoints for
        # the events
        firing_slices = [
                            slice( event[0]+start_frame, event[0]+stop_frame ) 
                            for event in self.events 
                            if event[1] == event_type
                    ]

        # now get the actual data of the event
        data = [ self.signal[ slice ] for slice in firing_slices ]
        data = np.array( data )

        self._data = data

        # store start and stop sec values
        # for later use in summary()
        self._start_sec = start_sec
        self._stop_sec = stop_sec

        return data

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
        if self._data is None: 

            start_sec = kwargs.pop( "start_sec", self._start_sec )
            stop_sec = kwargs.pop( "stop_sec", self._stop_sec )
            self.extract( start_sec = start_sec, stop_sec = stop_sec, **kwargs )

        data = list( self._data ) 
        signals = list(self.n_events.keys())
        n = len(data)

        start_sec, stop_sec = self._start_sec, self._stop_sec

        # generate a new figure
        fig, ax = plt.subplots(n,n)

        # now first plot the individual signals
        # on their own on diagonal plots
        for i in range(n):

            # only the last subplot should make a legend
            make_legend = i == n-1 
            plot_signal(
                    data[i], 
                    self.sampling_frequency, 
                    start_sec, stop_sec, 
                    x_scale, y_scale,
                    make_legend = make_legend,
                    ax = ax[i,i] )
                
            ax[i,i].set_title(f"Signal {signals[i]}")

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

            difference_plot( 
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

        fig.tight_layout()
        
        if output is None:
            plt.show()
            return fig
        plt.savefig(output, bbox_inches = "tight" )

    def _set_n_events(self) -> None:
        """
        Sets up a dictionary of the different event types
        found in the events data.
        """

        event_types = {event[1] for event in self.events}
        self.n_events = {event_type: len([event for event in self.events if event[1] == event_type]) for event_type in event_types}            

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

    parser = argparse.ArgumentParser(prefix_chars='-')
    parser.add_argument(
                            "--eeg_path", "--eeg", 
                            type=str, required=True, 
                            help = f"A file containing EEG signal data. Supported filetypes are {supported_filetypes}" 
                    )
    parser.add_argument(
                            "--event_path", "--event", 
                            type=str, required=True,
                            help = "A file containing event metadata for the signal file. Supported filetypes are {supported_filetypes}"
                    )
    parser.add_argument(
                            "--output", "-o", 
                            type=str, default = None,
                            help = "An output file into which the output summary figure should be saved. If none is provided (default) the Figure will simply be shown."
                    )
    parser.add_argument(
                            "--sampling_frequency", "--freq", "-f", 
                            type=float, required=True,
                            help = "The frequency at which the EEG signal data was recorded (in Hertz)."
                    )
    parser.add_argument(
                            "--p_value", "-p", 
                            type=float, default = 0.05, 
                            help = "The significance threshold at which to accept two signals being significantly different at a T-Test comparison. Default is 0.05."
                    )
    parser.add_argument(
                            "--start_sec", "--start", "-s", 
                            type=float, required=True,
                            help = "The upstream time-padding for event extraction (in seconds)."
                    )
    parser.add_argument(
                            "--stop_sec", "--stop", "-e", 
                            type=float, required=True,
                            help = "The downstream time-padding for event extraction (in seconds)."
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

    args = parser.parse_args()
    
    # the main program (reading datafiles, extracting, and summarizing)
    data = EEGData(args.eeg_path, args.event_path, args.sampling_frequency)
    data.extract( args.start_sec, args.stop_sec )
    data.summary(
                    significance_level = args.p_value,
                    x_scale = args.x_scale,
                    y_scale = args.y_scale,
                    output = args.output
                )

    if args.output is not None: 
        print( f"Output saved successfully to: '{args.output}'" ) 

if __name__ == "__main__":
    main()