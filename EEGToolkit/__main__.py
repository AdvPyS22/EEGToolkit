"""
This script defines an interactive 
web-app for the EEGData package.
"""
import sys
import os

import streamlit as st 
import pandas as pd

abspath = os.path.abspath(__file__)
dname = os.path.dirname(os.path.dirname(abspath))
sys.path.append(dname)

from EEGData import supported_filetypes
from auxiliary import Session, stEEGData

if __name__ == "__main__":
        session = Session()

        # =================================================================
        #                          Header Section
        # =================================================================


        st.markdown( 
                        """
        # EEGToolKit Viewer 
                        """, 
                )

        # =================================================================
        #                        Data Input Section
        # =================================================================

        file_container = st.container()
        signal_col, events_col = file_container.columns( 2 )

        # add a signal datafile
        signal_file = signal_col.file_uploader( 
                                                "EEG Signal Data",
                                                help = f"Upload here a datafile specifying a 1D array of EEG signal data.",
                                                type = supported_filetypes
                                        )
        if signal_file is not None: 
                session.add( "signal_filename", signal_file.name )
                session.add( "signal_file", signal_file, export = False )

        # add an events metadata file
        events_file = events_col.file_uploader(
                                                "Events Metadata",
                                                help = f"Upload here a datafile specifying a 2-column array of events metadata for the signal data file.",
                                                type = supported_filetypes
                                        )
        if events_file is not None:
                session.add( "events_filename", events_file.name )
                session.add( "events_file", events_file, export = False )


        # =================================================================
        #                      Data Processing Section
        # =================================================================

        # In case anybody is wondering why the weird column setup here. The answer
        # is that streamlit does not (yet) support nested columns nor stretching objects
        # over multiple columns. Hence, one column, one object, if they should "appear" 
        # in side by side, one above the other etc. layouts we need to have multiple
        # column blocks...

        controls = st.container()
        upper_ctrl_col1, upper_ctrl_col2 = controls.columns( (1, 3) )
        mid_ctrl_col1, mid_ctrl_col2, mid_ctrl_col3 = controls.columns( (2, 3, 3 ) )
        lower_ctrl_col1, lower_ctrl_col2, lower_ctrl_col3 = controls.columns( (2, 3, 3 ) )

        compute_baseline = upper_ctrl_col2.checkbox( 
                                                        "Compare Baseline",
                                                        help = "Compute the baseline of a given signal type / event and compares the signal to the baseline through positition-wise T-Tests.",
                                                        value = True
                                        )

        significance_level = upper_ctrl_col2.number_input( 
                                                "Significance Level", 
                                                help = "Choose a significance threshold level for signal-signal and signal-baseline comparisons",
                                                min_value = 1e-5, max_value = 1.0, value = 0.05, step = 1e-3, format = "%e"
                                        )
        session.add( "significance_level", significance_level )

        frequency = upper_ctrl_col2.number_input( 
                                                "Sampling Frequency", 
                                                help = "The sampling frequency in `Hertz` at which the signal was recorded",
                                                min_value = 1, max_value = 100000, value = 500, step = 100, format = "%d", 
                                        )
        session.add( "frequency", frequency )


        start_sec = mid_ctrl_col2.number_input( 
                                                "Upstream Buffer", 
                                                help = "The time buffer pre-event to extract (in seconds).",
                                                min_value = 1e-4, max_value = 10.0, value = 0.01, step = 0.001, format = "%e", 
                                        )
        session.add( "upstream_buffer", start_sec )

        stop_sec = mid_ctrl_col3.number_input( 
                                                "Downstream Buffer", 
                                                help = "The time buffer post-event to extract (in seconds).",
                                                min_value = 1e-4, max_value = 10.0, value = 1.0, step = 0.001, format = "%e", 
                                        )
        session.add( "downstream_buffer", stop_sec )

        xscale = lower_ctrl_col2.number_input( 
                                                "Timeframe Scale", 
                                                help = "A scaling factor to adjust x-axis time scale. E.g. `1000` to adjust data from seconds to milliseconds.",
                                                min_value = 1, max_value = 100000, value = 1000, step = 100, format = "%d", 
                                        )
        session.add( "timescale", xscale )

        yscale = lower_ctrl_col3.number_input( 
                                                "Signal Scale", 
                                                help = "A scaling factor to adjust y-axis signal scale. E.g. `1000` to adjust from volts to millivolts.",
                                                min_value = 1, max_value = 100000, value = 1000, step = 100, format = "%d", 
                                        )
        session.add( "signalscale", yscale )

        # slightly ugly splitting of the text, but the only way to adjust the text to the column layout I want...
        mid_ctrl_col1.markdown( "Once your setup is done, press the `Compute` button below to " )
        lower_ctrl_col1.markdown( "commence data evaluation." )
        run_button = lower_ctrl_col1.button( 
                                                "Compute",
                                                help = "Perform computations and output a figure",
                                        )

        # =================================================================
        #                       Computational Section
        # =================================================================

        figure_container = st.container()
        if run_button:
                if signal_file is None: 
                        st.error( "No Signal File was provided so far!" )
                        st.stop()
                if events_file is None:
                        st.error( "No Events File was provided so far!" )
                        st.stop()
                eegdata = stEEGData( 
                                        signal_path = signal_file, 
                                        event_path = events_file, 
                                        sampling_frequency = frequency 
                                )
                eegdata.extract( -start_sec, stop_sec )

                if compute_baseline:
                        eegdata.baseline()
                fig = eegdata.summary( 
                                        x_scale = xscale, 
                                        y_scale = yscale, 
                                        significance_level = significance_level 
                                )
                
                figure_container.pyplot( fig )