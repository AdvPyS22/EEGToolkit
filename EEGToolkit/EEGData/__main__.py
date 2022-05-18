import sys
import os 

import subprocess
import inspect
import argparse
import matplotlib.pyplot as plt

abspath = os.path.abspath(__file__)
dname = os.path.dirname(os.path.dirname(abspath))
sys.path.append(dname)

from EEGStats import plot_signal
from EEGData import EEGData, supported_filetypes


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
        main_file = os.path.join(directory, "__main__.py")
        #main_file = f"{directory}/__main__.py"

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