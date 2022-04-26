"""
TODO
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from EEGStats import _plot_, _difference_plot_

class EEGData():
    """
    TODO
    """

    def __init__(self,
                 signal_path:str,
                 event_path:str,
                sampling_frequency:float) -> None:
        """
        TODO
        """

        assert os.path.isfile(signal_path), f"{signal_path} is not a file!"
        assert os.path.isfile(event_path), f"{event_path} is not a file!"
        assert sampling_frequency > 0, f"Sampling frequency {sampling_frequency} should be a positive value"

        self.signal = np.load(signal_path)
        self.n_frames = len(self.signal)
        self.events = np.load(event_path)
        self.set_n_events()
        self.sampling_frequency = sampling_frequency

    def extract(self,
                start_sec:float,
                stop_sec:float,
                p_event_type=None) -> np.ndarray:
        """
        TODO
        """

        if p_event_type is None:
            return (self.extract(start_sec, stop_sec, event_type) for event_type in self.n_events.keys())
        if isinstance(p_event_type, (np.ndarray, list)):
            return (self.extract(start_sec, stop_sec, event_type) for event_type in p_event_type)

        start_frame = int(start_sec*500)
        stop_frame = int(stop_sec*500)

        firing_slices = [slice(event[0]+start_frame, event[0]+stop_frame) for event in self.events if event[1] == p_event_type]

        return np.array([self.signal[slice] for slice in firing_slices])

    def summary(self,
                significance_level:float,
                start_sec:float,
                stop_sec:float,
                x_scale:float,
                y_scale:float,
                output:str) -> None:
        """
        TODO
        """

        data = list(self.extract(start_sec, stop_sec))
        signals = list(self.n_events.keys())
        n = len(data)

        fig, ax = plt.subplots(n,n)

        for i in range(n):
            _plot_(ax[i,i], data[i], self.sampling_frequency, start_sec, stop_sec, x_scale, y_scale)
            ax[i,i].set_title(f"Signal {signals[i]}")

        for i,j in [(i,j) for i in range(n) for j in range(i)]:
            _difference_plot_(ax[i,j], data[i], data[j], significance_level, self.sampling_frequency, start_sec, stop_sec, x_scale, y_scale)
            ax[i,j].set_title(f"Significant {signals[j]} vs {signals[i]}")

        fig.tight_layout()

        if output is None:
            plt.show()
            return

        plt.savefig(output)

    def set_n_events(self) -> None:
        """
        TODO
        """

        event_types = {event[1] for event in self.events}
        self.n_events = {event_type: len([event for event in self.events if event[1] == event_type]) for event_type in event_types}            

def main():
    """
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
    parser.add_argument("--eeg_path", type=str, required=True)
    parser.add_argument("--event_path", type=str, required=True)
    parser.add_argument("--sampling_frequency", type=float, required=True)
    parser.add_argument("--p_value", type=float, required=True)
    parser.add_argument("--start_sec", type=float, required=True)
    parser.add_argument("--stop_sec", type=float, required=True)
    parser.add_argument("--x_scale", type=float, required=True)
    parser.add_argument("--y_scale", type=float, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    
    EEGData(args.eeg_path, args.event_path, args.sampling_frequency).summary(args.p_value,
                                                                             args.start_sec,
                                                                             args.stop_sec,
                                                                             args.x_scale,
                                                                             args.y_scale,
                                                                             args.output)

if __name__ == "__main__":
    main()