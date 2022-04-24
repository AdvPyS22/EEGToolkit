"""
    TODO
"""

import os
import numpy as np

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

    def set_n_events(self) -> None:
        """
            TODO
        """

        event_types = {event[1] for event in self.events}
        self.n_events = {event_type: len([event for event in self.events if event[1] == event_type]) for event_type in event_types}
        