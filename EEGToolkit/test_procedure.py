"""
This script performs some unit testing.
"""

import unittest
from EEGToolkit.EEGData import EEGData

signal_path = 'data/eeg.npy'
event_path = 'data/events.npy'
sampling_frequency = 500
data = EEGData(signal_path, event_path, sampling_frequency)

# General class for performing some unit tests
class Testclass(unittest.TestCase):

    def test_adjust_timsteps(self):
        self.assertTupleEqual(data._adjust_timesteps(-0.5, 1), (-250, 500), "Start frame and stop frame should be -250 and 500")

    def test_set_nevents(self):
        self.assertListEqual(data.events, [0, 1], "The list should be [0, 1]")

    def test_suffix(self):
        self.assertTrue(data._filesuffix(signal_path) == 'npy', "The suffix should be 'npy'")
        self.assertTrue(data._filesuffix(event_path) == 'npy', "The suffix should be 'npy'")


if __name__ == "__main__":
    unittest.main()