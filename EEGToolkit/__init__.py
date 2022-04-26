"""

This package is designed for quick analysis of EEG data for 
event reaction-time delay experiments. It works with separate data
files for EEG signal data (in volts) and event metadata (describing event 
timepoints). 

Supported file types are:
- `npy`
- `txt`     ( space-separated for events datafiles )
- `tsv`
- `csv`     (both `,` and `;` separated )   

### Example Usage
To use this module for data analysis, only three steps are necessary,
(1st) setup of the `EEGData` object, (2nd) event data extraction, and (3rd)
data summary (which performs signal comparison).

```python
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

```bash
python3 EEGData.py --eeg "./data/eeg.npy" --event "./data/events.npy" --sampling_frequency 500 --p_value 0.05 --start_sec -0.3 --stop_sec 1.0 --x_scale 1000 --y_scale 10000 --output "./test_output.pdf"
```
"""
from .EEGData import *
from .EEGStats import *