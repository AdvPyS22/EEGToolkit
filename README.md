# AP-EEG

## Project description
> ## EEG data analysis
> Students who choose this task ill be provided with the raw EEG recording of one channel, sampled at $500 [Hz]$ which was recorded from a participant presented with auditory stimuli. The students will also receive an events file describing when during the recording one out of two possible sounds were presented to the subject.
> ### Task 1
> Extract the time-period from $0.1 [sec]$ before to $1 [sec]$ after each event for the EEG signal. Store these time-periods in two arrays, one for each of the two event-/sound-types.
> ### Task 2
> Plot the average response for each type of event, with shading for the standard error. On the plot, also indicate $t=0$, the amplitude of the signal in microvolts (Y-axis), and the time relative to sound onset in $[ms]$ (X-axis).
> ### Task 3
> For each sound type, for each time-point after sound onset, test wheter it is significantly different across trials from the baseline. Also test at what time-points after sound onset the EEG signals for the two sound-types are different from each other. Produce one plot containing the average EEG signal of each sound type, along with lines indicating the significant time-periods.

## Milestones
> - **Part 1: Project management and first deliverable** (*Deadline: April 28*)
> - **Part 2: GitHub and second deliverable** (*Deadline: May 5*)
> - **Part 3: Classes and refactoring** (*Deadline: May 12*)
> - **Part 4: Unit tests and issues** (*Deadline: May 19*)
> - **Part 5: Third deliverable, presentation and virtual environment** (*Deadline: June 2*)

## Project roadmap

### Main Objective
The goal is to develop a small package with CLI to evaluate EEG data for reaction delays given different auditory stimuli.
While the outset data is given as a binary stimulus experiment, extending to an arbitrary number of different stimuli would be 
a desirable extention. 


### Implementation

#### General TODOs 
- [x] proper **Docstrings** and comments in general ...
- [x] re-factoring for private methods ... 
- [ ] requirements.txt
- [ ] setup.py 
    - [ ] (optional but kinda cool) adding EEGData to PYTHONPATH to allow direct CLI calling...
- [ ] `testpypi` distribution

#### CLI TODOs:
- [x] Core is pretty much finished already...
- [x] Change the `argparse` settings to allow defaults for parameters like x_scale to allow easier usage

#### Task TODOs:
##### Task 1: 
- [x] Event extraction is pretty much finished already 
> Possible additions: <br>
> - [x] split reading the files and extracting data, thereby allowing different filetypes to be processed. In case this is done, we might add support for `txt` / `tsv` / `csv` files for more general applicability of the software. 


##### Task 2: 
- [x] plotting is pretty much finished already...
- [ ] Add units to the figure axes labels
- [x] (optional) adjust color scheme ...

##### Task 3: 
- [ ] What about the intra-signal baseline-comparison ? 
- [x] signal-comparison is finished already...
