
Data for: Theta rhythmicity governs the timing of behavioural and hippocampal responses in humans specifically during memory-dependent tasks

Marije ter Wal, Juan Linde Domingo, Julia Lifanov, Frederic Roux, Luca Kolibius, Stephanie Gollwitzer, Johannes Lang, Hajo Hamer, David Rollings, Vijay Sawlani, Ramesh Chelvarajah, Bernhard Staresina, Simon Hanslmayr, Maria Wimber

Correspondence: Marije ter Wal - m.j.terwal@bham.ac.uk

Last updated: 30/10/2020

--------------------------------------------------

* HOW TO USE

This data is available for public use under a Creative Commons CC BY-NC-SA license, meaning the data can be used for non-commercial purposes, leading to publicly available results or publications, as long as the original data repository is referenced.

Please refer to the data repository as: 
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, R., Staresina, B., Hanslmayr, S., Wimber, M., (2020), Data for: Theta rhythmicity governs the timing of behavioural and hippocampal responses in humans specifically during memory- dependent tasks. figshare. Collection. DOI: 10.6084/m9.figshare.c.5192567

Please cite the manuscript presenting the results from this dataset as: 
Ter Wal, M., Linde Domingo, J., Lifanov, J., Roux, F., Kolibius, L., Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, R., Staresina, B., Hanslmayr, S., Wimber, M., (in preparation), Theta rhythmicity governs the timing of behavioural and hippocampal responses in humans specifically during memory- dependent tasks.


--------------------------------------------------

* SOFTWARE 

All data in this dataset are stored as .MAT files and can be opened using MATLAB (The Mathworks). 
Please contact the corresponding author if you require a different data format.

Scripts and functions that were used to analyse the data and obtain the results in Ter Wal et al., (in preparation) can be found here:
https://github.com/marijeterwal/behavioral-oscillations


--------------------------------------------------

* DESCRIPTION OF THE DATA

The dataset contains the following data:

1) behavioral data from 95 participants performing visual tasks (experiments 1-4)
2) behavioral data from 226 participants performing memory tasks (experiments 5-13; including 10 iEEG patients)
3) PPC data from hippocampal microwire recordings from 10 iEEG patients, as well as 100 time shuffled and 100 trial-label shuffled PPC datasets

The PPC data have been fully anonymized to protect the identity of the patients: information that could reveal the implantation scheme of specific patients has been removed. To this end, the original channel labels have been removed and replaced by bundle IDs. 


----- BEHAVIORAL DATA -----

* FILE NAME CONVENTION

BehavioralData_Experiment[experimentID].mat

Experiment IDs are the same as those used in the manuscript (see Tables S1 & S2 for details). 


* FIELD NAME CONVENTION

All behavioral data files contain a table 'dataBehavior', with every row of the table representing a trial (i.e. one button press) from one participant. Data from all participants of an experiment were combined.
The columns of the table contain all relevant information about the trial:
- subjID: ID of the subject, with the first digit(s) representing the experiment ID and the last two digits representing a unique identifier of the subject;
- ExpID: experiment ID, same as in the file name;
- TrialID: unique identifier of the object (visual) or cue-object pair (memory) for this subject. For memory tasks, encoding, retrieval and catch trials for the same cue-object pair have the same trialID;
- TrialType: 'visual', for visual tasks, and 'encoding','retrieval','ret&catch' or 'catch' for memory tasks;
- RTs: reaction time in seconds
- Acc: performance on the catch questions related to this trial: 0 for incorrect and 1 for correct. Where two catchs questions were asked, a trial was only deemed incorrect if both questions were answered incorrectly.


----- PPC DATA -----

* FILE NAME CONVENTION

PPC_SubjID-[subjID]_[taskphase]_[acc]_[datatype].mat

Where:
- subjID is an unique identifier of the iEEG patient. This identifier is also used in the behavioral data files
- taskphase: encoding, retrieval or catch trials
- acc: correct or incorrect trials
- datatype: original data (no datatype given), timeshuffle or trialshuffle


* FIELD NAME CONVENTION

All PPC data files contain a structure 'PPC' containing the PPC and baseline data as well as the necessary metadata. The structure contains two or three sets of data, locked to difference events in the trial, such as cue or stimulus onset and response (i.e. button press).
The structure's fields contain the following information:

- taskphase: 'encoding', 'retrieval', or 'catch'; same as the file name;
- subjID: unique subject identifier; same as in the file name;
- eventinfo: description of the types of event the data was locked to, for example 'cue', when locked to cue onset, 'response', when locked to the patient's button press, or 'catchonset', when locked to appearance of the catch question;
- bundle: array with a bundle ID for each of the recorded microwires. Channels with the same bundle ID came from the same bundle and hence share a reference;
- fsample: sampling rate of the PPC data (number of data points per second)
- freq: frequency axis corresponding to the PPC data (in Hz);
- time: cell containing time axes in seconds corresponding to the PPC data. Every cell contains the axis relative to the event-locking given in field eventinfo, with negative values meaning time before the event, and positive time point after the event;
- PPC: cell containing the PPC data locked to the different events in the trial as given in eventinfo. Within a cell, an array is given of size: channels x frequencies x timepoints;
- PPC_baseline: Pre-cue baseline PPC data, same convention as PPC.






