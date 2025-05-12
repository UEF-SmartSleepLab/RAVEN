%% Full night EEG analysis from .edf/.npy file


% PLEASE FIND THE INFORMATION BELOW FOR YOUR CONVENIENCE.


% INTRO:

% The RAVEN performs a comprehensive analysis of the 
% entire night's electroencephalography (EEG) signal for sleep 
% microstructures. Analysis is performed for spindles, k-complexes, 
% delta-waves (Slow Wave Sleep, SWS), and cyclic altering pattern (CAP).

% The function analyses these microstructures regardless of the sleep
% stage. Therefore, hypnogram is not needed for the RAVEN to operate
% correctly. 
 
% Analysis for spindles, K-complexes, and delta waves require atleast one 
% EEG channel. CAP analysis is based on the simultaneous analysis of two 
% different EEG channel signals to provide the most accurate event scoring 
% possible.



% KNOW YOUR SIGNALS:

% The function filters and analyses the EEG data from the original signal.
% The signal must be in European Data Format (EDF) or SleepLab Format (SLF).
% Before analysing, you need to be aware of the signal names in your data
% as well as the sampling frequencies used in measurements of the data. 




% HOW TO USE (INPUT):

% The RAVEN analysis can be executed directly by running this section
% (Main_code) that opens the RAVEN GUI. All of the code files must be in
% same folder as the Main_code.m. This includes RAVEN.m and RAVEN_UI.m

% NOTE: For CAP analysis you will need TWO signals. No more, no less. 

% NOTE: If you are using SleepLab Format, use the folder that includes the
% subject(s) folder(s) as input directory. RAVEN will automatically opens
% and search the correct data from the subfolders.



% ANALYSIS (OUTPUT):

% The algorithm saves the scored data into a specified location as separate 
% files in .mat format. These files are stored in a folder created by 
% RAVEN based on the file name (e.g., 001_EEG_analysis). The saved file 
% names follow the format 'Patient_ID_event', such as '001_spindles'. 
% Each scored event is saved as a separate entry within these files. 
% To prevent overwriting existing files for the same patient, RAVEN appends
% a timestamp to the folder name. This timestamp represents the current 
% time when the analysis starts, formatted in hours, minutes, and seconds 
% (e.g., 001_EEG_analysis_12_05_30).

% EVENT INFORMATION

% Spindles, K-complexes, delta waves (SWS): 
% "File_ID", your patient ID.
% "name", alternative annotation name for different GUI (Value = "UNSURE").
% "event_name", name of the scored event (microstructure).
% "input_channel", channels' name used for analysis.
% "sampling_frequency", your sampling frequency used in analysis.
% "indices", direct indices from the beginning of your used data.
% "start_sec", event starting time from the beginning of the data (seconds).
% "duration", duration of the scored event from starting time.
% "amplitude", max amplitude of scored event (not for K-complex). 

% CAP:
% A phases, and cycles.
% A phase:
% "File_ID", your patient ID.
% "name", alternative annotation name for different GUI (Value = "UNSURE").
% "event_name", name of the scored event.
% "sampling_frequency", your sampling frequency used in analysis.
% "subtype", subtype of the A phase (A1, A2, A3).
% "start_sec", event starting time from the beginning of the data (seconds).
% "duration", duration of the scored event from starting time.
% "indices", direct indices from the beginning of your used data.


% Cycles:
% "File_ID", your patient ID.
% "name", alternative annotation name for different GUI (Value = "UNSURE").
% "event_name", name of the scored event.
% "sampling_frequency", your sampling frequency used in analysis.
% "start_sec", event starting time from the beginning of the data (seconds).
% "duration", duration of the scored event from starting time.
% "whole_cycle", phases belonging to the cycle.


% Sequences:
% "File_ID", your patient ID.
% "name", alternative annotation name for different GUI (Value = "UNSURE").
% "event_name", name of the scored event.
% "sampling_frequency", your sampling frequency used in analysis.
% "start_sec", event starting time from the beginning of the data (seconds).
% "duration", duration of the scored event from starting time.
% "Sequence", Whole sequence.


% By running RAVEN_UI a graphical user interface is opened that you can use
% to run the RAVEN analysis. When adding the text/values to the fields
% dont use space, only comma (ex. C3,C4). 
% DO NOT CLOSE THE INTERFACE DURING THE ANALYSIS!


% THRESHOLDS:
% Used thresholds can be modified through the GUI to be more suitable.

    RAVEN_UI;

