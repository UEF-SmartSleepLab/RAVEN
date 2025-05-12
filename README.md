# RAVEN

A software for automated sleep microstructure analysis from electroencephalogram (EEG).

The RAVEN software is a MATLAB based software which provides automatic analysis of sleep spindles, K-complexes, delta waves, and cyclic alternating pattern (CAP) from raw EEG signal. 
RAVEN accepts raw EEG signals in European Data Format or SleepLab format and is operated through a user friendly GUI.

# How to use?
To use RAVEN you will need MATLAB and you have to download all three files (Main_code, RAVEN, and RAVEN_UI) into a same directory. RAVEN is then initialised by running the Main_code through MATLAB which opens the GUI. 
More detailed instructions are written down in the Main_code file.

# Version 1.0
This release includes the initial implementation of the core algorithm. Please note that this version has not yet been validated against expert-scored reference data. Validation is currently in progress, and results—along with any necessary adjustments—will be incorporated in the next version update.

Users are encouraged to do case-specific validation analyses and fine-tuning of the thresholds prior to the utilisation of the RAVEN. This is a foundational release and will evolve as validation and benchmarking continue.
