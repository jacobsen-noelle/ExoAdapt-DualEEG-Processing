# ExoAdapt-DualEEG-Processing
This pipeline performs preprocessing on dual-electrode electroencphalography data and source analysis to investigate gait-related changes in electrocortical activity during adaptation to a powered ankle exoskeleton. 
All local data processing is written for MATLAB and relies heavily on the EEGLAB toolbox (see Toolboxes.txt for additional external tools or EEGLAB plug-ins used). Some of the steps (independent component analysis, custom head modeling, source localization) are designed to be used with supercomputing resources.

These are the steps for dual-EEG processing:
1) Dual-EEG Preprocessing
   This step uses dual-electrode EEG and neck electromyography to isolate and remove motion and muscle artifacts from EEG [1]. It also rejects noisy channels and time windows.
2)  Add AMICA results
   After EEG preprocessing, Adaptive Mixture Indpendent Component Analysis is run on supercomputing resources extract independent components from "clean" scalp EEG data (see AMICA folder). This steps pulls results from the supercompter and adds the ICA matrix to the EEG dataset.
3) Add DIPFIT results
   Source localization using DIPFIT and custom head models is performed using supercomputing resources. The *custom-head-model* folder provides a pipeline for creating 5-layer volume conducting head models with the finite element method, following the fieldtrip-SIMBIO pipleine.
The *source-localization* folder provides code for running DIPFIT with the custom head models and implements parallel processing to improve computation speed. This step adds source model results from DIPFIT located on the supercomputer back into the local EEG dataset
4) IC rejection
   Automatic rejection of non-brain independent components based on a spatial and temporal features
5) Epoching and Timewarping
   Epochs and timewarps data based on gait events
6) Group Analysis
   Precomputes measures (e.g. power spectral density, event-related spectral pertubations) and performs component clustering
7) Plot results
   Perform statistical analyses and plot results

[1] N. Richer, R. J. Downey, W. D. Hairston, D. P. Ferris, and A. D. Nordin, “Motion and Muscle Artifact Removal Validation Using an Electrical Head Phantom, Robotic Motion Platform, and Dual Layer Mobile EEG,” IEEE Transactions on Neural Systems and Rehabilitation Engineering, vol. 28, no. 8, pp. 1825–1835, Aug. 2020, doi: 10.1109/TNSRE.2020.3000971. 
[2] Ryan J Downey, Daniel P Ferris, "iCanClean Removes Motion, Muscle, Eye, and Line-Noise Artifacts from Phantom EEG," Sensors. Oct. 2023
