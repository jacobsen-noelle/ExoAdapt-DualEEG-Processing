% Author: Ryan Downey, Univeristy of Florida
function [EEG EEG_chans Noise_chans neckEMG_chans Other_chans All_chans badEEGch badEMGch badNoiseCh] = autoChRej(EEG)
%% define channels
All_chans = 1:EEG.nbchan;
getchantypes;
Other_chans = setdiff(All_chans, unique([EEG_chans, Noise_chans, neckEMG_chans]));

%% Calculate things
%EEGchRange = prctile(EEG.data(EEG_chans,:)',99.5) - prctile(EEG.data(EEG_chans,:)',0.5);
%EMGchRange= prctile(EEG.data(neckEMG_chans,:)',99.5) - prctile(EEG.data(neckEMG_chans,:)',0.5);
%NoiseChRange = prctile(EEG.data(Noise_chans,:)',99.5) - prctile(EEG.data(Noise_chans,:)',0.5);

% figure; stem(chRange);
% badChans = find(chRange>200);


EEG_SD = nanstd(EEG.data(EEG_chans,:)');
EMG_SD = nanstd(EEG.data(neckEMG_chans,:)');
Noise_SD = nanstd(EEG.data(Noise_chans,:)');
thres = 3;
badEEGch = EEG_chans( find(EEG_SD > thres*median(EEG_SD)) );
badEMGch = neckEMG_chans( find(EMG_SD > thres*median(EMG_SD)) );
badNoiseCh = Noise_chans( find(Noise_SD > thres*median(Noise_SD)) );
%% reject things
EEG = pop_select( EEG,'nochannel',sort( unique([badEEGch, badEMGch, badNoiseCh, Other_chans]) ) );

%% redefine channels
All_chans = 1:EEG.nbchan;
getchantypes;

