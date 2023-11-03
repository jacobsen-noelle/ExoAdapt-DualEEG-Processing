EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
neckEMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
Fz_chans = find(strcmpi({EEG.chanlocs.type},'Fz') | strcmpi({EEG.chanlocs.type},'GRF'));
legEMG_chans =  find(strcmpi('leg-EMG',{EEG.chanlocs.type}));