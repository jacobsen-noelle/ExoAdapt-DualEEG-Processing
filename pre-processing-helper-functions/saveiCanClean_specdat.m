function cleanEEG = saveiCanClean_specdat(cleanEEG, rawEEG)

[EEG_spectra_in,FREQ,compeegspecdB_in,resvar_in,specstd_in] = spectopo(rawEEG.data(EEG_chans,timeWindow), length(timeWindow), rawEEG.srate, 'freqfac',1, 'plot','off');
            [EMG_spec_in,FREQ,compeegspecdB_in,resvar_in,specstd_in] = spectopo(rawEEG.data(EMG_chans,timeWindow), length(timeWindow), rawEEG.srate, 'freqfac',1, 'plot','off');
            [EMG_spec_out,FREQ,compeegspecdB_in,resvar_in,specstd_in] = spectopo(cleanEEG.data(EMG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off');
            [noise_spec_in] = spectopo(rawEEG.data(Noise_chans,timeWindow), length(timeWindow), rawEEG.srate, 'freqfac',1, 'plot','off');
            EEG_spectra_diff  = EEG_spectra_out-EEG_spectra_in; %difference between EEG spectra-in and spectr-out
            EMG_spectra_diff  = EMG_spec_out-EMG_spec_in; %difference between EMG spectra-in and spectr-out
            cleanEEG.etc.specdata.Noise_spectra_in = noise_spec_in;
            cleanEEG.etc.specdata.EMG_spectra_in = EMG_spec_in;
            cleanEEG.etc.specdata.EMG_spectra_out = EMG_spec_out;
            cleanEEG.etc.specdata.EEG_spectra_of_tdiff = EEG_spectra_of_tdiff;
            cleanEEG.etc.specdata.EEG_spectra_raw = EEG_spectra_in;
            cleanEEG.etc.specdata.EEG_spectra_clean = EEG_spectra_out;
            cleanEEG.etc.specdata.freq = FREQ;
            cleanEEG.etc.specdata.EEG_spectra_diff= EEG_spectra_diff;
            cleanEEG.etc.specdata.EMG_spectra_diff= EMG_spectra_diff;
end