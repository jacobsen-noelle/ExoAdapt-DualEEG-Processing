
% remove time intervals from EEG.dataset using code from clean_windows.m,  Copyright (C) Christian Kothe, SCCN, 2010, ckothe@ucsd.edu
% determine intervals to retain
function EEG = rm_timeintervals(EEG,sample_mask)
% apply selection
retain_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
% test = {};
% for i = 1:size(retain_data_intervals,1)
%     test{i,1} = [retain_data_intervals(i,1) retain_data_intervals(i,2)];
% end
%retain_data_intervals =retain_data_intervals./EEG.srate;
try
signal = EEG;
    signal = pop_select(signal, 'point', retain_data_intervals);
catch e
    if ~exist('pop_select','file')
        disp('Apparently you do not have EEGLAB''s pop_select() on the path.');
    else
        disp('Could not select time windows using EEGLAB''s pop_select(); details: ');
        hlp_handleerror(e,1);
    end
    %disp('Falling back to a basic substitute and dropping signal meta-data.');
    warning('Falling back to a basic substitute and dropping signal meta-data.');
    signal.data = signal.data(:,sample_mask);
    signal.pnts = size(signal.data,2);
    signal.xmax = signal.xmin + (signal.pnts-1)/signal.srate;    
    [signal.event,signal.urevent,signal.epoch,signal.icaact,signal.reject,signal.stats,signal.specdata,signal.specicaact] = deal(signal.event([]),signal.urevent([]),[],[],[],[],[],[]);
end

%update EEG dataset
signal.urevent = EEG.originaldata.event; %keep original event structure
signal.originaldata = EEG.originaldata;
signal.etc = EEG.etc;
EEG = signal;
end