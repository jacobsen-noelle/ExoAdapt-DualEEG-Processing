% Authoer: Ryan Downey, University of Florida
function EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankBool)
%re-reference cort to cort, noise to noise, and 

if fullRankBool
    disp('Rereferencing cort to cort, noise to noise, and emg to emg (with ZERO loss of rank!)');
else
    disp('Rereferencing cort to cort, noise to noise, and emg to emg (with a loss of rank of 1 for each rereferencing)');
end

%% define channels
% EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
% Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
% EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type})); 
%note: may need to add more lines in future if other external sensors are
%used with a different channel type than 'EMG' (e.g. could use an array of
%EOG sensors and want to re-reference them to themselves)

getchantypes; %NJ
%% re-ref cort to cort and noise to noise and emg to emg 
% EEG = pop_reref( EEG, [],'exclude',Noise_chans );
% EEG = pop_reref( EEG, [],'exclude',sort([EEG_chans, EMG_chans]) );

chansin = Noise_chans;
        nchansin = length(chansin);
        if fullRankBool
            refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%new full rank method 
        else
            refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method)
        end
        chansout = chansin;
        EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
chansin = EEG_chans;
        nchansin = length(chansin);
        if fullRankBool
            refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%new full rank method 
        else
            refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method)
        end
        chansout = chansin;
        EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
chansin = neckEMG_chans; %NJ
        nchansin = length(chansin);
        if fullRankBool
            refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%new full rank method 
        else
            refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method)
        end
        chansout = chansin;
        EEG.data(chansout,:) = refmatrix*EEG.data(chansin,:);
        
end