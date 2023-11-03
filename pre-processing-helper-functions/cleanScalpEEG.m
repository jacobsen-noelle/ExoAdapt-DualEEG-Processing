% Perform pre-processing
% Cleanline on EMG and scalp to remove 60 HZ noise
% All other steps just done on scalp EEG
%  -remove bridged electrdoes
%  -unbiased rereferencing (temp remove bad channels)
%  -channel rejrection
%  -time window rejection
function EEG = cleanScalpEEG(EEG,outputFolder,cleanline_flag, check_chan_std_flag,chanlocsFolder)
%Setup directories
outputFolder2= strcat(outputFolder,'\\Step3a-cleanline');
outputFolder3= strcat(outputFolder,'\\scalp only\\cleanline');
outputFolder4= strcat(outputFolder,'\\scalp only\\cleanline-unbiased reref');
outputFolder5= strcat(outputFolder,'\\scalp only\\cleanline-unbiased reref-rejchan');
outputFolder6= strcat(outputFolder,'\\scalp only\\cleanline-unbiased reref-rejchan rejwind');
if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder);mkdir(outputFolder2); mkdir(outputFolder3);
    mkdir(outputFolder4); mkdir(outputFolder5); mkdir(outputFolder6)
end

getchantypes;
%check to make sure no channel type indicies are empty
if any(isempty(EEG_chans) || isempty(neckEMG_chans) || isempty(Noise_chans) || isempty(Fz_chans))
    error('Missing channel data'); return;
end
nonEEGchanidx = setdiff(1:EEG.nbchan,EEG_chans);
%% Add channel locations and labels
if isempty(EEG.chanlocs(1).theta)
    fprintf('File does not have channel locations, attempting to load locations...');
    EEG = chanlocs_add_masked(EEG,chanlocsFolder); %add channel locations, taking into account which were removed during processsing, loads chanloc file from subject folder
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG,EEG, CURRENTSET );
    %Optimze head center for electrode locations
    EEG=pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
    EEG.filename= strcat(extractBefore(EEG.filename,'.set'),'_chanlocs.set');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
end
og_filename = extractBefore(EEG.filename,'.set');

if cleanline_flag ==1
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[EEG_chans neckEMG_chans] ,'computepower',1,'linefreqs',60,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG.comments = pop_comments(EEG.comments,'','-cleanline 60Hz',1);
    EEG.setname = [EEG.subject,' Cleanline'];
    %EEG = eeg_checkset(EEG);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
   % EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_cleanline.set'),'filepath',outputFolder2);
    %eeglab redraw
end

% check standard deviation of all channels and visually inspect for outliers; save figure
if check_chan_std_flag ==1
    std_EEG = std(EEG.data([EEG_chans],:),0,2); %standard deviation of all channels
    std_EMG = std(EEG.data([neckEMG_chans],:),0,2);
    std_Noise = std(EEG.data([Noise_chans],:),0,2);
    stdDataBefore = [mean(std_EEG), mean(std_EMG), mean(std_Noise)]; %plot later
end
allchan=EEG;

%% Scalp electrode pre-processing
%extract EEG channels from dataset
getchantypes;%find index with scalp channels
EEG = pop_select( EEG,'channel',[EEG_chans]); %select scalp channels
EEG.setname='Scalp channels';
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
%EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_scalp'),'filepath',outputFolder3);

% check for electrode bridging
fprintf('\nChecking for electrode bridging')
[EB, ~] = eBridge(EEG);%bad channels should be in EB
bridgedchan = [];
%remove one of the electrodes in the bridged pair
for pairi = 1:size(EB.Bridged.Pairs,2)
    bridgedchan = [bridgedchan, find(strcmpi(EB.Bridged.Pairs{2,pairi},{EEG.chanlocs.labels}))];
end
bridgedchan = unique(bridgedchan);
if size(bridgedchan,2)> 16 % if there are too many bridged electrodes ,...
    % we can't use this subject. 16 is abritrary cuttoff--...
    % really depends on the distribution of the bridged electrdoes.
    warning('Too many bridged electrodes; Skipping this subject');
    return;
end
if ~isempty(bridgedchan)
    EEG = pop_select( EEG,'nochannel',bridgedchan); %select scalp channels
    EEG.etc.eBridge = EB;
else
    fprintf('No channel bridging was found')
end

%% Robust unbiased rereference- EEG only
% Identify bad channels and calculate average voltage across all
% "good" channels. Subtract mean from ALL channels. Then in next step
% do actual channel rejection
originalEEG = EEG; % Keep original EEG.
% identify bad channels and exclude their values for unbiased rereference
fprintf('\nComputing initial channel rejection for robust referencing\n')
[unbiasedEEG] = channelrejection(originalEEG,'WindowCriterion', 'off');  %Auto channel rejection
fprintf('\nChannel rejection for unbiased reref \n\t-pop_rejchan \n\t\tkurtosis thres=5 , prob thres = 5 \n\t-trimOutlier \n\t\tstd thres = 500')
fprintf('\n\t-clean_rawdata \n\t\t arg_burst =-1; arg_highpass = [0.25 0.75]; arg_window = -1; arg_channel= 0.8; arg_noisy = 4; arg_burstRefMaxBadCh = 0.75;');

%find differences in channels kept from original
for index1 = 1:originalEEG.nbchan
    a = originalEEG.chanlocs(index1).labels;
    if isempty(find(strcmpi(a,{unbiasedEEG.chanlocs.labels})))
        fprintf('%s was excluded from unbiased rereference\n',a);
    end
end
unbiased_avg = mean(unbiasedEEG.data,1);
reref_EEG = originalEEG;
reref_EEG.data = originalEEG.data-unbiased_avg;

%% Rereference using "good" channels only
if useGPU ==1
    g = gpuArray(unbiasedEEG.data);
    unbiased_avg = mean(g,1);
    reref_EEG = originalEEG;
    reref_EEG.data = originalEEG.data-unbiased_avg;
    reref_EEG.data = gather(reref_EEG.data); %transfer gpuArray to local workspace
else
    unbiased_avg = mean(unbiasedEEG.data,1);
    reref_EEG = originalEEG;
    reref_EEG.data = originalEEG.data-unbiased_avg;
end
reref_EEG.ref = 'common';
[reref_EEG.chanlocs.ref] = deal('common');
reref_EEG.setname='Scalp channels- unbiased avg';
reref_EEG.comments = pop_comments(EEG.comments,'','Average re-referenced with the average taken after channel rejection (unbiased average). Average only subtracted from cortical channels',1);
reref_EEG = eeg_checkset( reref_EEG );
fprintf('\nUnbiased average re-reference complete\n');

%% Reject bad EEG channels
fprintf('Starting channel rejection using clean_rawdata...\n')
nchan = reref_EEG.nbchan;
[reref_chanrej_EEG, reref_chanrej_timerej] = channelrejection(reref_EEG);
reref_chanrej_timerej.comments = pop_comments(EEG.comments,'',...
    'reject bad channels using at least two of the following methods: kurtosis, probability, standard deviation, clean_artifacts \n time window rej',1);
fprintf('%i scalp channel(s) removed\n',originalEEG.nbchan-reref_chanrej_EEG.nbchan);

%identify which channels were removed
for index1 = 1:length(originalEEG.chanlocs)
    a = originalEEG.chanlocs(index1).labels;
    if isempty(find(strcmpi(a,{reref_chanrej_EEG.chanlocs.labels}), 1))
        fprintf('%s was removed \n',a);
    end
end

% Determine number of time samples rejected
frames_before= length(reref_EEG.times);
frames_after=length(reref_chanrej_timerej.times);
p_frames_kept = (frames_after/frames_before)*100;
p_frames_rej = 100-p_frames_kept;
disp('Channel and time window rejection complete');

fprintf('\nChannel cleaning parameters');
for i= 1:length(EEG.etc.clean_artifacts.parameters)
    fprintf('\n\t%s: %s',EEG.etc.clean_artifacts.parameters(i).crit,num2str(EEG.etc.clean_artifacts.parameters(i).value));
end

%% Add back other channels
EEGOUT= reref_chanrej_EEG; %copy EEG structure, then subsequently replace
EEGOUT.originaldata = allchan;
EEGOUT.data =[]; %remove data so we can add combined time
EEGOUT.pnts = reref_chanrej_EEG.pnts;

if p_frames_rej~=0
    clean_sample_mask = reref_chanrej_timerej.etc.clean_artifacts.clean_sample_mask; % 1=keep frame, 0= rej frame
    EEGOUT.etc.percent_frames_rej= p_frames_rej;

elseif p_frames_rej ==0 %if no time windows were rejected with clean_artifacts(), just use trial trim mask
    clean_sample_mask  =true(1,size(EEGOUT.data,2));% 1 = keep frame
end
EEGOUT.data = [reref_chanrej_EEG.data(:,:); allchan.data(nonEEGchanidx,:)];%add all non-EEG channels back to EEG dataset

[EEGOUT.chanlocs(end+1:end+length(nonEEGchanidx))] = allchan.chanlocs(nonEEGchanidx);
EEGOUT.nbchan = size(EEGOUT.data,1);

%% Reject time window
if size(clean_sample_mask,2) == size(EEGOUT.data,2)
    EEGOUT.urevent = allchan.event;
    EEGOUT = rm_timeintervals(EEGOUT,clean_sample_mask); % func from cleanrawdata
    EEGOUT.etc.clean_artifacts = reref_chanrej_timerej.etc.clean_artifacts;
    if size(find(clean_sample_mask ==1),2) ~= size(EEGOUT.data,2)
        error('Size of good samples in mask and output EEG.data do not match'); return;
    end
else
    error('Size of time window rejection mask and EEG.data do not match');return;
end
%%
%make sure no channel types were lost
[EEG_chans, neckEMG_chans, Noise_chans, Fz_chans] = chantypes(EEGOUT);
if isempty(Noise_chans) || isempty(neckEMG_chans) || isempty(EEG_chans) || isempty(Fz_chans)
    errors('Missing one of channel types (EEG,EMG,noise)');
end

%one last reref with all good channels
EEGOUT = rerefC2CN2NExt2Ext_func(EEGOUT,1);

%store
EEGOUT.setname = [EEGOUT.subject,' fully pre-processed'];
EEGOUT.comments = pop_comments(EEGOUT.comments,'', 're-reference cort to cort, noise to noise, and ext to ext',1);
EEGOUT.comments = pop_comments(EEGOUT.comments,'', ['fullRankBool=',num2str(fullRankBool)],1);
%[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG,EEG, CURRENTSET);
EEGOUT = eeg_checkset(EEGOUT);

% check standard deviation of all channels and visually inspect for outliers; save figure
if check_chan_std_flag ==1
    stdData = std(EEGOUT.data(EEG_chans,:),0,2); %check standard deviation of channels
    fig= figure; bar(stdData)
    hold on;
    yline(std(stdData,0,1),'--','1 SD, all chan');
    yline(3*std(stdData,0,1),'--','3 SD, all chan');
    th = title('Cleaned Scalp Channels');  th.FontSize = 18;
    yl = ylabel("Standard Deviation"); yl.FontSize = 16;
    xticks(1:EEGOUT.nbchan);
    set(gca, 'XTickLabel',{EEGOUT.chanlocs(EEG_chans).labels},'FontSize',10);
    savethisfig(fig,strcat(EEGOUT.subject,'_EEGchan_std'), strcat(outputFolder,'\Figures'),'jpg');
    close;
end

end

function [EEG_chans, neckEMG_chans, Noise_chans, Fz_chans] = chantypes(EEG)
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
Noise_chans = find(strcmpi('Noise',{EEG.chanlocs.type}));
neckEMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
Fz_chans = find(strcmpi({EEG.chanlocs.type},'Fz') | strcmpi({EEG.chanlocs.type},'GRF'));
end