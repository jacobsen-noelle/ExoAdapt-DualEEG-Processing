% run_iCanClean ()
% Use dual-electrode EEG and neck EMG channels to attenuate motion and
% muscle artifacts in scalp EEG channels. Uses canonical correlation
% analnum_cycles sis to find mutal components between two sets of channels
%
% Inputs:
%   EEG                     EEGLAB struct
%   visualizeFinalResults   [0|1] plot data with and without artifacts removed
%   iCCmethod               [1|2|3] iCanClean Method
%                           1) clean EEG w noise channels
%                           2) clean EEG w EMG channels
%                           3) clean EEG w noise then EMG channels
%   mainOutFolder           Main outputfolder for datasets, notes, spec
%                           subfolder 1-motion artifact remover, subfolder 2-motion and EMG
% Outputs:
%   cleanEEG        EEG with motion and muscle components removed from
%                   EEG.data
%

%
% Author: Noelle Jacobsen, Universitnum_cycles  of Florida
% Created: 2023-Sept-14
% Last updated: 2023-Sept-14

%%MAKE SURE FORCEPLTE DATA IS RETAINED
%input test name
% n = 0;
% for x = 1:length(testwindowLengths)
%     for xx=1:length(testRhoSqThres_source)
%         n=n+1;
%         testparams(n).rhoSqThres_source = testRhoSqThres_source(xx);
%         testparams(n).windowLength = testwindowLengths(x);
%     end
% end
%  %for testnum= 1:length(testparams)




function cleanEEG = run_iCanClean(EEG,ALLEEG,CURRENTSET,mainOutFolder, opt)
%load params
tStart=tic;
testname = opt.testname;
iCCmethod = opt.iCCmethod;
calc_psd = opt.calc_psd; %to comparem PSD before and after iCC
visualizeFinalResults = opt.visualizeFinalResults;
removeChan_flag = opt.removeChannels;
% Channel cleaning methods
% 1) clean EEG w noise channels
% 2) clean EEG w EMG channels
% 3) clean EEG w noise then EMG channels
%% setup file and notes export
%outputfolder for iCanClean output data (main output for notes, spec
%data, subfolder 1-motion artifact remover, subfolder 2-motion and EMG
%removed
subFolder1= [mainOutFolder,'\\RemoveMotionArtifact'];
subFolder2= [mainOutFolder,'\\PSD'];

if ~exist(subFolder2, 'dir') %check if dir exists
    mkdir(subFolder2)
end
if ~exist(subFolder1, 'dir') %check if dir exists
    mkdir(subFolder1)
end
cd(mainOutFolder)
% if isfile('SPECTDATA.mat') %big file with all spec data
%     load('SPECDATA.mat')
% end
%Text file to write data processing notes to
notesFile = 'Notes.txt';
if isfile(notesFile) %check if file exists
    notes = fopen(notesFile,'a+'); % if it exists, append
else
    notes = fopen(notesFile,'w');
    fprintf(notes,'%s\n',subFolder2);
    fprintf(notes,'%s\n',date);
end

fprintf(notes,[EEG.filename,'\n']);
fprintf(notes,'1st Pass- motion removed from EEG & EMG, moving window\n');
fprintf(notes,'2nd Pass- motion removed from EEG & EMG, inf window\n');
fprintf(notes,'3rd Pass- muscle removed from EEG, moving window\n');


% Add channel types if missing
for i= 1:length(EEG.chanlocs)
    if contains({EEG.chanlocs(i).labels}, 'EMG')
        EEG.chanlocs(i).type = 'EMG';
    elseif contains({EEG.chanlocs(i).labels}, 'EXG')
        EEG.chanlocs(i).type = 'EMG';
    elseif contains({ EEG.chanlocs(i).labels},'2-')
        EEG.chanlocs(i).type = 'Noise';
    elseif contains({ EEG.chanlocs(i).labels},'1-') && ~contains({EEG.chanlocs(i).labels}, 'EXG')
        EEG.chanlocs(i).type = 'EEG';
    elseif contains({EEG.chanlocs(i).labels},'GRF')
        EEG.chanlcos(i).type = 'GRF';
    else
        EEG.chanlcos(i).type = 'Exo';
    end
end
%eeglab redraw;
%     %% Add channel types
%     for i= 1:size(EEG.chanlocs,1)
%         if contains({EEG.chanlocs(i).labels}, 'Fz')
%             EEG.chanlocs(i).type = 'Fz';
%         elseif contains({EEG.chanlocs(i).labels}, 'EMG')
%             EEG.chanlocs(i).type = 'EMG';
%         elseif contains({EEG.chanlocs(i).labels}, 'EXG')
%             EEG.chanlocs(i).type = 'EMG';
%         elseif contains({ EEG.chanlocs(i).labels},'2-')
%             EEG.chanlocs(i).type = 'Noise';
%         elseif contains({ EEG.chanlocs(i).labels},'1-') && ~contains({EEG.chanlocs(i).labels}, 'EXG')
%             EEG.chanlocs(i).type = 'EEG';
%         end
%     end
%% Add chan locs if not present
if isempty(EEG.chanlocs(1).theta)
    fprintf('File does not have channel locations, attempting to load locations...');
    tempEEG = chanlocs_add_masked(EEG,chanlocsFolder); %add channel locations, taking into account which were removed during processsing, loads chanloc file from subject folder
    getchantypes;
    oldlocs = EEG.chanlocs;
    EEG.chanlocs = tempEEG.chanlocs;
    [EEG.chanlocs([neckEMG_chans Noise_chans]).labels] = oldlocs([neckEMG_chans Noise_chans]).labels;
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end
%% Find channel index
All_chans = [1:EEG.nbchan];

getchantypes;
back_head = [10 11 12 13 14 15 23 24 25 26 27 28 39 40 41];
top_head = [97 98 1 2 65 66 33 34];
front_head = [ 91 92 93 81 82 83 84 78 79 80];
rawEEG=EEG;
chanidx = [EEG_chans, neckEMG_chans, Noise_chans];
otherchanidx = setdiff([1:EEG.nbchan],chanidx);
%% Quick pre-processing
%HPF Filter all data -1 Hz
fc=1;% cut off frequency (Hz)
fn=EEG.srate/2; %nyquivst frequency = sample frequency/2;
order = 2; %2nd order filter, high pass
[b, a]=butter(order,(fc/fn),'high');
dataOut = filtfilt(b,a,double(EEG.data(chanidx,:))');
EEG.data(chanidx,:) = dataOut';
EEG.comments = pop_comments(EEG.comments,'',strcat('-HPF ','', num2str(fc),'','Hz, ',num2str(order),'-order Butterworth, zero-phase'),1);
fprintf(notes,'1Hz HPF, 2nd Order Butterworth, zero-phase\n');


%Conservative channel rejection
% remove channels STD>3
fullRankBool = 1;
EEG_quickreref = rerefC2CN2NExt2Ext_func(EEG,fullRankBool);% quick reference before channel rejection, but don't ultimately use this EEG data

%quick conservative channel rejection before iCanClean
[~, EEG_chans, Noise_chans, neckEMG_chans, Other_chans, ~, badEEGch, badEMGch, badNoiseCh,] = autoChRej(EEG_quickreref);
if ~isempty(badEEGch) || ~isempty(badEMGch) || ~isempty(badNoiseCh)
    fprintf('\nRemoving channels >3SD: ');
    if ~isempty(badEEGch)
        fprintf('\n\t-Bad EEG channel(s): [%s]\n',join(string(badEEGch), ','))

    end
    if ~isempty(badEMGch)
        fprintf('\t-Bad EMG channel(s): [%s]\n',join(string(badEMGch), ','))
    end
    if ~isempty(badNoiseCh)
        fprintf('\t-Bad noise channel(s): [%s]\n',join(string(badNoiseCh), ','))
    end
end

% fprintf('Interpolating bad EEG channels')
% EEG_scalp = pop_select( EEG,'nochannel',sort( unique([badEEGch, badEMGch, badNoiseCh, Other_chans]) ) ); %remove bad channels from data before rereference so only rereferencing once
% EEG = pop_interp(EEG, badEEGch, method);

fprintf('Temporarily removing extra channels\n')
EEG = pop_select( EEG,'nochannel',sort( unique([badEEGch, badEMGch, badNoiseCh, Other_chans]) ) ); %remove bad channels from data before rereference so only rereferencing once
EEG = rerefC2CN2NExt2Ext_func(EEG,fullRankBool); %rereference again after channel removal
getchantypes;
filtEEG= EEG;
if isempty(EEG_chans)
    error('No EEG chans selected')
end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%% Calculate PSD of input data
if calc_psd ==1
    %check to see if EEG spectra already been calculated
    timeWindow = find(EEG.times/1000 >= 0 & EEG.times/1000 <= inf); %arbitrarily calculating PSD no smaller window so we don't have to wait forever
    if ~isfield(EEG,'etc.specdata.EEG_spectra_raw')
        [EEG_spectra_in,FREQ,~,~,~] = spectopo(filtEEG.data(EEG_chans,timeWindow), length(timeWindow), filtEEG.srate, 'freqfac',1, 'plot','off','percent',10);
        specdata.EEG_spectra_raw = EEG_spectra_in;
        data_change =1;
    else
        EEG_spectra_in = specdata.EEG_spectra_raw;
    end
    %check to see if noise spectra already been calculated
    if ~isfield(EEG,'etc.specdata.Noise_spectra_in')
        [noise_spec_in] = spectopo(filtEEG.data(Noise_chans,timeWindow), length(timeWindow), filtEEG.srate, 'freqfac',1, 'plot','off','percent',10);
        EEG.etc.specdata.Noise_spectra_in = noise_spec_in;
        data_change =1;
    else
        noise_spec_in = specdata.Noise_spectra_in;
    end
    %check to see if EMG spectra already been calculated
    if ~isfield(EEG,'etc.specdata.EMG_spectra_in')
        [EMG_spec_in,FREQ,~,~,~] = spectopo(filtEEG.data(neckEMG_chans,timeWindow), length(timeWindow), filtEEG.srate, 'freqfac',1, 'plot','off','percent',10);
        EEG.etc.specdata.EMG_spectra_in = EMG_spec_in;
        data_change =1;
    else
        EMG_spec_in = specdata.EMG_spectra_in;
    end
    %store info
    if data_change ==1
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        %         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',EEG.subject,'overwrite','on','gui','off');
    end
    iCC_specdata.(testname).HPFdata = specdata; %save raw data for meta analysis
end
%% Run iCanClean
if iCCmethod ==1 % default cleaning (remove motion artifact from EEG channels)
    cleanEEG = iCanClean( EEG, [EEG_chans], [Noise_chans], visualizeFinalResults, params);
    %vis_artifacts(cleanEEG,filtEEG,'ChannelSubset',[EEG_chans(1:2:end) Noise_chans(1:4:end)]);
    %vis_artifacts(test5,test1,'ChannelSubset',[EEG_chans(1:2:end)]);
elseif iCCmethod ==2  % one pass removal of motion artifacts and neck muscle from EEG
    cleanEEG = iCanClean( EEG, [EEG_chans], [Noise_chans, neckEMG_chans], visualizeFinalResults, params);
elseif iCCmethod ==3
    % two or three step removal of motion artifacts and neck muscle from EEG
    params = [];
    pass=1;
   % params(pass).windowLength = 1; %length of data that gets cleaned (set to inf for entire dataset, set to some number for moving window)
    params(pass).cleanWindow = 1;
    params(pass).rhoSqThres_source = 0.95;
    params(pass).cleanYBool = false;
    params(pass).giveCleaningUpdates=0;
    params(pass).X = '[EEG_chans, neckEMG_chans]';
    params(pass).Y = '[Noise_chans]';
    params(pass).useBoundary =1;
    params(pass).plotStatsOn= 0;
    %         params.extraTime_pre = 0.95;  %... adding extra data on front and back
    %         params.extraTime_post = 0.95;
    %         params.stepSize = 0.1;
    %first pass - remove motion artifact, moving window
    fprintf('\niCC run #1')
    cleanEEG1 = iCanClean( EEG, [EEG_chans, neckEMG_chans], [Noise_chans], visualizeFinalResults, params(pass));
    %iCC_specdata.(testname)(f).firstpass.params = params; %keep track of params used
    if visualizeFinalResults
        savethisfig(gcf,[extractBefore(EEG.filename,'_chanlocs'),'_pass1'],[subFolder1,'\vis_artifacts'],'fig');
    end
    output(pass).numNoiseCompsRemovedPerWindow = cleanEEG1.etc.iCanClean.numNoiseCompsRemovedPerWindow;
    output(pass).numNoiseCompsRemovedOnAvg = cleanEEG1.etc.iCanClean.numNoiseCompsRemovedOnAvg;
    output(pass).params = params(pass);
   
    %second pass- remove motion artifact, inf window
    pass=2;
    params(pass).cleanWindow = inf;
    %params(pass).windowLength = inf; %length of data that gets cleaned (set to inf for entire dataset, set to some number for moving window)
    params(pass).rhoSqThres_source = 0.05;
    params(pass).rhoSqThres_source = 0.05;
    params(pass).cleanYBool = false;
    params(pass).giveCleaningUpdates=0;
    params(pass).X = '[EEG_chans, neckEMG_chans]';
    params(pass).Y = '[Noise_chans]';
    params(pass).plotStatsOn= 0;
    params(pass).useBoundary = 1;
    fprintf('\niCC run #2')
    cleanEEG2 = iCanClean(cleanEEG1, [EEG_chans, neckEMG_chans], [Noise_chans], visualizeFinalResults, params(pass));
    cleanEEG2.iCanClean.params = params;%concantenate all parameters used
    %iCC_specdata.(testname)(f).secondpass.params = params;%keep track of params used
     if visualizeFinalResults
        savethisfig(gcf,[extractBefore(EEG.filename,'_chanlocs'),'_pass2'],[subFolder1,'\vis_artifacts'],'fig');
     end
    output(pass).numNoiseCompsRemovedPerWindow = cleanEEG2.etc.iCanClean.numNoiseCompsRemovedPerWindow;
    output(pass).numNoiseCompsRemovedOnAvg = cleanEEG2.etc.iCanClean.numNoiseCompsRemovedOnAvg;
    output(pass).params = params(pass);
  

    %third pass- remove muscle artifacts,moving window
    pass = 3;
    fc=20;% cut off frequency (Hz);
    order = 2; %2nd order filter, high pass
    [b, a]=butter(order,(fc/fn),'high');
    dataOut = filtfilt(b,a,double(cleanEEG2.data(neckEMG_chans,:))'); %filter neck EMG so iCC can focus on muscle activity
    cleanEEG2_filt = cleanEEG2;
    cleanEEG2_filt.data(neckEMG_chans,:) = dataOut';

    %params(pass).windowLength = 1; %length of data that gets cleaned (set to inf for entire dataset, set to some number for moving window)
    params(pass).cleanWindow = 1; %length of data that gets cleaned (set to inf for entire dataset, set to some number for moving window)
    params(pass).rhoSqThres_source = 0.5; %MiM used 0.25
    params(pass).cleanYBool = false;
    params(pass).giveCleaningUpdates=0;
    params(pass).X = '[EEG_chans]';
    params(pass).Y = '[neckEMG_chans]';
    params(pass).plotStatsOn= 0;
    params(pass).useBoundary = 1;
    fprintf('\niCC run #3')
    cleanEEG3 = iCanClean(cleanEEG2_filt, EEG_chans, neckEMG_chans, visualizeFinalResults, params(pass));
     if visualizeFinalResults
        savethisfig(gcf,[extractBefore(EEG.filename,'_chanlocs'),'_pass3'],[mainOutFolder,'\vis_artifacts'],'fig');
     end
    output(pass).numNoiseCompsRemovedPerWindow = cleanEEG3.etc.iCanClean.numNoiseCompsRemovedPerWindow;
    output(pass).numNoiseCompsRemovedOnAvg = cleanEEG3.etc.iCanClean.numNoiseCompsRemovedOnAvg;
    output(pass).params = params(pass);

end

%iCC_specdata.(testname)(f).thirdpass.params = params;%keep track of params used
%close;
%% visualize time series results
%vis_artifacts(cleanEEG3,EEG,'ChannelSubset',[EEG_chans(back_head) neckEMG_chans(1:2:end)]);

%         %vis_artifacts(cleanEEG,filtEEG,'ChannelSubset',[10 11 12 14 23 24 25 27 39 41 119   128   138   148   158   168   178   188   198   208 ]);
%         %vis_artifacts(cleanEEG,EEG,'ChannelSubset',[EEG_chans(front_head) neckEMG_chans]);
%         %vis_artifacts(cleanEEG,filtEEG,'ChannelSubset',[EEG_chans(1:2:end) neckEMG_chans Ychans(1:4:end)]);
%         %vis_artifacts(cleanEEG,EEG,'ChannelSubset',[Xchans(1:2:end) Ychans]);
%
%% store EEG dataset
%store  data from second pass
[ALLEEG, cleanEEG2, CURRENTSET] = eeg_store( ALLEEG, cleanEEG2, CURRENTSET );
cleanEEG2 = eeg_checkset( cleanEEG2 );
cleanEEG2.setname = strcat(cleanEEG2.subject,' iCC');
if ~isempty(otherchanidx)  %keep external channels
    cleanEEG2.data = [cleanEEG2.data; rawEEG.data(otherchanidx,:)];
    cleanEEG2.nbchan = size(cleanEEG2.data,1); %update # of chans
    numNewSignals = length(otherchanidx);
    %% Update channel labels and type - Forceplates,exo, leg EMG
    for newsig_i = 1:numNewSignals
        cleanEEG2.chanlocs(end+1).labels = rawEEG.chanlocs(otherchanidx(newsig_i)).labels;
        cleanEEG2.chanlocs(end).type = rawEEG.chanlocs(otherchanidx(newsig_i)).type;
    end
end
[ALLEEG, cleanEEG2, CURRENTSET] = eeg_store( ALLEEG, cleanEEG2, CURRENTSET );
cleanEEG2 = pop_saveset(EEG, 'filename',strcat(extractBefore(cleanEEG2.filename,'.set'),'_iCC'),'filepath',subFolder1,'version', '7.3');

[ALLEEG,cleanEEG3, CURRENTSET] = eeg_store( ALLEEG, cleanEEG3, CURRENTSET );

%% visualize spectra results
if calc_psd ==1
    %calculate spectra of cleaned EEG data
    %from 2nd pass (motion artifacts removed)
    cleanEEG=cleanEEG2;
    timeWindow = find(cleanEEG.times/1000 >= 0 & cleanEEG.times/1000 <= inf); %arbitrarily calculating PSD no smaller window so we don't have to wait forever
    [EEG_spectra_out,FREQ,compeegspecdB_out,resvar_out,specstd_out] = spectopo(cleanEEG.data(EEG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off','percent',10);
    %[EEG_spectra_of_tdiff] =  spectopo( cleanEEG.data(EEG_chans,timeWindow)-filtEEG.data(EEG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off');

    %calulate other spectras
    if iCCmethod ==1
        [EMG_spec_out,FREQ,~,~,~] = spectopo(cleanEEG.data(neckEMG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off','percent',10);
        EEG_spectra_diff  = EEG_spectra_out-EEG_spectra_in; %difference between EEG spectra-in and spectr-out
        EMG_spectra_diff  = EMG_spec_out-EMG_spec_in; %difference between EMG spectra-in and spectr-out
        specdata.Noise_spectra_in = noise_spec_in;
        specdata.EMG_spectra_in = EMG_spec_in;
        specdata.EMG_spectra_out = EMG_spec_out;
        specdata.EEG_spectra_of_tdiff = EEG_spectra_of_tdiff;
        specdata.EEG_spectra_raw = EEG_spectra_in;
        specdata.EEG_spectra_clean = EEG_spectra_out;
        specdata.freq = FREQ;
        specdata.EEG_spectra_diff= EEG_spectra_diff;
        specdata.EMG_spectra_diff= EMG_spectra_diff;
        iCC_specdata.(testname).firstpass.specdat = specdata;
        iCC_specdata.(testname).subject = filtEEG.subject;
        iCC_specdata.(testname).filename = filtEEG.filename;
    elseif iCCmethod ==3
        [EMG_spec_out,FREQ,~,~,~] = spectopo(cleanEEG.data(neckEMG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off','percent',10);
        EMG_spectra_diff  = EMG_spec_out-EMG_spec_in; %difference between EMG spectra-in and spectr-out
        EEG_spectra_diff  = EEG_spectra_out-EEG_spectra_in;
        %cleanEEG.etc.specdata.spectra_of_tdiff = EEG_spectra_of_tdiff;
        specdata.spectra_clean = EEG_spectra_out;
        specdata.freq = FREQ;
        specdata.spectra_diff= EEG_spectra_diff;
        specdata.EMG_spectra_out =EMG_spec_out;
        iCC_specdata.(testname).secondpass.specdat = specdata;
        iCC_specdata.(testname).subject = filtEEG.subject;
        iCC_specdata.(testname).filename = filtEEG.filename;
    end

    %         % rnk = rank(double(cleanEEG.data)); % Calculate rank of output data (component removal reduces rank)
    %         % fprintf(['CCA EEG rank out: ' int2str(rnk) '\n']
    disp('Spectras have been calculated')

    %% Plot PSD
    %
    %plotting colors
    blue = [0 0.4470 0.7410];
    orange = [0.8500, 0.3250, 0.0980];
    purple = [0.4940, 0.1840, 0.5560];
    green =  [0.4660, 0.6740, 0.1880];
    cyan = [0.3010, 0.7450, 0.9330];
    red = [0.6350, 0.0780, 0.1840];
    yellow = [0.9290, 0.6940, 0.1250];

    spec_plot= figure;%('WindowState','maximized'); % Plot spectra of mean channel input and output
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 8 6]);
    set(gcf,'Color','w');
    mytitle= title({'Median Spectral Power Density',...
        strcat(extractBefore(filtEEG.filename,'.set'),' ',testname,' ||  Noise Removed')},'Interpreter', 'none');
    hold on;
    plot(FREQ, median(EEG_spectra_in,1),'Color',cyan,'LineWidth', 2); plot(FREQ, prctile(EEG_spectra_in,25), '--','Color',cyan,'LineWidth', .5); plot(FREQ, prctile(EEG_spectra_in,75), '--','Color',cyan,'LineWidth', .5);
    set(gca,'FontSize',12);
    plot(FREQ, median(EEG_spectra_out,1),'Color',green,'LineWidth', 2); plot(FREQ, prctile(EEG_spectra_out,25), '--','Color',green,'LineWidth', .5); plot(FREQ, prctile(EEG_spectra_out,75), '--','Color',green,'LineWidth', .5);
    plot(FREQ, median(EEG_spectra_diff,1), 'k','LineWidth',1);
    legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Difference between EEG IN vs OUT- Q2','Location', 'North', 'Orientation', 'horizontal');

    if iCCmethod ==1
        plot(FREQ, median(noise_spec_in,1) ,'Color',red,'LineWidth',1);
        legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Difference between EEG IN vs OUT- Q2','Noise Spectra','Location', 'North', 'Orientation', 'horizontal');
    elseif iCCmethod==2
        plot(FREQ, median(noise_spec_in,1),'Color',red,'LineWidth',1);
        plot(FREQ, median(EMG_spec_in,1) ,'m','LineWidth',1);
        plot(FREQ, median(EMG_spec_out,1) ,'Color',purple,'LineWidth',1);
        legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Difference between EEG IN vs OUT- Q2','Raw Noise Spectra','EMG IN','EMG OUT/Clean','Location', 'North', 'Orientation', 'horizontal');
    elseif iCCmethod==3
        plot(FREQ, median(noise_spec_in,1) ,'Color',red,'LineWidth',1);
        plot(FREQ, median(EMG_spec_in,1) ,'Color','m','LineWidth',1);
        plot(FREQ, median(EMG_spec_out,1) ,'Color',purple,'LineWidth',1);
        legend('EEG IN, Q2', 'Q1','Q3', 'EEG OUT/CLEAN, Q2','Q1','Q3','Difference between EEG IN vs OUT- Q2','Raw Noise Spectra','EMG IN','EMG OUT/Clean','Location', 'North', 'Orientation', 'horizontal');
    end

    %             legend boxoff;
    xlabel('Frequency (Hz)');
    ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
    hold off;
    xlim([.5 100]);  set(gca,'Box','off');
    %YLIM = get(gca,'YLim');
    ylim([-30 30])
    set(gca,'XTick', [0.5,4,8,13,30,60,100]);
    set(gca, 'XTickLabel', {'','4','8','13','30','60','100'});

    %save figure to

    savethisfig(spec_plot,strcat(extractBefore(filtEEG.filename,'.set'),' ',testname,'-',mytitle.String{1}), subFolder1,'jpg');
    close;
    %______________________________________________
    %     %% Plot Mean Spectra
    %     spectra_in = EEG_spectra_in;
    %     spectra_out = EEG_spectra_out;
    %     spectra_diff = EEG_spectra_of_tdiff;
    %     spec_plot= figure; % Plot spectra of mean channel input and output
    %     set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 12 6]);
    %     set(gcf,'Color','w');
    %     mytitle=title(strcat(testname,' Mean Power spectra'));
    %     hold on;
    %     plot(FREQ, mean(spectra_in,1),'Color',cyan,'LineWidth',2); plot(FREQ, prctile(spectra_in,25), '--','Color', cyan,'LineWidth', .5); plot(FREQ, prctile(spectra_in,75), '--','Color',cyan,'LineWidth', .5);
    %     set(gca,'FontSize',12);
    %     plot(FREQ, mean(spectra_out,1), 'Color', green,'LineWidth', 2); plot(FREQ, prctile(spectra_out,25), '--','Color', green,'LineWidth', .5); plot(FREQ, prctile(spectra_out,75), '--','Color', green,'LineWidth', .5);
    %     plot(FREQ, mean(spectra_diff,1), 'k','LineWidth',1);
    %     legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Diff Q2','Location', 'North', 'Orientation', 'horizontal');
    %
    %     if Y ==1
    %         plot(FREQ, mean(noise_spec_in,1) ,'Color',red,'LineWidth',1);
    %         legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Diff Q2','Noise Spectra','Location', 'North', 'Orientation', 'horizontal');
    %     elseif Y==2 || Y==3
    %         plot(FREQ, mean(EMG_spec_in,1) ,'Color',purple,'LineWidth',1);
    %         plot(FREQ, mean(noise_spec_in,1) ,'Color',red,'LineWidth',1);
    %         legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Diff Q2','EMG Spectra','Noise Spectra','Location', 'North', 'Orientation', 'horizontal');
    %     end
    %     legend boxoff;
    %     xlabel('Frequency (Hz)');ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
    %     hold off;
    %     xlim([.5 70]);  set(gca,'Box','off');
    %     YLIM = get(gca,'YLim');
    %     set(gca,'XTick', [0.5,4,8,13,30,60],'XTickLabel', {'','4','8','13','30','60'});
    %
    %     %save figure
    %     addpath('C:\Users\jacobsen.noelle\Documents\GitHub\EEG_Processing\helper-functions');
    %     savethisfig(spec_plot,mytitle.String, strcat(inputFolder,'\Figures'));
    %     close;
    %
    %     %plot EEG spec before/after iCanClean w/ shaded error bars
    %     spec_plot=figure;
    %     set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 12 6]);
    %     shadedErrorBar(FREQ, EEG_spectra_in, {@mean,@std}, 'lineprops', '-')
    %     xlim([0 60])
    %     hold on;
    %     shadedErrorBar(FREQ, EEG_spectra_out, {@mean,@std}, 'lineprops', '-')
    %     xlabel('Frequency (Hz)');
    %     ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
    %     mytitle=title('EEG Power Spectra');
    %     legend('Before CCA cleaning','After CCA cleaning')
    %     set(gcf,'Color','w','FontSize',20,'XTick', [0.5,4,8,13,30,60],'XTickLabel', {'','4','8','13','30','60'})
    %
    %     %save figure to
    %     addpath('C:\Users\jacobsen.noelle\Documents\GitHub\EEG_Processing\helper-functions');
    %     savethisfig(spec_plot,mytitle.String, strcat(inputFolder,'\Figures'));
    %     close;
    %     %% store the cleaned data
    %     if Y ==1
    %         cleanEEG.filename = strcat(extractBefore(filtEEG.filename,'.set'),'_CCA_EEGEMGcleanwNoise.set');
    %         cleanEEG.comments = pop_comments(EEG.comments,'','1Hz HPF; autoChRej (percentiles); CCA; EEG and EMG channels cleaned with noise channels, 2s window',1);
    %     elseif Y==3
    %         cleanEEG.filename = strcat(extractBefore(filtEEG.filename,'.set'),'_CCA_EEGEMGcleanwnoise_EEGcleanwEMG.set');
    %         cleanEEG.comments = pop_comments(EEG.comments,'','CCA; EEG and EMG channels cleaned noise channels, then EEG cleaned with EMG channels, 2s window, cleanYbool= true',1);
    %     end

    %% repeat for third pass
    %calculate spectra of cleaned EEG data
    %from 3rd pass (muscle artifacts removed)
    timeWindow = find(cleanEEG.times/1000 >= 0 & cleanEEG.times/1000 <= inf); %arbitrarily calculating PSD no smaller window so we don't have to wait forever
    [EEG_spectra_out,FREQ,~,~,~] = spectopo(cleanEEG.data(EEG_chans,timeWindow), length(timeWindow), cleanEEG.srate, 'freqfac',1, 'plot','off','percent',10);
    EEG_spectra_diff  = EEG_spectra_out-EEG_spectra_in;
    %cleanEEG.etc.specdata.spectra_of_tdiff = EEG_spectra_of_tdiff;
    specdata.spectra_clean = EEG_spectra_out;
    specdata.freq = FREQ;
    specdata.spectra_diff= EEG_spectra_diff;
    iCC_specdata.(testname).thirdpass.specdat = specdata;
    iCC_specdata.(testname).subject = filtEEG.subject;
    iCC_specdata.(testname).filename = filtEEG.filename;
    %% plot median PSD, third pass
    spec_plot= figure;%('WindowState','maximized'); % Plot spectra of mean channel input and output
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 12 6]);
    set(gcf,'Color','w');
    mytitle= title({'Median Spectral Power Density',...
        strcat(extractBefore(filtEEG.filename,'.set'),' ',testname,' ||  Noise and EMG Removed')},'Interpreter', 'none');
    hold on;
    plot(FREQ, median(EEG_spectra_in,1),'Color',cyan,'LineWidth', 2); plot(FREQ, prctile(EEG_spectra_in,25), '--','Color',cyan,'LineWidth', .5); plot(FREQ, prctile(EEG_spectra_in,75), '--','Color',cyan,'LineWidth', .5);
    set(gca,'FontSize',12);
    plot(FREQ, median(EEG_spectra_out,1),'Color',green,'LineWidth', 2); plot(FREQ, prctile(EEG_spectra_out,25), '--','Color',green,'LineWidth', .5); plot(FREQ, prctile(EEG_spectra_out,75), '--','Color',green,'LineWidth', .5);
    plot(FREQ, median(EEG_spectra_diff,1), 'k','LineWidth',1);
    legend('EEG in Q2', 'Q1','Q3', 'EEG out Q2','Q1','Q3','Difference between EEG IN vs OUT- Q2','Location', 'North', 'Orientation', 'horizontal');
    plot(FREQ, median(noise_spec_in,1) ,'Color',red,'LineWidth',1);
    plot(FREQ, median(EMG_spec_in,1) ,'Color','m','LineWidth',1);
    plot(FREQ, median(EMG_spec_out,1) ,'Color',purple,'LineWidth',1);
    legend('EEG IN, Q2', 'Q1','Q3', 'EEG OUT/CLEAN, Q2','Q1','Q3',...
        'Difference between EEG IN vs OUT- Q2','Raw Noise Spectra','EMG IN',...
        'EMG OUT/Clean','Location', 'northeast', 'Orientation', 'vertical');
    %             legend boxoff;
    xlabel('Frequency (Hz)');
    ylabel('Power 10*log_{10} (\muV^{2}/Hz)');
    hold off;
    xlim([.5 100]);  set(gca,'Box','off');
    %YLIM = get(gca,'YLim');
    ylim([-30 30])
    set(gca,'XTick', [0.5,4,8,13,30,60,100]);
    set(gca, 'XTickLabel', {'','4','8','13','30','60','100'});
    %% save figure
    cd(subFolder2)
    savethisfig(spec_plot,strcat(extractBefore(filtEEG.filename,'.set'),' ',testname,'-',mytitle.String{1}), subFolder2,'jpg');
    close;
    save([subFolder2,'\',strcat(extractBefore(EEG.filename,'_chanlocs'),'_iCC_specdata.mat')], 'iCC_specdata')
end

if removeChan_flag 
    cleanEEG = cleanEEG3;
%store  data from third pass
cleanEEG  = eeg_checkset( cleanEEG  );
cleanEEG.setname = strcat(cleanEEG.subject,' iCC');
if ~isempty(otherchanidx)  %keep external channels
    cleanEEG.data = [cleanEEG.data; rawEEG.data(otherchanidx,:)];
    cleanEEG.nbchan = size(cleanEEG.data,1); %update # of chans
    numNewSignals = length(otherchanidx);
    %Update channel labels and type - Forceplates,exo, leg EMG
    for newsig_i = 1:numNewSignals
        cleanEEG.chanlocs(end+1).labels = rawEEG.chanlocs(otherchanidx(newsig_i)).labels;
        cleanEEG.chanlocs(end).type = rawEEG.chanlocs(otherchanidx(newsig_i)).type;
    end
end

else
%replace good channels with clean data, but don't reject channels yet bc
%then you can't merge dataset
cleanEEG = rawEEG;
for chani =  1:cleanEEG.nbchan
    if any(strcmp(cleanEEG.chanlocs(chani).labels, {cleanEEG3.chanlocs.labels}))
        cleanEEG.data(chani,:) = cleanEEG3.data(strcmp(cleanEEG.chanlocs(chani).labels, {cleanEEG3.chanlocs.labels}),:) ;
    end
end
%store rejected channels
cleanEEG.reject.rejch = false([cleanEEG.nbchan,1]);
cleanEEG.reject.rejch([badEEGch,badEMGch],1) = 1;
end

cleanEEG.etc.iCanClean = output;
cleanEEG.iCanClean =[];

%remove noise channels, no longer needed
EEG = cleanEEG;
getchantypes;
fprintf('\nRemoving noise channels')
cleanEEG = pop_select(cleanEEG,'nochannel',sort( unique([Noise_chans]) ) ); %remove bad channels from data before rereference so only rereferencing once



fclose(notes);
%Measure elapsed time
tEnd = toc(tStart);
elapsed_min= tEnd/60; %minutes
elapsed_hr = 0;
if elapsed_min>120
    R = rem(elapsed_min,60);
    elapsed_hr = floor(elapsed_min/60); %hrs
    elapsed_min= round(R);
end
fprintf('iCC run time = %i hours and %.f minutes\n',elapsed_hr, elapsed_min)
diary OFF
end

% %% send email when finished
% msg = strcat('iCanClean Finished!! Total run time = ',num2str(elapsed_hr),...
%     ' hrs and ', num2str(elapsed_min), ' minutes. Files located at>> ',subFolder2);
% send_email('MISSION ACCOMPLISHED \(^o^)/', msg)
