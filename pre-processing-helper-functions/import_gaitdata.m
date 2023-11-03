% import_gaitdata()
% 




function EEG = import_gaitdata(EEG,syncB_Folder,gaitdataFolder)
%Parameters
syncB_columnHeader = "Dev1_ai0"; % column header name for sync signal in .csv from motive
autoload_syncB = 1; % [1|0] automatically load csv with sync signal from motive. 0 = manaully select file
sync_freq = 0.25; %sync pulse frequency (Hz)
sync_duty = 50; % sync duty cycle (%)
sync_amp = 3.3; %sync pulse amplitude
GRF_thresh =10;
%% Load Gait Events (DAQB) (pre-labeled from GRF data by Rachel, stored in table)
%find subject and condition index in data structure
condition = extractAfter(EEG.filename,[EEG.subject,'_']);
condition = extractBefore(condition,'.set');
if contains(condition,'pow') || contains(condition,'noExo')
  EEG.condition = condition(1:5);
else
  EEG.condition = condition;
end



if strcmp(condition,'noExo')%adjust condition names to match what's in Exo_Gait_Data
    gaitdata_cond = 'no_exo';
elseif contains(condition,'pow') && ~contains(condition,'unpow')
    gaitdata_cond = ['exo_',condition];
else
    gaitdata_cond  = condition;
end

gaitdata_filename = [EEG.subject,'_',gaitdata_cond,'.mat'];
gaitdata_fileList = dir(gaitdataFolder);
filei = find(strcmp({gaitdata_fileList.name},gaitdata_filename));
if isempty(filei)
    error('Can''t find gait events file for this subject: %s',gaitdata_filename)
else
    gaitdata = load([gaitdataFolder,'\',gaitdata_fileList(filei).name]);
    gaitdata = gaitdata.dat;
end
%% Get EEG-Motive Sync Events
if isfield(gaitdata,'tSyncEvents_EEG') && isfield(gaitdata,'tSyncEvents_Motive') %check if time mapping was already done and saved
    tSyncEvents_EEG = gaitdata.tSyncEvents_EEG;
    tSyncEvents_Motive = gaitdata.tSyncEvents_Motive;

elseif isfield(gaitdata,'tsyncA_EEG') && isfield(gaitdata,'tsyncB_Motive') %check if time mapping was already done and saved
    tSyncEvents_EEG = gaitdata.tsyncA_EEG; %rename fields
    tSyncEvents_Motive = gaitdata.tsyncB_Motive;
    tsyncA_EEG =[];
    tsyncB_Motive =[];
else
    % Find EEG Sync Events (DAQA)
    %find trigger/sync events that were labeled in fix_popbiosig_eventlabels.m
    %just use findpeaks so that EEG events are the same as motive

    % EEG_sync_lat=[EEG.event(find(strcmp({EEG.event.type},'TRIG') | strcmp({EEG.event.type},'SYNC'))).latency]; % finds the latency of the trigger indicating the start of the trial in EEG samples (probably 512 Hz
    [~,EEG_sync_lat]=findpeaks(abs(EEG.trigger), 'MinPeakHeight', .2, 'MinPeakDistance', 150 ); %Matches Rachel's way of finding events in motive sync
    EEG_sync_events_t = EEG.times(EEG_sync_lat)/1000; %convert from ms to sec

    % quick quality check
    [EEG_sync_events_t] = syncQC(EEG_sync_events_t,sync_freq);

    % Find Motive Sync Events (DAQB)
    gaitevents = table2struct(gaitdata.GaitEvents);
    syncB_events_t = [gaitevents(strcmp({gaitevents.EventType},'sync')).Time];
    fprintf('size of syncBevents from gait event struct: %i\n',size(syncB_events_t,2))

    %load sync signal from system B just to be sure
    if ~autoload_syncB
        %Upload motive trigger file, make sure you have the current folder that contains the excel sheets open in Matlab
        FILTERSPEC = '*.csv';
        TITLE = 'Choose excel file with trigger data (ends in 29751840)';
        [syncB_filename, syncB_Folder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
    else
        fileList = dir(syncB_Folder);
        try %deal with file naming convention issues
            syncB_filename = fileList(strcmpi({fileList.name},strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_', extractBefore(EEG.filename, '.set'),'__29751840.csv') ) & [fileList.isdir]==0).name;
        catch
            try
                condnum = extractAfter(condition, 'pow_');
                cond_prefix = extractBefore(condition,'_');
                tmpname = strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_',EEG.subject,'_',cond_prefix,condnum,'__29751840.csv');
                syncB_filename = fileList(strcmpi({fileList.name},tmpname) & [fileList.isdir]==0).name;

            catch
                fprintf('Can''t find motive file: %s\n', strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_', extractBefore(EEG.filename, '.set'),'__29751840.csv') )
                fprintf('Please manually select motive trigger file for %s', EEG.filename)
                FILTERSPEC = '*.csv';
                TITLE = 'Choose excel file with trigger data (ends in 29751840)';
                cd(syncB_Folder)
                [syncB_filename, syncB_Folder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
            end
        end
    end
    cd(syncB_Folder)
    syncB_table = readtable(syncB_filename);
    deviceFrame = syncB_table.DeviceFrame;
    MocapTime = syncB_table.MocapTime;
    DAQB_Fs = deviceFrame(MocapTime==1)- deviceFrame(MocapTime==0); %calculate sampling rate
    DAQB_Fs = DAQB_Fs(end);
    fprintf('\nMotive_Fs = %i Hz\n',DAQB_Fs)
    syncB_t = syncB_table.DeviceFrame/DAQB_Fs;
    try
        syncB = syncB_table.(syncB_columnHeader);
    catch
        try
            syncB = syncB_table.Dev1_ai1;
        catch
            syncB = syncB_table.Sync;
        end
    end

    if syncB_t(end)<300 %check trial length
        warning('Motive trial seems too short')
    end
    clear syncB_table MocapTime deviceFrame

    % quick quality check
    if ~isempty(syncB_events_t)
        syncB_events_t = syncQC( syncB_events_t ,sync_freq);
    end

    % check sync event sizes from system A&B
    if size( syncB_events_t,2) == size(EEG_sync_events_t,2) && size(EEG_sync_events_t,2)>=100 %need at least two events, 100 just to be sure

        disp('Great! Sync events already match')
        tSyncEvents_Motive = syncB_events_t ;
        tSyncEvents_EEG = EEG_sync_events_t;
        tshift = EEG_sync_events_t(1)-syncB_events_t(1); %shift DAQB time for plotting purposes

        figure; tiledlayout(2,1)
        ax1=nexttile;
        plot(EEG.times/1000, EEG.trigger);
        hold on; stem(tSyncEvents_EEG, 1.2*ones(size(tSyncEvents_EEG)), 'k');
        ylim([-0.5 1.5]);
        xlabel('Time (s)'); ylabel('Digital Signal'); title('EEG Sync')
        legend({'signal','sync events'},'Location','southeast')

        ax2 = nexttile; plot(syncB_t+tshift, syncB); hold on;
        stem(tSyncEvents_Motive+tshift, 4*ones(size(tSyncEvents_Motive)), 'k'); hold off;
        linkaxes([ax1 ax2],'x')
        ylim([-1 5]);
        xlabel('Time (s)'); ylabel('Voltage'); title('Motive Sync')
        legend({'sync signal','sync events'},'Location','southeast')

    else    %deal with wonky motive signal
        fprintf('Number of sync events automatically found don''t match. All hands on deck! (manual inspection)\n')
        [tSyncEvents_EEG, tSyncEvents_Motive] = findSyncBuddies(EEG, syncB, syncB_t, sync_freq, sync_amp, sync_duty);
    end

    if size(tSyncEvents_Motive,1) ~=1; tSyncEvents_Motive = tSyncEvents_Motive';end
    if size(tSyncEvents_EEG,1) ~=1; tSyncEvents_EEG = tSyncEvents_EEG';end

    % save
    savethisfig(gcf,[EEG.subject,'_',condition,'_syncbuddies'],[gaitdataFolder,'\sync_qc\syncbuddies\fig'],'fig')
    savethisfig(gcf,[EEG.subject,'_',condition,'_syncbuddies'],[gaitdataFolder,'\sync_qc\syncbuddies\jpg'],'jpg')
    close;
    gaitdata.tSyncEvents_EEG =  tSyncEvents_EEG;
    gaitdata.EEG_Fs = EEG.srate;
    gaitdata.tSyncEvents_Motive = tSyncEvents_Motive;
    gaitdata.Motive_Fs = DAQB_Fs;
    dat = gaitdata;
    save([gaitdataFolder,'\',gaitdata_fileList(filei).name],"dat",'-append');
end

% %% Load COP daata
%     %load sync signal from system B just to be sure
%     if ~autoload_COP
%         %Upload motive trigger file, make sure you have the current folder that contains the excel sheets open in Matlab
%         FILTERSPEC = '*.csv';
%         TITLE = 'Choose excel files wih forceplate data ';
%         [syncB_filename, syncB_Folder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
%     else
%         fileList = dir(syncB_Folder);
%         try %deal with file naming convention issues
%             syncB_filename = fileList(strcmpi({fileList.name},strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_', extractBefore(EEG.filename, '.set'),'__29751840.csv') ) & [fileList.isdir]==0).name;
%         catch
%             try
%                 condnum = extractAfter(condition, 'pow_');
%                 cond_prefix = extractBefore(condition,'_');
%                 tmpname = strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_',EEG.subject,'_',cond_prefix,condnum,'__29751840.csv');
%                 syncB_filename = fileList(strcmpi({fileList.name},tmpname) & [fileList.isdir]==0).name;
% 
%             catch
%                 fprintf('Can''t find motive file: %s\n', strcat('PHD',string(str2num(extractAfter(EEG.subject,'S'))+14),'_', extractBefore(EEG.filename, '.set'),'__29751840.csv') )
%                 fprintf('Please manually select motive trigger file for %s', EEG.filename)
%                 FILTERSPEC = '*.csv';
%                 TITLE = 'Choose excel file with trigger data (ends in 29751840)';
%                 cd(syncB_Folder)
%                 [syncB_filename, syncB_Folder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
%             end
%         end
%     end
%     cd(syncB_Folder)
% 
% 
%    fp1 = readtable(COP_filename);
%    fp2 = readtable('PHD15_S01_deadapt_forceplate_2');
%    cop1 = sqrt(fp1.Cx.^2 + fp1.Cy.^2);
%    cop2 = sqrt(fp2.Cx.^2 + fp2.Cy.^2);
%    x  = 1:length(cop1);
%    x = x';
%    x1 = 60000;
%    cop1_dydx = diff(cop1)./diff(x);
%    cop2_dydx = diff(cop2)./diff(x);
%    x = size(cop1_dydx);
%    cop1_dy2dx = diff(cop1_dydx)./diff(x);
%    cop2_dy2dx = diff(cop2_dydx)./diff(x);
% 
%    GRF_L = fp1.Fz;
%    GRF_R = fp2.Fz;
% 
% 
%     figure; tiledlayout(2,1);
%     ax1 = nexttile;
%     plot(GRF_L(1:x1)); hold on;
%     plot(cop1(1:x1))
%     plot(cop1_dydx(1:x1))
%     plot(cop1_dy2dx(1:x1))
%     title('Left')
%     
%     ax2 = nexttile;
%      plot(GRF_R(1:x1)); hold on;
%     plot(cop2(1:x1))
%     plot(cop2_dydx(1:x1))
%     plot(cop2_dy2dx(1:x1))
%      linkaxes([ax1, ax2],'xy')
%     title('Right')
%     legend({'Fz','cop','cop_dydx','cop_dy2dx'}, 'interpreter','none')

% =======================================================================
% Add Gait Data to EEG
% =======================================================================
%
%% Add addition data channels to EEG (forceplates, exo)
if size(tSyncEvents_Motive,1) ~=1
    tSyncEvents_Motive =tSyncEvents_Motive';
end
%Check sync events
tSyncEvents_EEG = sort(tSyncEvents_EEG); %time of EEG sync events in seconds
tSyncEvents_Motive = sort(tSyncEvents_Motive);%time of Motive sync events in seconds
if size( tSyncEvents_EEG) ~= size(tSyncEvents_Motive)
    error('Number of sync events for EEG and Motive do not match')
end
% storesync info in EEG
EEG.etc.TimeMapping(1).Condition = EEG.condition;
EEG.etc.TimeMapping(1).tSyncEvents_EEG = tSyncEvents_EEG;
EEG.etc.TimeMapping(1).tSyncEvents_Motive =tSyncEvents_Motive;

%load  data
DAQB_data = table2struct(gaitdata.Biomech_Data_Tseries);
DAQB_t = [DAQB_data.Time];
fields = fieldnames(DAQB_data);

datain =[];
for fieldi = 1:length(fields)
    if ~strcmp(fields{fieldi},'Time')
       % datain = [datain; DAQB_data.(fields{fieldi})];
        EEG.chanlocs(end+1).labels = fields{fieldi};
        if  contains(fields{fieldi},'Filt') ||  contains(fields{fieldi},'Raw')
            EEG.chanlocs(end).type = 'Leg-EMG';
        elseif contains(fields{fieldi},'GRF')
            EEG.chanlocs(end).type = 'GRF';
        else
            EEG.chanlocs(end).type = 'Exo';
        end

    end
end
dataout = AlignDAQ2DAQ(EEG.srate, EEG.times/1000, DAQB_t, tSyncEvents_EEG, tSyncEvents_Motive, datain); %both time vectors must be the same units
EEG.data = [EEG.data; dataout];
EEG.nbchan = size(EEG.data,1);
EEG = eeg_checkset(EEG);

%check that GRF data was zeroed
GRF_L = EEG.data(strcmp({EEG.chanlocs.labels},'Left_GRF'),:);
GRF_R = EEG.data(strcmp({EEG.chanlocs.labels},'Right_GRF'),:);
if min(GRF_L)<-50 %rezero if minimum is less than -50N
    %GRF_L = GRF_L-min(GRF_L);
    GRF_L = GRF_L-prctile(GRF_L,5);
    fprintf('\nRezeroing left GRF')
    EEG.data(strcmp({EEG.chanlocs.labels},'Left_GRF'),:) = GRF_L;
end
if min(GRF_R)<-50 %rezero if minimum is less than -50N
    %GRF_R = GRF_R-min(GRF_R);
    GRF_R = GRF_R-prctile(GRF_R,5);
    fprintf('\nRezeroing right GRF')
    EEG.data(strcmp({EEG.chanlocs.labels},'Right_GRF'),:) = GRF_R;
end

%% Load gait gait events
% get gait events stored in external structure
%gaitevents(strcmpi([gaitevents.EventType],'sync')) = []; %remove, not using these events
% t_DAQb_in_DAQa_base = interp1(tSyncEvents_Motive,tSyncEvents_EEG, [gaitevents.Time],'linear','extrap');
% %verify mapping is correct
% figure;
% plot(tSyncEvents_Motive,tSyncEvents_EEG,'o',[gaitevents.Time],t_DAQb_in_DAQa_base,':.');
% title('Mapping between EEG time and Motive time');
% xlabel('Motive local time (s)'); ylabel('EEG local time (s)'); legend('Events','Interpolated');
% savethisfig(gcf,[EEG.subject,'_',condition],[gaitdataFolder,'\sync_qc\AlignDAQ2DAQ'],'jpg')
% close;
%% Add events to EEG.events
% for eventi = 1:size(gaitevents,1)
%     if ~isempty(gaitevents(eventi).EventType)
%     gaitevents(eventi).EEG_Time = t_DAQb_in_DAQa_base(eventi);
%     if gaitevents(eventi).EEG_Time >=0 && gaitevents(eventi).EEG_Time <=EEG.times(end)/1000
%     EEG.event(end+1).latency = t_DAQb_in_DAQa_base(eventi)*EEG.srate; %units should be in samples/frames
%     EEG.event(end).type = gaitevents(eventi).EventType{1,1}; %note that we use end here instead of end+1 since the length of EEG.event has been updated
%     EEG.event(end).cond = EEG.condition;
%     end
%     end
% end

%% get gait events from raw data
cd(gaitdataFolder)
gaitevents = getgaitevents(EEG,GRF_thresh, gaitdataFolder);

%% add events to EEG.event
for eventi = 1:length(gaitevents)
    if ~isempty(gaitevents(eventi).type)
    EEG.event(end+1).latency = round(gaitevents(eventi).latency); %units should be in samples/frames
    EEG.event(end).type = gaitevents(eventi).type; %note that we use end here instead of end+1 since the length of EEG.event has been updated
    EEG.event(end).cond = EEG.condition;
    end
end

EEG.event(strcmpi({EEG.event.type},'sync')) = []; %remove, not using these events
for eventi = 1:size(tSyncEvents_EEG,2)
    EEG.event(end+1).latency = tSyncEvents_EEG(eventi)*EEG.srate; %units should be in samples/frames
    EEG.event(end).type = 'SYNC'; %note that we use end here instead of end+1 since the length of EEG.event has been updated
    EEG.event(end).cond = EEG.condition;
end
EEG = eeg_checkset(EEG,'eventconsistency');

%% quality check
plot_gaitevents(EEG,'time');
savethisfig(gcf,[EEG.subject,'_',condition],[gaitdataFolder,'\gaitevents\fig'],'fig')
close;

EEG.comments = pop_comments(EEG.comments,'',strcat('-extracted gait events with','', num2str(GRF_thresh),' (N) threshold'),1);
EEG.setname = [EEG.subject,'-','gait events'];
%[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
end


%% sync event quality check
function [goodsyncEvents] = syncQC(syncEvents,sync_freq)
%sync quality check-- remove sync events that happened too close together
goodsyncEvents = syncEvents(1);
for eventi = 2:length(syncEvents) %assume first event is accurate
    timediff = syncEvents(eventi)-goodsyncEvents(end);
    if timediff > (1/sync_freq - (1/sync_freq)*0.05) && timediff < (1/sync_freq + (1/sync_freq)*0.05) %check that event interval is reasonable given frequency (+ 5%)
        goodsyncEvents = [goodsyncEvents, syncEvents(eventi)];
    end
end
end