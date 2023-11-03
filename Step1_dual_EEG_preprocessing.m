%% Exoskeleton Adapatation Pre-processing
%
% Author: Noelle Jacobsen, University of Florida
% Created 9/8/23, last updated
%
% Inputs:
%
%
% Ouputs:
%
%
% elecLocPath                   folder path containing all subject .txt
%                               files with digitized electrode locations
% codeRepo                      Code repository path

%% initialize paths
%set global variables
global mydir Exo_Gait_Events

mydir = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation';
eeglabpath = 'C:\Users\jacobsen.noelle\Desktop\eeglab2022.0';
codeRepo = 'R:\Ferris-Lab\jacobsen.noelle\Matlab_Scripts\ExoAdapt-EEG-Processing\';
addpath(genpath(codeRepo))
eleclocsFolder = [mydir,'\Data\raw_data\electrode_locations'];
EEGinputFolder= [mydir,'\Data\processed_data\BATCH-07-Mar-2023-Step1a'];
savePath = [mydir,'\Data\processed_data']; %subfolders will be created and added

%load Exo_Gait_Events table
% [ADD]



rawdatafolderList = dir([mydir,'\Data\raw_data']);
% outputFolder_trim= strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1d-Trim');
% if ~exist(outputFolder_trim, 'dir') %check to see if output folder exists, if not make new one
%     mkdir(outputFolder_trim);
% end
%
%
% cd(savePath)
% FILTERSPEC = '*.set';
% TITLE = 'Load EEG dataset';
% [ EEGfileList,EEGinputFolder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
% if iscell(EEGfileList) %array of chars
%     numfiles = length(EEGfileList);
% else
%     numfiles = size(EEGfileList,1);
% end

% %% Condition level processing
% clc;close all
% fin = []; count = 0;
% if ~exist('badsets','var')
%     badsets ={};
% end
% cd(savePath)
% diary([date,'-diary'])
% diary on
% fprintf('================================================================\n')
% fprintf('Step1_dual_EEG_preprocessing\n%s\n', date);
% fprintf('================================================================\n')
%
% for filei = 1:numfiles
%     mytic = tic;
%     cd(EEGinputFolder)
%     %file name depends on format of fileList
%     if iscell(EEGfileList) %array of chars
%         EEGfile = EEGfileList {filei};
%     elseif ischar(EEGfileList) %one char
%         EEGfile = EEGfileList;
%     elseif isstruct(EEGfileList)
%         EEGfile = EEGfileList(filei).name;
%     else
%         disp('Error in file type selected');
%     end
%     %  try
%     fprintf('\n\n%s',EEGfile)
%     fprintf('\n================================================================\n')
%     %% Load EEG dataset
%     if ~exist('ALLCOM')
%         addpath(eeglabpath)
%         [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %start EEGlab
%     end
%     EEG = pop_loadset('filename',EEGfile,'filepath',EEGinputFolder);
%     %add subject number
%     if isempty(EEG.subject)
%         subject= extractBefore(EEGfile,'_');
%         EEG.subject = subject;
%     end
%
%     %subject folder
%     subFolder = fullfile(rawdatafolderList(1).folder,rawdatafolderList(contains({rawdatafolderList.name},EEG.subject) & [rawdatafolderList.isdir]==1).name);
%     %% 1b) Add channel locations
%     if isempty(EEG.chanlocs(1).theta)
%         fprintf('File does not have channel locations, attempting to load locations...');
%         outputFolder= strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1b-Add_chanlocs');
%         if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
%             mkdir(outputFolder);
%         end
%         EEG = add_chanlocs(EEG,eleclocsFolder);
%         [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%         EEG = eeg_checkset( EEG );
%         EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_chanlocs.set'),'filepath',outputFolder);
%     end
%
%     %% 1c) Import gait events and other data channels
%      gaitdataFolder = [mydir,'\Data\processed_data\gait_data']; %folder containing gait data for each subject and condition
%     if isempty([EEG.event(strcmpi({EEG.event.type},'RHS'))])
%         outputFolder_gaitdat = strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1c-Import_gait_data');
%         if ~exist(outputFolder_gaitdat, 'dir') %check to see if output folder exists, if not make new one
%             mkdir(outputFolder_gaitdat);
%         end
%         syncB_Folder = [subFolder,'\Motive']; %folder containing system b .csv with sync signal, exported from Motive
%         newfilename = EEG.filename; %gait events .mat files were stored w/o '_chanlocs'
%         EEG.filename = [extractBefore(EEG.filename,'_chanlocs'),'.set'];
%         EEG = import_gaitdata(EEG,syncB_Folder, gaitdataFolder);
%         [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%         EEG.filename = newfilename;
%         EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_gaitdat.set'),'filepath',outputFolder_gaitdat);
%     end
%     %         %% 1d) Attenuate motion and muscle artifacts with iCanClean
%     %         opt.testname= 'EEGExo';
%     %         opt.iCCmethod =3;
%     %         opt.calc_psd= 1; %calculate PSD for input data
%     %         opt.visualizeFinalResults =1;
%     %          opt.removeChannels; %remove bad emg and eeg channels-- turn off if you haven't merged datasets yet.
%     %         cleanEEG = run_iCanClean(EEG,ALLEEG,CURRENTSET,outputFolder_iCC, opt);
%     %         [ALLEEG, EEG] = eeg_store(ALLEEG, cleanEEG, CURRENTSET);
%     %         EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_iCC.set'),'filepath',outputFolder_iCC);
%
    %% 1e) Crop data using gait events (gait data was already cropped to when
    % exo when recording, so gait events only occured then
    % add trigger as channel
%     EEG.data(end+1,:) = EEG.trigger;
%     EEG.chanlocs(end+1).labels = 'Sync';
%     EEG.chanlocs(end).type = 'Sync';
%     EEG.nbchan = size(EEG.data,1);
%
%     fprintf('\nCropping data to first and last RHS event')
%     RHS = [EEG.event(strcmpi({EEG.event.type},'RHS')).latency];
%     buffer = 5*EEG.srate; %5 second buffer before and after for time-freq analysis so no weird things happen at the edges
%     rejwin = [2 RHS(1)+buffer; RHS(end)+buffer EEG.pnts];
%     EEG = eeg_eegrej_2023( EEG, rejwin);
%     RHS = [EEG.event(strcmpi({EEG.event.type},'RHS')).latency];
%     plot_gaitevents(EEG,'time') %double check events because of event latency warning from eeg_eegrej
%     savethisfig(gcf,[EEG.subject,'_', EEG.condition,'_cropped'],[gaitdataFolder,'\gaitevents'],'fig');
%     if strcmp(EEG.condition,'pow_1')
%         xlim([EEG.times(1) EEG.times(10000)]/1000)
%     else
%         xlim([(RHS(end)/EEG.srate)-15 RHS(end)/EEG.srate])
%     end
%     savethisfig(gcf,[EEG.subject,'_', EEG.condition],[gaitdataFolder,'\gaitevents'],'jpg');
%     close;
%     [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%     EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_trim.set'),'filepath',outputFolder_trim);
%
%
%     %         %or crop manually
%     %         ogEEG =EEG;
%     %         EEG.data(end+1,:) = EEG.trigger;
%     %         EEG.chanlocs(end+1).labels = 'Sync';
%     %         EEG.nbchan = size(EEG.data,1);
%     %         pop_eegplot(EEG)
%
%     %         fprintf('\n\niCanClean Parameters:\n');
%     %         paramList = fields(cleanEEG.iCanClean.params);
%     %         for p = 1:length(paramList)
%     %             fprintf(notes,'\t-%s = %i\n',paramList{p},cleanEEG.iCanClean.params.(paramList{p}));
%     %         end
%
%     %
%     %     catch ME
%     %         fprintf('Skipping this subject due to issues: ')
%     %         fprintf('%s\n',ME.identifier)
%     %         badsets{end+1} = EEGfile;
%     %     end
%
%     %estimate time remaining
%     fprintf('\nFinished file %i/%i\n', filei,numfiles);
%     t_remaining(mytic,fin,count,numfiles)
%     clear cursor_info
%     close all;
% end
%
% completedfiles = dir(fullfile(outputFolder_trim,'*.set'));
% badseti= []; badsets ={};
% %other way to check for badsets
% for filei = 1:numfiles
%     filename = extractBefore(EEGfileList{filei},'.set');
%     if ~any(contains({completedfiles.name},filename))
%         badsets{end+1} = filename;
%         badseti = [badseti;filei];
%     end
% end
%
% if ~isempty(badsets)
%     fprintf('\nThese datasets were skipped:')
%     arrayfun(@(x) fprintf('%s\n',x{1,1}),badsets);
% end
% diary off

%% Subject level processing
inputFolder = outputFolder_trim;
outputFolder_merge = strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1e-Merged');
outputFolder_iCC = strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1f-iCanClean');
outputFolder_clean= strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1g-Fully_Preprocessed');
if ~exist(outputFolder_iCC, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder_iCC);
end
if ~exist(outputFolder_merge, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder_merge)
end
if ~exist(outputFolder_clean, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder_clean)
end

%     Get file related to each subject
subList = {};
badsubs ={};
EEGfileList = dir(fullfile(inputFolder,'*.set'));

for filei = 1:size(EEGfileList,1)
    subList{filei} = extractBefore(EEGfileList(filei).name,'_');
end
subList = unique(subList);

%%
cd(savePath)
diary 'preprocessing'
fprintf('================================================================\n')
fprintf('\t\tStep1_dual_EEG_preprocessing (steps e-f)\n\t\t%s\n', date);
fprintf('================================================================\n\n')
count = 0; fin = []; jobList ={};
for subi = 1:length(subList)
    try
        try
            fprintf('\n================================================================\n')
            fprintf('%s\n',subList{subi});
            fprintf('================================================================\n')
            subFiles = {EEGfileList(contains({EEGfileList.name},subList{subi})).name};

            %order conditions
            files_noexo = sort(subFiles(contains(subFiles,'noExo')));
            files_unpow = sort(subFiles(contains(subFiles,'unpow')));
            files_pow = sort(subFiles(contains(subFiles,'pow')& ~contains(subFiles,'un')));
            files_deadapt = sort(subFiles(contains(subFiles,'deadapt')));

            if isempty(files_noexo) || isempty(files_unpow) || isempty(files_pow) || isempty(files_deadapt)
                disp(newline)
                warning(['Skipping this subject (',subList{subi},') because they are missing conditions!!']); %warn user
                disp(newline)
                badsubs{end+1} = subList{subi};
                continue
            end

            subFiles =[files_noexo, files_unpow,files_pow,files_deadapt];
            numfiles = size(subFiles,2);

            eeglab;
            for filei = 1:numfiles
                %Load datasets
                fileName = subFiles{filei};
                EEG = pop_loadset('filename',fileName,'filepath',inputFolder);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); %append to new EEG to ALLEEG
                EEG.etc.iCanClean = []; %TMP
                EEG.urevent = [];
                [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = pop_saveset( EEG, 'savemode','resave','filename',EEG.filename,'filepath',EEG.filepath);
                %rejch(:,filei) = EEG.reject.rejch;%channel rejection mask, 0=keep, 1=reject
            end

            %% 1e) Merge files
            fprintf('Merging files in this order:')
            arrayfun(@(x) fprintf('%s\n',x{1,1}),subFiles);
            EEG = pop_mergeset( ALLEEG,1:numfiles, 0);
            EEG = eeg_checkset(EEG);
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = pop_saveset( EEG, 'filename',strcat(EEG.subject,'_allcond.set'),'filepath',outputFolder_merge);
            %     rejch = sum(rejch,2); %reject a channel if it was rejected in any condition file
            %     rejch = find(rejch>=1);
            %     if ~isempty(rejch)
            %     EEG = pop_select(EG,'nochannel',sort(unique(rejch))); %remove bad channels from data before rereference so only rereferencing once
            %     end

            %% 1f) Attenuate motion and muscle artifacts with iCanClean
            opt.testname= 'EEGExo';
            opt.iCCmethod =3; %3 cycles
            opt.calc_psd= 1; %calculate PSD for input data
            opt.visualizeFinalResults =1;
            opt.removeChannels = 1; %remove bad emg and eeg channels-- turn off if you haven't merged datasets yet.
            cleanEEG = run_iCanClean(EEG,ALLEEG,CURRENTSET,outputFolder_iCC, opt);
            cleanEEG = eeg_checkset(cleanEEG);
            [ALLEEG, EEG] = eeg_store(ALLEEG, cleanEEG, CURRENTSET);
            EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_iCC.set'),'filepath',outputFolder_iCC);

            %% 1g) Additional pre-processing
            %   clean line noise,rereference, channel and time window rejection
            cleanline_flag =1; %cleanline noise ON/OFF
            check_chan_std_flag =1; %check standard deviation of all channels and store figure
            useGPU = 0;
            EEG = cleanEEG_func(EEG,outputFolder_clean,cleanline_flag, check_chan_std_flag,useGPU,eleclocsFolder);
            EEG = eeg_checkset(EEG);
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'_iCC.set'),'_clean.set'),'filepath',outputFolder_clean);

            %% 1g) Prepare AMICA for Hipergator
            %parameters
            params = [];
            params.num_time = '00:45:00'; % Wall time hh:mm:ss
            params.pcaOverride =0 ; % override # of components to reduce using PC
            params.max_iter = 2000;
            params.channels = 'EEG';
            mainAMICAdir = '\\exasmb.rc.ufl.edu\blue\dferris\jacobsen.noelle';
            outputdir_unix =  '/blue/dferris/jacobsen.noelle/AMICA/';

            prepare_HPG_AMICA_func(EEG,extractBefore(EEG.filename,'.set'),...
                mainAMICAdir,outputdir_unix,...
                params,'jacobsen.noelle@ufl.edu')

            subDirNum = date;
            shortfileName = EEG.subject;
            AMICAdir_unix =[outputdir_unix,subDirNum,'/',shortfileName,'/'];
           jobList{end+1} = ['cd ' AMICAdir_unix];
            jobList{end+1} = 'sbatch runAMICA_hipergator.sh';

        catch ME
            fprintf('Skipping this subject due to issues: ')
            fprintf('%s\n',ME.identifier)
            badsubs{end+1} = subList{subi};
        end

        %estimate time remaining
        fprintf('\nFinished file %i/%i\n', subi,length(subList));
        t_remaining(mytic,fin,count,length(subList))
        close all;

    catch
        pause(60)
    end
end

if ~isempty(badsubs)
    fprintf('\nThese subjects were skipped:')
    arrayfun(@(x) fprintf('%s\n',x{1,1}),badsubs);
end
diary off


fprintf('\nJob List:\n')
arrayfun(@(x) fprintf('\t%s\n',x{1,1}),jobList);
disp('Done creating AMICA stuff');
%  setpref('Internet','E_mail','jacobnoe.reports@gmail.com');
%  sendmail('jacobsen.noelle@ufl.edu','Calculation complete.')



%=====================================================================
% for subjects that were halfway done...
%=====================================================================
%% just steps after merge
inputFolder = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation\Data\processed_data\2023-09-20_Step1f-iCanClean';
inputFolder = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation\Data\processed_data\2023-09-21_Step1g-Fully_Preprocessed';
outputFolder_iCC = strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1f-iCanClean');
outputFolder_clean= strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step1g-Fully_Preprocessed');
if ~exist(outputFolder_iCC, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder_iCC);
end
if ~exist(outputFolder_clean, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder_clean)
end

cd(inputFolder)
FILTERSPEC = '*.set';
TITLE = 'Load EEG dataset';
[ EEGfileList,EEGinputFolder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
if iscell(EEGfileList) %array of chars
    numfiles = length(EEGfileList);
else
    numfiles = size(EEGfileList,1);
end

%EEGfileList = dir(fullfile(inputFolder,'*.set'));
if ~exist('ALLCOM')
    addpath(eeglabpath)
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %start EEGlab
end
jobList = {};
%%
for filei =1:numfiles
    %file name depends on format of fileList
    if iscell(EEGfileList) %array of chars
        EEGfile = EEGfileList {filei};
    elseif ischar(EEGfileList) %one char
        EEGfile = EEGfileList;
    elseif isstruct(EEGfileList)
        EEGfile = EEGfileList(filei).name;
    else
        disp('Error in file type selected');
    end
    EEG = pop_loadset('filename',EEGfile,'filepath',EEGinputFolder);
    ALLEEG = EEG; CURRENTSET =1;
    %% 1f) Attenuate motion and muscle artifacts with iCanClean
    if ~contains(EEG.filename,'iCC') &&  ~contains(EEG.filename,'clean')
        opt.testname= 'EEGExo';
        opt.iCCmethod =3; %3 cycles
        opt.calc_psd= 1; %calculate PSD for input data
        opt.visualizeFinalResults =1;
        opt.removeChannels = 1; %remove bad emg and eeg channels-- turn off if you haven't merged datasets yet.
        cleanEEG = run_iCanClean(EEG,ALLEEG,CURRENTSET,outputFolder_iCC, opt);
        cleanEEG = eeg_checkset(cleanEEG);
        [ALLEEG, EEG] = eeg_store(ALLEEG, cleanEEG, CURRENTSET);
        EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_iCC.set'),'filepath',outputFolder_iCC);
    end
    %% 1g) Additional pre-processing
    %   clean line noise,rereference, channel and time window rejection
    if ~contains(EEG.filename,'clean')
        cleanline_flag =1; %cleanline noise ON/OFF
        check_chan_std_flag =1; %check standard deviation of all channels and store figure
        useGPU = 0;
        EEG = cleanEEG_func(EEG,outputFolder_clean,cleanline_flag, check_chan_std_flag,useGPU,eleclocsFolder);
        EEG = eeg_checkset(EEG);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'_iCC.set'),'_clean.set'),'filepath',outputFolder_clean);
    end
    %% 1g) Prepare AMICA for Hipergator
    %parameters
    params = [];
    params.num_time = '00:45:00'; % Wall time hh:mm:ss
    params.pcaOverride =0 ; %PCA value to use
    params.max_iter = 2000;
    params.channels = 'EEG';
    mainAMICAdir = '\\exasmb.rc.ufl.edu\blue\dferris\jacobsen.noelle';
    outputdir_unix =  '/blue/dferris/jacobsen.noelle/AMICA/';

    prepare_HPG_AMICA_func(EEG,extractBefore(EEG.filename,'.set'),...
        mainAMICAdir,outputdir_unix,...
        params,'jacobsen.noelle@ufl.edu')

    subDirNum = date;
    shortfileName = EEG.subject;
    AMICAdir_unix =[outputdir_unix,subDirNum,'/',shortfileName,'/'];
    jobList{end+1} = ['cd ' AMICAdir_unix];
    jobList{end+1} = 'sbatch runAMICA_hipergator.sh';
end

fprintf('\nJob List:\n')
arrayfun(@(x) fprintf('\t%s\n',x{1,1}),jobList);
disp('Done creating AMICA stuff');
