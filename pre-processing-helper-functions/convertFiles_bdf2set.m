%% Converting of .bdf files into .set and .fdt
% Step 1: Import Data
    %   Starting EEGlab dataset from scratch (.bdf files)
    %   Noelle Jacobsen
    %   3/18/21
    %  
    %   This function loads BDF files and exports datasets to a new local folder

%clc; clear all; close all;
%% Enter Paramters
eeglabpath = 'C:\Users\jacobsen.noelle\Desktop\eeglab2022.0';
codepath = 'R:\Ferris-Lab\jacobsen.noelle\Matlab_Scripts\ExoAdapt-EEG-Processing';
addpath(codepath)
inputFolder = uigetdir('')
outputFolder= 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation\Batch Processing\';
bad_sets = {}; %keep track of datasets with errors
tStart = tic;
subjects = input('Enter Subject number: '); %Manually input subject number (with no '0')
%subjects= []; %numeric range or leave empty to proceess all subjects
folders = dir(inputFolder);
folderi = find(contains({folders.name}, 'PHD'));
count = 1;

%process all subject folders that start with PhDXX_SXX
if isempty(subjects)
    fprintf('Folder is empty ');
end
	
%% Load bdf files based on subject number
for subject_num = subjects
    tic
    clear ALLCOM EEG ALLEEG LASTCOM PLUGINLIST STUDY CURRENTSET CURRENTSTUDY
    phd_num = subject_num+14;
    subject_num= num2str(subject_num);
    phd_num = num2str(phd_num);
 
if (size(subject_num,2)==1)
    subject_num= strcat('PHD',phd_num,'_S0',subject_num);
elseif (size(subject_num,2)==2)
    subject_num= strcat('PHD',phd_num,'_S',subject_num);
end
fprintf('Starting subject %s\n', subject_num);
%inputFolder = strcat('R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation\Data\',subject_num,'\EEG');
cd(inputFolder)
fileList = dir('*.bdf'); %loads all .bdf files in folder
%OR Load files manually
FILTERSPEC = '*.bdf';
TITLE = 'Choose a dataset to load';
[fileList inputFolder] = uigetfile(FILTERSPEC, TITLE,'MultiSelect', 'on');
cd(inputFolder)
%% Check # of files selected
if iscell(fileList) == 1
    n = length(fileList);
elseif isstruct(fileList)
    n = size(fileList,1);
else 
    n=1;
end
%% Start up EEGLAB if not running already
if ~exist('ALLCOM')
    addpath(eeglabpath)
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %start EEGlab
end

%% Create datasets from BDF files
for i=1:n
    %file name depends on format of fileList
    if iscell(fileList) %array of chars
        file = fileList{i};
    elseif ischar(fileList) %one char
        file = fileList;
    elseif isstruct(fileList)
        file = fileList(i).name;
    else
        disp('Error in file type selected');
    end
fprintf('Reading file %i/%i: ',i,n);
disp(file); fprintf('  ...');
cd(inputFolder);
[EEG, command, H] = pop_biosig(strcat(inputFolder,'\',file)); %must modify pop_biosig so that you save variable "H"
EEG.H = H;
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
EEG = eeg_checkset( EEG );
%eeglab redraw;
%% Fix event labels 
%Plug-ins can mislabel events when importing BDFs. This fixes event labels
%by reading the 24bit data from Biosemi USB trigger interface
EEG.event=[];
EEG.urevent=[];
try
[EEG] = fix_popsig_eventlabels(EEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG,0);
catch
    bad_sets = {bad_sets; file};
end
%figure; plot (EEG.trigger); title(EEG.filename);
% %% Remove excess channels (only keep 128 scalp, 8 EMG, and 128 noise)
% EEG = eeg_checkset( EEG );
% EEG = pop_select( EEG,'channel',{'1-A1' '1-A2' '1-A3' '1-A4' '1-A5' '1-A6' '1-A7' '1-A8' '1-A9' '1-A10' '1-A11' '1-A12' '1-A13' '1-A14' '1-A15' '1-A16' '1-A17' '1-A18' '1-A19' '1-A20' '1-A21' '1-A22' '1-A23' '1-A24' '1-A25' '1-A26' '1-A27' '1-A28' '1-A29' '1-A30' '1-A31' '1-A32' '1-B1' '1-B2' '1-B3' '1-B4' '1-B5' '1-B6' '1-B7' '1-B8' '1-B9' '1-B10' '1-B11' '1-B12' '1-B13' '1-B14' '1-B15' '1-B16' '1-B17' '1-B18' '1-B19' '1-B20' '1-B21' '1-B22' '1-B23' '1-B24' '1-B25' '1-B26' '1-B27' '1-B28' '1-B29' '1-B30' '1-B31' '1-B32' '1-C1' '1-C2' '1-C3' '1-C4' '1-C5' '1-C6' '1-C7' '1-C8' '1-C9' '1-C10' '1-C11' '1-C12' '1-C13' '1-C14' '1-C15' '1-C16' '1-C17' '1-C18' '1-C19' '1-C20' '1-C21' '1-C22' '1-C23' '1-C24' '1-C25' '1-C26' '1-C27' '1-C28' '1-C29' '1-C30' '1-C31' '1-C32' '1-D1' '1-D2' '1-D3' '1-D4' '1-D5' '1-D6' '1-D7' '1-D8' '1-D9' '1-D10' '1-D11' '1-D12' '1-D13' '1-D14' '1-D15' '1-D16' '1-D17' '1-D18' '1-D19' '1-D20' '1-D21' '1-D22' '1-D23' '1-D24' '1-D25' '1-D26' '1-D27' '1-D28' '1-D29' '1-D30' '1-D31' '1-D32' '1-EXG1' '1-EXG2' '1-EXG3' '1-EXG4' '1-EXG5' '1-EXG6' '1-EXG7' '1-EXG8' '2-A1' '2-A2' '2-A3' '2-A4' '2-A5' '2-A6' '2-A7' '2-A8' '2-A9' '2-A10' '2-A11' '2-A12' '2-A13' '2-A14' '2-A15' '2-A16' '2-A17' '2-A18' '2-A19' '2-A20' '2-A21' '2-A22' '2-A23' '2-A24' '2-A25' '2-A26' '2-A27' '2-A28' '2-A29' '2-A30' '2-A31' '2-A32' '2-B1' '2-B2' '2-B3' '2-B4' '2-B5' '2-B6' '2-B7' '2-B8' '2-B9' '2-B10' '2-B11' '2-B12' '2-B13' '2-B14' '2-B15' '2-B16' '2-B17' '2-B18' '2-B19' '2-B20' '2-B21' '2-B22' '2-B23' '2-B24' '2-B25' '2-B26' '2-B27' '2-B28' '2-B29' '2-B30' '2-B31' '2-B32' '2-C1' '2-C2' '2-C3' '2-C4' '2-C5' '2-C6' '2-C7' '2-C8' '2-C9' '2-C10' '2-C11' '2-C12' '2-C13' '2-C14' '2-C15' '2-C16' '2-C17' '2-C18' '2-C19' '2-C20' '2-C21' '2-C22' '2-C23' '2-C24' '2-C25' '2-C26' '2-C27' '2-C28' '2-C29' '2-C30' '2-C31' '2-C32' '2-D1' '2-D2' '2-D3' '2-D4' '2-D5' '2-D6' '2-D7' '2-D8' '2-D9' '2-D10' '2-D11' '2-D12' '2-D13' '2-D14' '2-D15' '2-D16' '2-D17' '2-D18' '2-D19' '2-D20' '2-D21' '2-D22' '2-D23' '2-D24' '2-D25' '2-D26' '2-D27' '2-D28' '2-D29' '2-D30' '2-D31' '2-D32'});
% %label channel types
%  for v= 1:length(EEG.chanlocs)
%             if contains({EEG.chanlocs(v).labels}, 'EMG')
%                 EEG.chanlocs(v).type = 'EMG';
%             elseif contains({EEG.chanlocs(v).labels}, 'EXG')
%                 EEG.chanlocs(v).type = 'EMG';
%             elseif contains({ EEG.chanlocs(v).labels},'1-') && ~contains({EEG.chanlocs(v).labels}, 'EXG')
%                 EEG.chanlocs(v).type = 'EEG';
%             elseif contains({ EEG.chanlocs(v).labels},'2-')
%                  EEG.chanlocs(v).type = 'noise';
%             end
%  end
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET);
% EEG = eeg_checkset( EEG );
%% Save dataset
   if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
       mkdir(outputFolder);
   end
   cd(outputFolder);
   %% extract filename -- exclude phdXXX extractAfter()
if file(1) == 'P'
   file = extractAfter(file, 6);
end
file = extractBefore(file,'.bdf');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',file,'gui','off');
%eeglab redraw;
fprintf('New dataset saved\n')
end
disp('All datasets for this subject saved')
%estimate time remaining
fin = toc;
est = fin*(length(subjects)-count)/60;
count = count+1; %keep track of how many subjects we've processed
fprintf('Finished subject %s\n', subject_num);
fprintf('Estimated time remaining: %i minutes\n', round(est));

end
 
if ~isempty(bad_sets)
fprintf('WARNING: Need to fix %i datasets', length(bad_sets)-1);
disp(bad_sets)
else
    fprintf('DONE\n')
end

%Measure elapsed time
tEnd = toc(tStart);
elapsed_min= tEnd/60; %minutes
elapsed_hr = 0;
if elapsed_min>120
    R = rem(elapsed_min,60);
    elapsed_hr = floor(elapsed_min/60); %hrs
    elapsed_min= round(R);
end
fprintf('Total run time = %i hours and %.f minutes',elapsed_hr, elapsed_min); 
