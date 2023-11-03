%% SAVE COMP DATA AS MAT FILE FOR CUSTOM HEAD MODEL FITTING
% Add AMICA results from Hipergator and save to new dataset in local folder
% Inputs
% mydir              (string) main directory where output folder will be stored
% subjects           (array) numeric array of subject numbers   
% downsample flag    (int) [0|1] option to downsample data
% newFs              (num) new sampling rate. Only works if flag is set ON

%set parameters
mydir = 'R:\Ferris-Lab\jacobsen.noelle\Split Belt Pilot Study\';
downsample_flag = 1; 
newFs = 256;
AMICA15_path = 'Z:\jacobsen.noelle\AMICA\AMICA_15';
subjects = [34];
outputFolder = strcat(mydir,'BATCH-',datestr(now, 'yyyy-mm-dd'),'-Step4-AMICA');


folderList ={}; 
TITLE = 'Load EEG folders from HPG';
[EEGfolder] = uigetdir('Z:\jacobsen.noelle\AMICA', TITLE);
x=0;
for subject_num = [subjects]  
  subject_num = num2str(subject_num);
    x= x+1;
    folderList{x} = strcat(EEGfolder,'\S',subject_num);
end

if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder);
end
addpath(AMICA15_path) 

% folderList ={};
% i=1;
% for subject_num = [42]
%     folderList{i} = strcat('Z:\jacobsen.noelle\AMICA\13-May-2021\S',num2str(subject_num),'_splitbelt\');
%   i=i+1;
% end
folderList = cellstr(folderList);
%% loop    
for X = 1:length(folderList)
%load EEG file
  if iscell(folderList) %array of chars
        folder = folderList {X};
    elseif ischar(folderList) %one char
        folder = folderList(X,:);
    elseif isstruct(folderList)
        folder = folderList(X).name;
    else
        disp('Error in file type selected');
  end
amicaFolder = convertCharsToStrings(folder);
%amicaFolder = strcat('M:\dferris\jacobsen.noelle\AMICA\',subject_num,'.',num);
%amicaFolder = strcat('Z:\jacobsen.noelle\AMICA\51221\',folder);
temp = dir(fullfile(amicaFolder,'*.set'));

if length(temp)==1
    fileName = temp.name;
    disp(fileName);
else
    disp('figure this out later'); %ask person to click on file if file doesn't exist or more than one exists
end

    %%  Open EEGlab (if not already open) and load EEG files
    if ~exist('ALLCOM')
        addpath('C:\Users\jacobsen.noelle\Desktop\eeglab2021.0');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    end
EEG = pop_loadset('filename',fileName,'filepath',temp.folder); 

%check downsample param
if downsample_flag
    if newFs> EEG.srate
        error('new sampling rate is larger than current sampling rate. Cannot downsample')
    end
end

%get rid of old ICA results from infomax
EEG = rmfield(EEG,'icaweights');                                                
EEG.icaweights = [];
EEG = rmfield(EEG,'icawinv'); 
EEG.icawinv = [];
EEG = rmfield(EEG,'icaact');
EEG.icaact = [];
EEG = rmfield(EEG,'icasphere');
EEG.icasphere = [];

%load up AMICA results
EEG.etc.amicaResultStructure = loadmodout15(temp.folder);


%Insert AMICA results into EEGLAB
EEG.icaweights = EEG.etc.amicaResultStructure.W;
EEG.icawinv = EEG.etc.amicaResultStructure.A;

num_pcs = EEG.etc.amicaResultStructure.num_pcs; %this is to account for doing PCA reduction 
EEG.icasphere  = EEG.etc.amicaResultStructure.S(1:num_pcs,:); %if you do PCA reduction and send the full S matrix,  you will get an error. Someone online said you should take the first X rows where X = your PCA reduction
num_components =  num2str(size(EEG.icaweights,1));

%make sure everything is ok and calculate the activation matrix
pop_editoptions('option_computeica',1); %randomly needed for some people to run even if you manually selected this checkbox with the GUI
EEG = eeg_checkset(EEG, 'ica'); %Note: I have no idea why it says "LLt not set" or what LLt is

%store updated EEG structure
AMICA = size(EEG.icaact,1);
AMICA = num2str(AMICA);
EEG.filename = strcat(extractBefore(EEG.filename,'.set'),'_AMICA',AMICA);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
eeglab redraw

%downsample
if downsample_flag
    if EEG.srate > newFs
    EEG = pop_resample(EEG,newFs);
    EEG.comments = pop_comments(EEG.comments,'',strcat('Downsampled to ',newFs,' Hz'),1);
    EEG.setname = strcat([EEG.subject,'','AMICA- downsamp']);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_DS','.set'),'filepath',outputFolderDS,'version','7.3');
    eeglab redraw;
    fprintf('Downsampled to %i',newFs);
    end

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );
eeglab redraw

%save comp data in fieldtrip format for use in custom head model later
comp = eeglab2fieldtrip(EEG, 'componentanalysis');
cd(amicaFolder)
save('comp.mat','comp'); clear comp

%save dataset
EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',outputFolder, 'version', '7.3');

%update GUI
%eeglab redraw

end
   
   
   