%% Start EEGLAB
eeglab;


%% Load EEG file
amicaFolder = 'M:\dferris\jacobsen.noelle\AMICA\9.12\';

temp = dir(fullfile(amicaFolder,'*.set'));

if length(temp)==1
    fileName = temp.name;
    disp(fileName);
else
    disp('figure this out later'); %ask person to click on file if file doesn't exist or more than one exists
end
EEG = pop_loadset('filename',fileName,'filepath',amicaFolder);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw

%% get rid of old ICA results from infomax
EEG = rmfield(EEG,'icaweights');
EEG.icaweights = [];
EEG = rmfield(EEG,'icawinv'); 
EEG.icawinv = [];
EEG = rmfield(EEG,'icaact');
EEG.icaact = [];
EEG = rmfield(EEG,'icasphere');
EEG.icasphere = [];

%load up AMICA results
EEG.etc.amicaResultStructure = loadmodout15(amicaFolder);


%Insert AMICA results into EEGLAB
EEG.icaweights = EEG.etc.amicaResultStructure.W;
EEG.icawinv = EEG.etc.amicaResultStructure.A;

num_pcs = EEG.etc.amicaResultStructure.num_pcs; %this is to account for doing PCA reduction 
EEG.icasphere  = EEG.etc.amicaResultStructure.S(1:num_pcs,:); %if you do PCA reduction and send the full S matrix,  you will get an error. Someone online said you should take the first X rows where X = your PCA reduction

%make sure everything is ok and calculate the activation matrix
pop_editoptions('option_computeica',1); %randomly needed for some people to run even if you manually selected this checkbox with the GUI
EEG = eeg_checkset(EEG, 'ica'); %Note: I have no idea why it says "LLt not set" or what LLt is

%store updated EEG structure
EEG.comments = pop_comments(EEG.comments,'','AMICA- PCA 72',1);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET );


%update GUI
eeglab redraw
disp('Done')

   
   
   