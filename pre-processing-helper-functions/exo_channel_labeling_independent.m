%% Electrode Channel Labeling
% Noelle Jacobsen, 4/13/21
% Manually label electrode locations and save locations as file called "SXX_chanlocs.txt""SXX_chanlocs.mat"in
% specified output folder 

% Close windows, clear workspace + command window
clc; close all; clear all;

%enter Z drive directory name
mydir= 'R:\Ferris-Lab\'; 
inputFolder = strcat(mydir,'jacobsen.noelle\EEG_Processing\Channel Locations');
addpath(inputFolder)
% outputFolder = input('Enter name of output folder:');
% or manually code output folder
outputFolder = strcat(mydir,'johnprieschl\Exo Electrode Locations\Corrected Chanlocs'); 
eeglabPath = 'R:\Ferris-Lab\johnprieschl\EEG Processing\EEGLab Plugins\eeglab2021.1';
addpath(eeglabPath)
rmpath(genpath('R:\F25erris-Lab\share\MindInMotion\eeglab2020_0\'))

%% Enter subject number and file paths

subject_num = input('Enter Subject number: S');
PHD_spec = subject_num + 14;
PHD_spec = num2str(PHD_spec);
subject_num= num2str(subject_num);
if (size(subject_num,2)==1)
    subject_num= strcat('0',subject_num);
end
s = input('Enter scanner number (1= Artec, 2= itSeez3d, 3= Structure Scanner): ');
if s == 1
   scanner = 'Artec'
elseif s==2
   scanner = 'itSeez3D'
elseif s==3
    scanner = 'Structure Scanner'
else
    error('Scanner name error: Please enter a number (1,2,3) that corresponds with the scanner used for this collection');
end 

%% Open electrode template image

cd(inputFolder)
%template = imread('cap_layout_128.jpg');
%image(template);

%% Load EEG file based on subject number

EEGfile = strcat('S',subject_num,'_stand_pre.set'); %only need to mark locations for one file per subject, as electrode locations are reused for the same subject
EEGFolder = strcat(mydir,'johnprieschl\Exo Electrode Locations\Data\','PHD',PHD_spec,'_S',subject_num,'\EEG'); %folder that contains the EEG dataset
scanPath = strcat(mydir,'johnprieschl\Exo Electrode Locations\Data\PHD',PHD_spec,'_S',subject_num,'\Head Scan'); %folder that contains the final scan with .obj, .jpg, and .mtl file
cd(EEGFolder);

%  Open EEGlab (if not already open) and load EEG files
if ~exist('ALLCOM')
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
end

EEG = pop_loadset('filename',EEGfile,'filepath',EEGFolder);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );

% remove all non-scalp channels 
EEG = pop_select(EEG,'channel',{'1-A1' '1-A2' '1-A3' '1-A4' '1-A5' '1-A6' '1-A7' '1-A8' '1-A9' '1-A10' '1-A11' '1-A12' '1-A13' '1-A14' '1-A15' '1-A16' '1-A17' '1-A18' '1-A19' '1-A20' '1-A21' '1-A22' '1-A23' '1-A24' '1-A25' '1-A26' '1-A27' '1-A28' '1-A29' '1-A30' '1-A31' '1-A32' '1-B1' '1-B2' '1-B3' '1-B4' '1-B5' '1-B6' '1-B7' '1-B8' '1-B9' '1-B10' '1-B11' '1-B12' '1-B13' '1-B14' '1-B15' '1-B16' '1-B17' '1-B18' '1-B19' '1-B20' '1-B21' '1-B22' '1-B23' '1-B24' '1-B25' '1-B26' '1-B27' '1-B28' '1-B29' '1-B30' '1-B31' '1-B32' '1-C1' '1-C2' '1-C3' '1-C4' '1-C5' '1-C6' '1-C7' '1-C8' '1-C9' '1-C10' '1-C11' '1-C12' '1-C13' '1-C14' '1-C15' '1-C16' '1-C17' '1-C18' '1-C19' '1-C20' '1-C21' '1-C22' '1-C23' '1-C24' '1-C25' '1-C26' '1-C27' '1-C28' '1-C29' '1-C30' '1-C31' '1-C32' '1-D1' '1-D2' '1-D3' '1-D4' '1-D5' '1-D6' '1-D7' '1-D8' '1-D9' '1-D10' '1-D11' '1-D12' '1-D13' '1-D14' '1-D15' '1-D16' '1-D17' '1-D18' '1-D19' '1-D20' '1-D21' '1-D22' '1-D23' '1-D24' '1-D25' '1-D26' '1-D27' '1-D28' '1-D29' '1-D30' '1-D31' '1-D32'});

% use getchanlocs to mark channel locations on 3D head image
EEG = get_chanlocs_mod(EEG, scanPath, 'deleteTxtOutput', 0,'grayTextures',0, 'moveElecInwards', 7.5,'templatePath',strcat(inputFolder,'\montageTemplate.mat'), 'saveName',strcat(outputFolder,'\PHD',PHD_spec,'S',subject_num,'_chanlocs'),'scannerAppName', scanner); %'templatePath',strcat(inputFolder,'\montageTemplate'));
cd(outputFolder);
file = strcat('PHD',PHD_spec,'_S',subject_num,'_chanlocs');
scan.chanlocs = EEG.chanlocs(1:128);
scan.chaninfo = EEG.chaninfo;
save(file, 'scan')

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
eeglab redraw

%plot locations
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];

figure; 
plot3(X,Y,Z,'o')
title(['S',subject_num,' channel locations'])
saveas(gcf,file,'jpg')
saveas(gcf,file,'fig')
close;
%% Look back at the output folder and its files...

fprintf('\n\n\tRemember to rename the files (correct underscores and spacing!)\n\n')
