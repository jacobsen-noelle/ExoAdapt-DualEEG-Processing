% Add channel locations
% Inputs:
% EEG               EEG structure
% eleclocsFolder    folder path that contains subject text file with
%                   channel locations in XYZ coordinate format
%                   labeled through EEGLAB GUI
%
% Note: this function searches for text file that contains "SXX_chanlocs.txt"


function EEG = add_chanlocs(EEG, eleclocsFolder)
%% Remove excess channels (only keep 128 scalp, 8 EMG, and 128 noise)
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG,'channel',{'1-A1' '1-A2' '1-A3' '1-A4' '1-A5' '1-A6' '1-A7' '1-A8' '1-A9' '1-A10' '1-A11' '1-A12' '1-A13' '1-A14' '1-A15' '1-A16' '1-A17' '1-A18' '1-A19' '1-A20' '1-A21' '1-A22' '1-A23' '1-A24' '1-A25' '1-A26' '1-A27' '1-A28' '1-A29' '1-A30' '1-A31' '1-A32' '1-B1' '1-B2' '1-B3' '1-B4' '1-B5' '1-B6' '1-B7' '1-B8' '1-B9' '1-B10' '1-B11' '1-B12' '1-B13' '1-B14' '1-B15' '1-B16' '1-B17' '1-B18' '1-B19' '1-B20' '1-B21' '1-B22' '1-B23' '1-B24' '1-B25' '1-B26' '1-B27' '1-B28' '1-B29' '1-B30' '1-B31' '1-B32' '1-C1' '1-C2' '1-C3' '1-C4' '1-C5' '1-C6' '1-C7' '1-C8' '1-C9' '1-C10' '1-C11' '1-C12' '1-C13' '1-C14' '1-C15' '1-C16' '1-C17' '1-C18' '1-C19' '1-C20' '1-C21' '1-C22' '1-C23' '1-C24' '1-C25' '1-C26' '1-C27' '1-C28' '1-C29' '1-C30' '1-C31' '1-C32' '1-D1' '1-D2' '1-D3' '1-D4' '1-D5' '1-D6' '1-D7' '1-D8' '1-D9' '1-D10' '1-D11' '1-D12' '1-D13' '1-D14' '1-D15' '1-D16' '1-D17' '1-D18' '1-D19' '1-D20' '1-D21' '1-D22' '1-D23' '1-D24' '1-D25' '1-D26' '1-D27' '1-D28' '1-D29' '1-D30' '1-D31' '1-D32' '1-EXG1' '1-EXG2' '1-EXG3' '1-EXG4' '1-EXG5' '1-EXG6' '1-EXG7' '1-EXG8' '2-A1' '2-A2' '2-A3' '2-A4' '2-A5' '2-A6' '2-A7' '2-A8' '2-A9' '2-A10' '2-A11' '2-A12' '2-A13' '2-A14' '2-A15' '2-A16' '2-A17' '2-A18' '2-A19' '2-A20' '2-A21' '2-A22' '2-A23' '2-A24' '2-A25' '2-A26' '2-A27' '2-A28' '2-A29' '2-A30' '2-A31' '2-A32' '2-B1' '2-B2' '2-B3' '2-B4' '2-B5' '2-B6' '2-B7' '2-B8' '2-B9' '2-B10' '2-B11' '2-B12' '2-B13' '2-B14' '2-B15' '2-B16' '2-B17' '2-B18' '2-B19' '2-B20' '2-B21' '2-B22' '2-B23' '2-B24' '2-B25' '2-B26' '2-B27' '2-B28' '2-B29' '2-B30' '2-B31' '2-B32' '2-C1' '2-C2' '2-C3' '2-C4' '2-C5' '2-C6' '2-C7' '2-C8' '2-C9' '2-C10' '2-C11' '2-C12' '2-C13' '2-C14' '2-C15' '2-C16' '2-C17' '2-C18' '2-C19' '2-C20' '2-C21' '2-C22' '2-C23' '2-C24' '2-C25' '2-C26' '2-C27' '2-C28' '2-C29' '2-C30' '2-C31' '2-C32' '2-D1' '2-D2' '2-D3' '2-D4' '2-D5' '2-D6' '2-D7' '2-D8' '2-D9' '2-D10' '2-D11' '2-D12' '2-D13' '2-D14' '2-D15' '2-D16' '2-D17' '2-D18' '2-D19' '2-D20' '2-D21' '2-D22' '2-D23' '2-D24' '2-D25' '2-D26' '2-D27' '2-D28' '2-D29' '2-D30' '2-D31' '2-D32'});

%% Add channel types
for i= 1:length(EEG.chanlocs)
    if contains({EEG.chanlocs(i).labels}, 'EMG')
        EEG.chanlocs(i).type = 'EMG';
    elseif contains({EEG.chanlocs(i).labels}, 'EXG')
        EEG.chanlocs(i).type = 'EMG';
    elseif contains({ EEG.chanlocs(i).labels},'2-')
        EEG.chanlocs(i).type = 'Noise';
    elseif contains({ EEG.chanlocs(i).labels},'1-') && ~contains({EEG.chanlocs(i).labels}, 'EXG')
        EEG.chanlocs(i).type = 'EEG';
    end
end

getchantypes; %
[non_EEG_chans] = [setdiff(1:length(EEG.chanlocs),EEG_chans)];
%% Read channel locations from text file
eleclocsFolder= char( eleclocsFolder);
fileList = dir([eleclocsFolder,'\*.txt']);
filei = find(contains({fileList.name},[EEG.subject,'_chanlocs.txt']));
filename = fileList(filei).name;
cd(eleclocsFolder);
old_locs = EEG.chanlocs;
fprintf('Importing locations with readlocs()...\n')
mychanlocs = readlocs(filename,'format',{'labels','X','Y','Z'}); %load channel locations from txt file

chan_ind=[];
for index = 1:length(EEG.chanlocs)
    c = EEG.chanlocs(index).labels;
    if ~isempty(find(strcmpi(c,{mychanlocs.labels})))
        x = find(strcmpi({mychanlocs.labels},c));
        chan_ind = [chan_ind; x];
    end
end
EEG.chanlocs = mychanlocs(chan_ind); %add EEG chan locs

for index = 1:length(non_EEG_chans) %add other channel labels that don't have locations
    EEG.chanlocs(end+1).labels = old_locs(non_EEG_chans(index)).labels;
end
EEG = eeg_checkchanlocs(EEG);

%% Optimize head center
EEG.chanlocs = pop_chancenter( EEG.chanlocs, [],[non_EEG_chans] ); %omit EMG/Noise channels from optimization

%% Verfiy that electrode locations look normal
%avg head length = 187mm (M),197mm (F)
%avg head breadth = 152mm (M), 144mm (F)
maxY = max([EEG.chanlocs.Y]);
minY = min([EEG.chanlocs.Y]);
maxX = max([EEG.chanlocs.X]);
minX = min([EEG.chanlocs.X]);
diffY = maxY-minY; %mm
diffX = maxX-minX; %mm
if diffX<172 | diffX>213 | diffY< 133 | diffY >165
    warning('Channel locations seem abnormal')
end

%% Add channel types

end

