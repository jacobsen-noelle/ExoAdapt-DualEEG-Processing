%% Fix Biosemi events labeled by pop_biosig()
% The pop_biosig plug-in can mislabel events when importing BDFs from Biosemi. This function fixes event labels
% by reading the decimal data output by pop_biosig 24bit data from Biosemi USB trigger interface
%
% EEG = fix_popsig_eventlabels (EEG) - fix events mislabeled by BIOSIG
% toolbox and update EEG.event structure
%
% Inputs:
%   EEG - EEGLAB data structure
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Noelle Jacobsen, University of Florida
% Created 18-March-2021

function [EEG] = fix_popsig_eventlabels (EEG)
%Covert event channel data from decimal to 24bit
channels = de2bi(EEG.H.BDF.ANNONS);%Convert decimal numbers in event channel to binary numbers
my_chans_ind =[];
count = 1;
%each bit is a separate digital channel. Look through each channel to find the
%channels that have event information based on them switching states
%between 0-1. If a channel has at least one state change (diff ~=0), then
%that channel index number will be recorded in the variable "my_chans_ind"
for i= 1:24
    events =  find(diff(channels(:,i))~=0);
    if size(events,1)~=0 %Check that there is at least one event (i.e. digital state change)
        my_chans_ind(count) = i;
        count= count+1;
    end
end
clear events

%Find events within the digital channel
for i=1:size(my_chans_ind,2)
    event_chan = channels(:,my_chans_ind(i)); %analog data from event channel
    good_events= find(diff(channels(:,my_chans_ind(i)))~=0); %identify state change (LOW/HIGH, 0/1)
    %Sort channels based on expected number of events (handful of button
    %presses, >100 trigger/sync signal events, and maybe one random event
    %at the beginning or end
    if size(good_events,1) >1 &&  size(good_events,1) <15
        button_event_lat = find(diff(channels(:,my_chans_ind(i)))==-1)+1; %find button press events, searching for '-1' because that's when button state is switching from HIGH->LOW  (1->0), default state is HIGH (1).  Adding 1 to index because we want the position where the state turns LOW (button pressed), not the latency of the last HIGH state
        %Add events to EEG events structure
        for x = 1:size(button_event_lat,1)
            EEG.event(end+1).latency = button_event_lat(x); %add new trigger events
            EEG.event(end).type= 'BUTTON';
        end
        fprintf('Found %i button press events\n', size(button_event_lat,1));
        if size(button_event_lat,1)>50
            error('Number of button press events seems abnormal. Check decoding port data from 24bit');
        end
    elseif size(good_events,1) >=15 %check for at least 1 min of data (0.25 Hz sync); may need to adjust this threshold-- sometimes there are garbage channels random events
  EEG.trigger = event_chan'; %add trigger data to EEG structure, transposed to match EEGlab format [channel x time];
        %Find trigger events. Note: The 6th bit codes for the port connected to signal generator
        trig_event_lat = find(diff(channels(:,my_chans_ind(i)))==1)+1; %find rising edge of square wave, searching for where diff == 1 because that means the port is switching from 0->1 (default state LOW?). Adding 1 to index because we want the position where the state turns HIGH, not the latency of the last LOW state
        %figure; plot(EEG.trigger'); hold on; stem(trig_event_lat,1.5*ones(size(trig_event_lat)));
        %Add events to EEG events structure
        for x = 1:size(trig_event_lat,1)
            EEG.event(end+1).latency = trig_event_lat(x); %add new trigger events
            EEG.event(end).type= 'SYNC';
        end
        fprintf('Found %i trigger events\n', size(trig_event_lat,1));
        if size(trig_event_lat,1)<5
            error('Number of trigger events seems abnormal. Check decoding port data from 24bit');
        end
    elseif size(good_events,1)==1
        other_event_lat= find(diff(channels(:,my_chans_ind(i)))~=0)+1;
        %Add random events to EEG events structure
        for x = 1:size(other_event_lat,1)
            EEG.event(end+1).latency = other_event_lat(x); %add new events
            eventName = num2str(EEG.H.BDF.ANNONS(other_event_lat(x)));
            EEG.event(end).type= eventName;
        end
        fprintf('Found %i other event(s)\n', size(other_event_lat,1));
    end
end

clear event_chan good_events my_chans_ind button_event_pos other_event_pos
end