% findSyncBuddies()
% Finds matching sync events from two systems that were send a sync/trigger
% signal. Plots sync events and asks user to use visual input to align,
% despite number of events from both systems being different
% Assumes sync signal is a square wave
%
% Inputs:
%   trigA_t          time vector of sync events in system A
%   trigB_t          time vector of sync events in system B
%   sync_freq        frequency of sync signal square wave
%   sync_amp         sync signal amplitude
%   sync_duty        sync signal duty cycle % (percentage)
%   DAQB_Fs         sampling frequency of system B
%   
% Outputs:
%   tsyncA          time vector of sync events in system A, with that has 
%                   a corresponding sync event in system B
%   tsyncB          time vector of sync events in system B, with that has 
%                   a corresponding sync event in system A
%
% Author: Noelle Jacobsen, University of Florida
% Created: 2023-Sept-10
% Last updated: 2023-Sept-12

function [tEvents_DAQa, tEvents_DAQb] = findSyncBuddies(EEG, syncB, syncB_t,sync_freq, sync_amp, sync_duty)
%check input parameters
syncA = EEG.trigger;
syncA_t = EEG.times/1000; %convert from ms to s

%find sync events
    [~,EEG_sync_lat] =findpeaks(abs(EEG.trigger), 'MinPeakHeight', .2, 'MinPeakDistance', 150 ); %Matches Rachel's way of finding events in motive syn
    syncB_events = diff(syncB)> sync_amp-sync_amp*0.05; %amplitude threshold +/- 5%
    if ~any(syncB_events) 
        syncB_events = diff(syncB)> 1; %slope for some was smaller than amplitude bc of sampling
    end
  syncB_events = diff(syncB)> 1; %sTEMP slope for some was smaller than amplitude bc of sampling


%quick quality check
 figure; tiledlayout(2,1)
    ax1=nexttile;
    plot(syncA_t, syncA); hold on;
     if ~isempty(EEG_sync_lat)
    EEG_sync_events_t = EEG.times(EEG_sync_lat)/1000; %convert from ms to sec
    stem(EEG_sync_events_t, 1.2*ones(size(EEG_sync_events_t)), 'k');
     end    
    title('EEG Sync Signal'); ylabel('Voltage')
    
  
    nexttile;
    plot(syncB_t, syncB); hold on;
    if ~isempty(syncB_events)
        syncB_events_t = syncB_t(syncB_events);
        stem(syncB_events_t, (sync_amp*1.1)*ones(size(syncB_events_t)), 'k');
    end
    title('DAQB Sync Signal'); ylabel('Voltage')

    realbad = input('Are sync signals really noisy? (1=yes, 0=no): ');
    close;

    %%
    if realbad ~=0
    % Find syncB events using template matching with square wave (0.25 Hz, 3.3 A)
    syncB_events = diff(syncB)> sync_amp-sync_amp*0.05; %amplitude threshold +/- 5%
    if ~any(syncB_events)
        syncB_events = diff(syncB)> 1; %slope for some was smaller than amplitude bc of sampling
    end
    syncB_events_t = syncB_t(syncB_events);

    template_signal = (sync_amp/2)*(1+square(syncB_t*2*pi*sync_freq,sync_duty));
    x = syncB_t; %shift template square wave to match

    %manually shift template to align with syncB
    % plot syncB with template squarewave overlaid
    fprintf('\n Shift template square wave to roughly align with raw signal. Click data tips for exact values.\n')
    figure; tiledlayout(2,1)
    ax1=nexttile;
    plot(syncB_t, syncB);
    title('DAQB Sync Signal'); ylabel('Voltage')

    ax2= nexttile; plot(x, template_signal);
    linkaxes([ax1 ax2],'x')
    title('template')
    xlabel('Time (s)'); ylabel('Voltage')
    shift = input('shift tmp signal(s): ');
    totalshift = 0;
    while shift~=0
        totalshift = totalshift+shift;
        ax2;
        plot(x+totalshift, template_signal);
        shift = input('shift tmp signal(s): ');
    end
    x = x+totalshift;
    close;
    fprintf('shift= %i\n',totalshift)

    %use template matching to get true syncB events, despite poor SNR
    tmp_syncB_events = diff(template_signal)> sync_amp-sync_amp*0.05; % threshold at sync amplitude (+/- %5)
    tmp_syncB_events_t = x(tmp_syncB_events);

    % clean up syncB events
    % loop through template events and find if true sync event exists in that
    % window
    syncB_mask = false(size(syncB_events_t)); %initialize mask
    window = 0.03; %10ms, half of time window (s) to look for sync buddies in
    %====================================================================
    %  modified from findsyncbuddies.m, Ryan Downey, UF 2021
    %====================================================================
    for eventi = 1:length(syncB_events_t) 
        SyncTime = syncB_events_t(eventi);
        [val, ~] = min(abs(tmp_syncB_events_t-SyncTime)); %find closest event in time
        if val < window % (10/1000 = +-10 ms) to look for buddies
            syncB_mask(eventi) = 1; %1 == keep event
        end
    end
    %====================================================================
    true_syncB_events_t = syncB_events_t(syncB_mask); % 1 == true event
    closeEvents = find(diff(true_syncB_events_t)<0.05);
    syncB_mask = true(size( true_syncB_events_t)); %initialize mask
    for eventi = 1:length(closeEvents)
    
        SyncTime1 = true_syncB_events_t(closeEvents(eventi));
        SyncTime2 = true_syncB_events_t(closeEvents(eventi)+1); 
        [val1, ~] = min(abs(tmp_syncB_events_t-SyncTime1)); %find closest template event in time
        [val2, ~] = min(abs(tmp_syncB_events_t-SyncTime2)); %find closest template event in time
        if val1 < val2 %check which event is closest to template event
            syncB_mask(closeEvents(eventi)+1) = 0; %remove spurious event that's farther from template event
        else
           syncB_mask(closeEvents(eventi)) = 0;
        end 
    end
   true_syncB_events_t = true_syncB_events_t(syncB_mask); % 1 == true event
 

%     %plot motive sync with template square wave
%     figure; plot(syncB_t, syncB); hold on;
%     plot(x, template_signal,'--');
%     stem(true_syncB_events_t, (sync_amp*1.1)*ones(size(true_syncB_events_t)), 'k');
%     xlabel('Time (s)'); ylabel('Voltage')
%     legend('DAQB Sync Signal','template','sync event')
    clear totalshift
    %% Clean up EEG sync events
   % disp(size(EEG_sync_events_t))
    figure; plot(syncA_t, syncA,'Color',[0 0.4470 0.7410])
    hold on; stem(EEG_sync_events_t , 1.2*ones(size(EEG_sync_events_t )), 'k');
    xlim([EEG_sync_events_t(1)-5,EEG_sync_events_t(10)]); title('EEG Sync')
    ylim([-0.5 1.5])
    remove_first_event = input('Do you want to remove the first event? 1 = YES, 0=NO      ');
    if remove_first_event == 1
        %manually remove first crappy events
        %badEvents = [1];
        %tempLogicStuff = true(1,length(EEG.event));
        %tempLogicStuff(badEvents) = 0;
        %EEG.event = EEG.event(tempLogicStuff); %keep only events not listed in "badEvents"
        EEG_sync_events_t  = EEG_sync_events_t (2:end); %simple code to remove just first one
    end
    close;
    %use template matching to get true syncB events, despite poor SNR
    template_signal = 0.5+(0.5)*square((syncA_t)*2*pi*sync_freq,sync_duty);%0.5 amplitude, 0.25 hz
    x = syncA_t; % time (s)

    % %automatically shift template to match --> can't get to work, will have to
    % %use manual
    % %shift in small increments for one full period and find where signal overal is best (min value
    % %of difference)
    % Fs = EEG.srate;
    % shift = 0:(1/sync_freq)*Fs; % # of frames to shift [0,F], where F = # frames in one period
    % sigDiff = zeros(size(shift));
    % Y1 = EEG.trigger;
    % template_signal = 0.5+(0.5)*square((EEG.times/1000+10)*2*pi*0.25);%0.5 amplitude, 0.25 hz
    % yDiff = [];
    % m = size(Y1,2);
    % for shifti = 1:length(shift)
    %     Y2 = template_signal(1+shift:end);
    %     Y2 = Y2(1:m); %crop to size;
    %     yDiff(shifti) = sum(Y1-Y2);
    % end
    % [val, shifti] = min(abs(yDiff));
    % shift_t = shift(shifti)/Fs;

    %manually shift template to match
    %plot EEG sync with template square wave
    fprintf('\n Shift template square wave to roughly align with raw signal. Click data tips for exact values.\n')
    figure; tiledlayout(2,1)
    ax1=nexttile;
    plot(syncA_t, syncA);
    ylim([-0.5 1.5])
    title('EEG Sync Signal')

    ax2= nexttile; plot(x, template_signal); hold off;
    ylim([-0.5 1.5])
    title('Template Square Wave')
    linkaxes([ax1 ax2],'x')
    shift = input('shift tmp signal(s): ');
    totalshift = 0;
    while shift~=0
        totalshift = totalshift+shift;
        ax2;
        plot(x+totalshift, template_signal);
        shift = input('shift tmp signal(s): ');
    end
    x = x+totalshift;
    fprintf('shift= %i\n',totalshift)
    close;

    %find template sync events
    tmp_syncAevents = diff(template_signal)> 0.9; % threshold at sync amplitude (+/- %5)
    tmp_syncAevents_t = x(tmp_syncAevents);

    % loop through templatate events and find if true sync event exists in that
    % window
    syncA_mask = false(size(EEG_sync_events_t ));
    window = 0.03; %30ms, half of time window (s) to look for sync buddies in
    %====================================================================
    %  modified from findsyncbuddies.m, Ryan Downey, UF 2021
    %====================================================================
    for tmpeventi = 1:length(tmp_syncAevents_t)
        tempSyncTime = tmp_syncAevents_t(tmpeventi);
        [val, eventi] = min(abs(EEG_sync_events_t -tempSyncTime)); %find closest event in time
        if val < window % (10/1000 = +-10 ms) to look for buddies
            syncA_mask(eventi) = true;
        end
    end
    %====================================================================
    true_syncA_events_t = EEG_sync_events_t (syncA_mask);
    
    %remove spurious events too close together
    closeEvents = find(diff(true_syncA_events_t)<0.05);
    syncA_mask = true(size(true_syncA_events_t)); %initialize mask
    for eventi = 1:length(closeEvents)
        SyncTime1 = true_syncA_events_t(closeEvents(eventi));
        SyncTime2 = true_syncA_events_t(closeEvents(eventi)+1); 
        [val1, ~] = min(abs(tmp_syncA_events_t-SyncTime1)); %find closest template event in time
        [val2, ~] = min(abs(tmp_syncA_events_t-SyncTime2)); %find closest template event in time
        if val1 < val2 %check which event is closest to template event
            syncA_mask(closeEvents(eventi)+1) = 0; %remove spurious event that's farther from template event
        else
           syncA_mask(closeEvents(eventi)) = 0;
        end 
    end
   true_syncA_events_t = true_syncA_events_t(syncA_mask); % 1 == true event
    disp(['Average EEG sync pulse time difference (s) = ',num2str(mean(diff(true_syncA_events_t)))]);

%     figure; plot(syncA_t, syncA); hold on;
%     stem(true_syncA_events_t, 1.2*ones(size(true_syncA_events_t)), 'k');
%     legend('sync','events')
%     title('All EEG Sync Events');
%     xlabel('Time (s)'); ylabel('Digital Signal');
    clear tmp_eventi eventi
    else
        true_syncA_events_t= EEG_sync_events_t;
        true_syncB_events_t = syncB_events_t;
    end

    %% Roughly align sync A and B
    % visually align to get a rough estimate of sync buddies
    
    fprintf('\nRougly align DAQB sync to EEG sync. Click data tips for exact datapoints \n')
    figure; tiledlayout(2,1)
    ax1=nexttile;
    plot(syncA_t, syncA); title('EEG')
    hold on; stem(true_syncA_events_t, 1.2*ones(size(true_syncA_events_t)), 'k');
    legend('sync','events')

    totalshift = true_syncA_events_t(1)-true_syncB_events_t(1); %estimate shift of DAQB time
    ax2 = nexttile; plot(syncB_t+totalshift, syncB); hold on;
    stem(true_syncB_events_t+totalshift, 3.35*ones(size(true_syncB_events_t)), 'k'); hold off;
    legend('sync','events'); title('DAQB')
    linkaxes([ax1 ax2],'x')
    shift = input('shift tmp signal(s): ');
    while shift~=0
        totalshift = totalshift+shift;
        ax2;
        plot(syncB_t+totalshift, syncB); hold on;
        stem(true_syncB_events_t+totalshift, 3.35*ones(size(true_syncB_events_t)), 'k'); hold off;
        shift = input('shift tmp signal(s): ');
    end
    fprintf('Shift (s)= %i\n',totalshift)
    close;
    %% Find sync buddies
    % loop through syncA events and find a matching syncB event in
    % specified time window
    syncA_mask1 = false(size(true_syncA_events_t));
    syncB_mask1 = false(size(true_syncB_events_t));
    true_syncB_events_tshift = true_syncB_events_t+totalshift; %shift time of syncB events to roughly align with system A, just for visualization purposes
    window = 0.01; %10ms, half of time window (s) to look for sync buddies in
    %====================================================================
    %  modified from findsyncbuddies.m, Ryan Downey, UF 2021
    %====================================================================
    for eventiA = 1:length(true_syncA_events_t)
        tempSyncTime = true_syncA_events_t(eventiA);
        [val, eventiB] = min(abs(true_syncB_events_tshift-tempSyncTime)); %find closest event in time
        if val < window % (10/1000 = +-10 ms) to look for buddies
            syncB_mask1(eventiB) = 1;
            syncA_mask1(eventiA) = 1;
        end
    end

    % loop through syncB events and find a matching syncA event
    syncA_mask2 = false(size(true_syncA_events_t));
    syncB_mask2 = false(size(true_syncB_events_t));
    for eventiB = 1:length(true_syncB_events_tshift)
        tempSyncTime = true_syncB_events_tshift(eventiB);
        [val eventiA] = min(abs(true_syncA_events_t-tempSyncTime)); %find closest event in time
        if val < window % (10/1000 = +-10 ms) to look for buddies
            syncA_mask2(eventiA) = 1;
            syncB_mask2(eventiB) = 1;
        end
    end

    %find matching events found in both systems
    syncA_mask = syncA_mask1 & syncA_mask2; %final mask
    syncB_mask = syncB_mask1 & syncB_mask2;
    tEvents_DAQa  = true_syncA_events_t(syncA_mask);
    tEvents_DAQb  = true_syncB_events_t(syncB_mask);
    
    %plot final sync matches
    figure; tiledlayout(2,1)
    ax1=nexttile;
    plot(syncA_t, syncA);
    hold on; stem(tEvents_DAQa, 1.2*ones(size(tEvents_DAQa)), 'k');
    ylim([-0.5 1.5]);
    xlabel('Time (s)'); ylabel('Digital Signal'); title('EEG Sync')
    legend({'signal','sync events'},'Location','southeast')

    ax2 = nexttile; plot(syncB_t+totalshift, syncB); hold on;
    stem(tEvents_DAQb+totalshift, 4*ones(size(tEvents_DAQb)), 'k'); hold off;
    linkaxes([ax1 ax2],'x')
    ylim([-1 5]);
    xlabel('Time (s)'); ylabel('Voltage'); title('Motive Sync')
    legend({'signal','sync events'},'Location','southeast')
end