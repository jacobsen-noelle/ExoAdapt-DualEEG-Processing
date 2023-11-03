function gaitevents = getgaitevents(EEG, thresh, outputfolder)
condition = extractAfter(EEG.filename,[EEG.subject,'_']);
condition = extractBefore(condition,'.set');
GRF_L = EEG.data(strcmp({EEG.chanlocs.labels},'Left_GRF'),:);
GRF_R = EEG.data(strcmp({EEG.chanlocs.labels},'Right_GRF'),:);

if size(GRF_R,1)>1
    GRF_R = GRF_R';
end
if size(GRF_L,1)>1
    GRF_L = GRF_L';
end

% Left forceplate gait events
% find 40 consecutive frames above 10N
idx= GRF_L> thresh;
ii1=strfind([0 idx 0],[0 1]);
ii2=strfind([0 idx 0],[1 0])-1;
ii=(ii2-ii1+1)>=40; % at least 40 consecutive frames
%     data_out=arrayfun(@(x,y) GRF_L(x:y),ii1(ii),ii2(ii),'un',0);
%     time_out=arrayfun(@(x,y) Time(x:y),ii1(ii),ii2(ii),'un',0);
LHS = ii1(ii); %ii1 =start index of threshold crossing
LTO = ii2(ii); %ii2 =start index of threshold crossing


% Right forceplate gait events
% find 40 consecutive frames above 10N
idx= GRF_R> thresh;
ii1=strfind([0 idx 0],[0 1]);
ii2=strfind([0 idx 0],[1 0])-1;
ii=(ii2-ii1+1)>=40; % at least 40 consecutive frames
RHS = ii1(ii);
RTO = ii2(ii);

%second thresh- backup estimate events
thresh2 = 300;
% Left forceplate gait events
% find 40 consecutive frames above 10N
idx= GRF_L> thresh2;
ii1=strfind([0 idx 0],[0 1]);
ii2=strfind([0 idx 0],[1 0])-1;
ii=(ii2-ii1+1)>=40; % at least 40 consecutive frames
%     data_out=arrayfun(@(x,y) GRF_L(x:y),ii1(ii),ii2(ii),'un',0);
%     time_out=arrayfun(@(x,y) Time(x:y),ii1(ii),ii2(ii),'un',0);
LHS_tmp = ii1(ii); %ii1 =start index of threshold crossing
LTO_tmp = ii2(ii); %ii2 =start index of threshold crossing


% Right forceplate gait events
% find 40 consecutive frames above 10N
idx= GRF_R> thresh2;
ii1=strfind([0 idx 0],[0 1]);
ii2=strfind([0 idx 0],[1 0])-1;
ii=(ii2-ii1+1)>=40; % at least 40 consecutive frames
RHS_tmp = ii1(ii);
RTO_tmp = ii2(ii);
% Estimate gait event  from abnormal GFR using higher threshold and linear
% fit
x2 = 50; %num x points to calulate slope from
y = -10; %y thresh crossing we're estimating with slope. -10 was more accurate than finding 0 or 10 N crossing
m = (GRF_L(LHS_tmp+x2)- GRF_L(LHS_tmp))/x2; %slope
LHS_est = LHS_tmp-( (GRF_L(LHS_tmp)-y)/m );
m = (GRF_L(LTO_tmp-x2)- GRF_L(LTO_tmp))/x2; %slope
LTO_est = LTO_tmp-( (y-GRF_L(LTO_tmp))/m );

m = (GRF_R(RHS_tmp+x2)- GRF_R(RHS_tmp))/x2; %slope
RHS_est = RHS_tmp-( (GRF_R(RHS_tmp)-y)/m );
m = (GRF_R(RTO_tmp-x2)- GRF_R(RTO_tmp))/x2; %slope
RTO_est = RTO_tmp-( (y-GRF_R(RTO_tmp))/m );


% only keep estimated gait events where there isn't already a marked event
% in the window
window = 0.25*EEG.srate; % 0.25s window, units are in frames
events_est = {'LHS_est','LTO_est','RHS_est','RTO_est'};
events = {'LHS','LTO','RHS','RTO'};
newevents = {};
%%
for typei = 1:length(events_est)
    est = eval(events_est{typei});
    GE = eval(events{typei});
    est_event_mask = false(size(est));
    for eventi = 1:length(est)
        EventTime = est(eventi);
        [val, ~] = min(abs(GE-EventTime)); %find closest event in time
        if val > window % if closest event is outside of this window, keep the estimated events
            est_event_mask(eventi) = 1; %1 == keep event
        end
    end
    newevents{typei} = est(est_event_mask);
end
LHS = sort([LHS, newevents{1,1}])';
LTO = sort([LTO,newevents{1,2}])';
RHS = sort([RHS, newevents{1,3}])';
RTO = sort([RTO, newevents{1,4}])';
%%

% figure; tiledlayout(2,1);
% ax1= nexttile;
% X = 1:EEG.pnts;
% plot( GRF_L,'Color',[0.7 0.7 0.7]); hold on;
% stem(LHS, mean(GRF_L)*ones(size(LHS)), 'c');
% stem(LTO,mean(GRF_L)*ones(size(LTO)), 'm')
% legend('forceplate','LHS','LTO')diarydiary
% ylabel('Ground Reaction Force (N)')
% title('Left')
% if strcmp(xax,'time')
%     xlabel('Time (s)');
%     fprintf('Event Type \t Average Time Diff (s)\n')
% else
%     xlabel('Latency')
% end
%
% ax2 = nexttile;
% plot(X, GRF_R,'Color',[0.7 0.7 0.7]); hold on;
% stem(RHS, mean(GRF_R)*ones(size(RHS)), 'c');
% stem(RTO,mean(GRF_R)*ones(size(RTO)), 'm');
% linkaxes([ax1, ax2],'x')
% legend('forceplate','RHS','RTO')
% ylabel('Ground Reaction Force (N)')
% title('Right')
condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO); %events x subcondition
%% Check gait events

if ~strcmp(EEG.condition,'pow_2') %don't need to double check this condition
    condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO); %events x subcondition
    for condi = 1:size(condmask.RHS,2)

        cursor_info = [];
        % left
        [newHSLocs,newTOLocs,remove_HSevent,remove_TOevent] = fixGaitEvents (LHS(condmask.LHS(:,condi)),...
            LTO(condmask.LTO(:,condi))',GRF_L); %need 30 strides for sub, 40 just for padding in cause some are rejected in preprocessing
        if ~isempty(remove_HSevent)
            for v=1:length(remove_HSevent)
                LHS(LHS == remove_HSevent(v)) = []; %remove events at this index (cursor locations user selected for manual event removal)
            end
        end
        if ~isempty(remove_TOevent)
            for v = 1:length(remove_TOevent)
                LTO(LTO == remove_TOevent(v)) = []; %remove event at this cursor location index
            end
        end
        LHS = sort([LHS;newHSLocs]); LTO = sort([LTO; newTOLocs]);

        condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO);
        %double check all events (manually add/remove events)
        fprintf(['\nLeft GRF, pt %i\n',condi])
        tmp = LHS(condmask.LHS(:,condi));
        twin = tmp(1):tmp(end);
        [LHSsubcon, LTOsubcon] = check_gaitevents(LHS(condmask.LHS(:,condi)),LTO(condmask.LTO(:,condi)),GRF_L); %#ok<ASGLU>
        LHS(condmask.LHS(:,condi)) =[];
        LTO(condmask.LTO(:,condi)) =[];
        LHS = sort([LHS;LHSsubcon]);
        LTO = sort([LTO;LTOsubcon]);



        %% right
        [newHSLocs,newTOLocs,remove_HSevent,remove_TOevent] = fixGaitEvents (RHS(condmask.RHS(:,condi)),...
            RTO(condmask.RTO(:,condi))',GRF_R);
        if ~isempty(remove_HSevent)
            for v=1:length(remove_HSevent)
                RHS(RHS == remove_HSevent(v)) = []; %remove events at this index (cursor locations user selected for manual event removal)
            end
        end
        if ~isempty(remove_TOevent)
            for v = 1:length(remove_TOevent)
                RTO(RTO == remove_TOevent(v)) = []; %remove event at this cursor location index
            end
        end

        RHS = sort([RHS; newHSLocs]); RTO = sort([RTO;newTOLocs]);

        condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO);
        %double check all events (manually add/remove events)
        fprintf(['\nRight GRF, pt %i\n',condi])
        tmp = RHS(condmask.RHS(:,condi));
        twin = tmp(1):tmp(end);
        [RHSsubcon, RTOsubcon] = check_gaitevents(RHS(condmask.RHS(:,condi)),RTO(condmask.RTO(:,condi)),GRF_R);
        RHS(condmask.RHS(:,condi)) =[];
        RTO(condmask.RTO(:,condi)) =[];
        RHS = sort([RHS;RHSsubcon]);
        RTO = sort([RTO;RTOsubcon]);

        %% check gait cycle sequence
       [RHS, LHS, RTO, LTO] = check_gaitcycle(RHS,LHS,RTO,LTO,EEG.srate);
 
        %% confirm events
        condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO); %events x subcondition
        tmp = RHS(condmask.RHS(:,condi));
        twin = round(tmp(1)):round(tmp(end));
        figure; tiledlayout(2,1);
        ax1= nexttile;
        plot(EEG.times(twin)/1000, GRF_L(twin),'Color',[0.7, 0.7,0.7]); hold on;
        stem(LHS(condmask.LHS(:,condi))/EEG.srate, mean(GRF_L)*ones(size(LHS(condmask.LHS(:,condi)))), 'c');
        stem(LTO(condmask.LTO(:,condi))/EEG.srate,mean(GRF_L)*ones(size(LTO(condmask.LTO(:,condi)))), 'm')
        legend('forceplate','LHS','LTO')
        xlabel('Time (s)'); ylabel('Ground Reaction Force (N)')
        title('Left')

        ax2 = nexttile;
        plot(EEG.times(twin)/1000, GRF_R(twin),'Color',[0.7, 0.7,0.7]); hold on;
        stem(RHS(condmask.RHS(:,condi))/EEG.srate, mean(GRF_R)*ones(size(RHS(condmask.RHS(:,condi)))), 'c');
        stem(RTO(condmask.RTO(:,condi))/EEG.srate,mean(GRF_R)*ones(size(RTO(condmask.RHS(:,condi)))), 'm');
        linkaxes([ax1, ax2],'yx')
        legend('forceplate','RHS','RTO')
        xlabel('Time (s)'); ylabel('Ground Reaction Force (N)')
        title('Right')
        xlim([EEG.times(twin([1,end]))/1000])

        if condi ==2
            sgtitle([EEG.subject,', ',condition],'interpreter','none')
            figtitle = [EEG.subject,'_',condition,'_pt2'];
        else
            sgtitle([EEG.subject,', ',condition],'interpreter','none')
            figtitle = [EEG.subject,'_',condition];
        end
        %%
        savethisfig(gcf,figtitle,[outputfolder,'\gaitevents_suncond\jpg'],'jpg');
        savethisfig(gcf,figtitle,[outputfolder,'\gaitevents_subcond\fig'],'fig');
        close;

        condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO);

    end

end


%% compile events intro structure
gaitevents = [];
eventtype = {'LHS','LTO','RHS','RTO'};
for i = 1:length(eventtype)
    dat = eval(eventtype{i});
    for eventi = 1:length(dat)
        gaitevents(end+1).type = eventtype{i};
        gaitevents(end).latency = dat(eventi);
    end
end
close;

end



function  [RHS_good, LHS_good, RTO_good, LTO_good] = check_gaitcycle(RHS,LHS,RTO,LTO, Fs)
RHS_good = []; LHS_good = []; RTO_good = []; LTO_good = [];
for ge = 1:size(RHS,1)
    wind_start = RHS(ge);
    if ge==size(RHS,1) %if last gait event, estimate window
        wind_end = RHS(ge)+2*Fs; %estimate next gait event (2s gait cycle)
    else
        wind_end = RHS(ge+1); %use next RHS event as window boundary
    end
    if ((wind_end-wind_start)/Fs)<2.5 %if window is greater than 2.5s, bad gait cycle, skip to next
        LHS_subi = find(LHS > wind_start & LHS< wind_end);
        LTO_subi = find(LTO > wind_start & LTO< wind_end);
        RTO_subi = find(RTO > wind_start & RTO< wind_end);
        if ~isempty(LHS_subi) && ~isempty(LTO_subi) && ~isempty(RTO_subi)
            RHS_good = [RHS_good; RHS(ge)];
            LHS_good = [LHS_good; LHS(LHS_subi(1))];
            LTO_good = [LTO_good; LTO(LTO_subi(1))];
            RTO_good = [RTO_good; RTO(RTO_subi(1))];
        end
    end
end

end

function condmask = getsubcondmask(EEG,RHS,RTO,LHS,LTO)
condmask.RHS = false(size(RHS));%theres probably a more efficient way to do this but I can't think rn
condmask.LHS = false(size(LHS));
condmask.RTO = false(size(RTO));
condmask.LTO = false(size(LTO));
if strcmp(EEG.condition ,'pow_1')
    condmask.RHS(1:30) = 1;
    condmask.LHS(1:30) = 1;
    condmask.RTO(1:30) = 1;
    condmask.LTO(1:30) = 1;
elseif strcmp(EEG.condition ,'deadapt')
    condmask.RHS = false([length(RHS),2]);
    condmask.RHS(1:30,1) = 1; %beginning strides
    condmask.RHS(end-30:end,2) = 1;
    condmask.LHS = false([length(LHS),2]); %cycle twice, for beginning and end of deadapt
    condmask.LHS(1:30,1) = 1; %beginning strides
    condmask.LHS(end-30:end,2) = 1;
    condmask.RTO = false([length(RTO),2]); %cycle twice, for beginning and end of deadapt
    condmask.RTO(1:30,1) = 1; %beginning strides
    condmask.RTO(end-30:end,2) = 1;
    condmask.LTO = false([length(LTO),2]); %cycle twice, for beginning and end of deadapt
    condmask.LTO(1:30,1) = 1; %beginning strides
    condmask.LTO(end-30:end,2) = 1;
else
    try
        condmask.RHS(end-30:end,:) =1;
        condmask.LHS(end-30:end,:) =1;
        condmask.RTO(end-30:end,:) =1;
        condmask.LTO(end-30:end,:) =1;
    catch
        condmask.RHS = true(size(RHS));%theres probably a more efficient way to do this but I can't think rn
        condmask.LHS = true(size(LHS));
        condmask.RTO = true(size(RTO));
        condmask.LTO = true(size(LTO));
    end

end
end


function [newHSLocs,newTOLocs,remove_HSevent,remove_TOevent] = fixGaitEvents (HeelStrikeLocs,ToeOffLocs,GRF)
cont = 'No, recheck for bad gait events';
% HeelStrikeLats = (HeelStrikeLocs*EEG.srate)'; %latency of events
% ToeOffLats = (ToeOffLocs*EEG.srate)';
newHSLocs = [];
newTOLocs =[];
remove_HSevent =[];
remove_TOevent=[];
if size(HeelStrikeLocs,2)>1
    HeelStrikeLocs = HeelStrikeLocs';
end
if size(ToeOffLocs,2)>1
    ToeOffLocs = ToeOffLocs';
end
while strcmp(cont, 'No, recheck for bad gait events')
    %% find bad HS
    clear badTOIndex badHSIndex badHS badTO
    ToeOffLocs = sort(ToeOffLocs); %sort by event latency, should be alread sorted but just in case you run this section again
    HeelStrikeLocs = sort(HeelStrikeLocs);
    %Difference between gait event latencies; can be used to easily spot outlier events
    TODiff = diff(ToeOffLocs);
    HSDiff = diff(HeelStrikeLocs);
    %Find gait event outliers based on those greater or less than (3) standard
    %deviations
    badTOIndex = find(TODiff >(mean(TODiff)+4*std(TODiff)) | TODiff <(mean(TODiff)-3*std(TODiff)));
    if isempty(badTOIndex)
        badTOIndex = find(TODiff >1024 | TODiff < 256); %check fo super close together or far apart events
    end
    badTO = ToeOffLocs(badTOIndex);
    badHSIndex = find(HSDiff >(mean(HSDiff)+4*std(HSDiff)) | HSDiff <(mean(HSDiff)-3*std(HSDiff)));
    if isempty(badHSIndex)
        badHSIndex = find(HSDiff >1024 | HSDiff < 256);
    end
    badHS = HeelStrikeLocs(badHSIndex);





    %% Fix Heel Strikes
    %Fix unmarked/mislabeled heel strikes
    newHS_EventEstimate = [];
    disp(['Number of bad heel strikes found: ', num2str(length(badHSIndex))]);
    for i = 1:length(badHSIndex)
        clear badHSwindow repeated_events
        previously_added_events = [];
        cursor_info = [];
        userinput= [];
        y = 0;
        if badHSIndex(i) <= 5
            badHSwindow = HeelStrikeLocs(1:badHSIndex(i)+5);%if first bad event is within the first 5 gait events, have window start at 1
            avgHSdiff = mean(diff(HeelStrikeLocs));  %if you don't have 5 previous events to take the difference of, just use the mean of all heel strike locations
        else
            %try
            badHSwindow = HeelStrikeLocs(badHSIndex(i)-5:(min([badHSIndex(i)+5,length(HeelStrikeLocs)])));
            avgHSdiff = mean(diff(HeelStrikeLocs(badHSIndex(i)-5:badHSIndex(i)))); %average latency between previous 5 heel strikes
            %             catch
            %                 badHSwindow = HeelStrikeLocs(1:badHSIndex(i)+5);%if first bad event is within the first 5 gait events, have window start at 1
            %                 avgHSdiff = mean(diff(HeelStrikeLocs));  %if you don't have 5 previous events to take the difference of, just use the mean of all heel strike locations
            %             end
        end
        figure; plot(round(badHSwindow(1)):round(badHSwindow(end)),GRF(round(badHSwindow(1)):round(badHSwindow(end))),'Color',[0.7 0.7 0.7]);
        hold on;
        stem(badHSwindow, 1000*ones(size(badHSwindow)),'k'); set(gcf, 'Position',[0 200 1500 650])
        xlabel('Frame'); ylabel('Ground Reaction Force (N)');

        newHS_EventEstimate = HeelStrikeLocs(badHSIndex(i))+avgHSdiff; %new estimated heel strike is the average time away from the previous HS
        stem(newHS_EventEstimate, 900*ones(size(newHS_EventEstimate,1)), ':m','LineWidth',2); %plots estimated heel strike, dotted green line
        for iv = 1:length(newHSLocs) %mark any new events added to previous windows so you don't accidentally mark a spot more than once
            if newHSLocs(iv) > badHSwindow(1) & newHSLocs(iv) < badHSwindow(end)
                stem(newHSLocs(iv), 1000, ':k')
                previously_added_events = [previously_added_events; newHSLocs(iv)]; %store these, will plot again later
            end
        end
        title(['Bad Heel Strike #', num2str(i),'/',num2str(length(badHSIndex))])

        disp('To remove a gait event, right-click on event and export cursor data to workspace')
        pause;
        userinput = input('\nIs the new Heel Strike event OK? 1=YES, 0= NO, 2= manual, 3=Skip: ');
        if isempty(userinput)
            warning('Please enter a valid input argument')
            userinput = input('\nIs the new toe-off event OK? 1=YES, 0= NO, 2= manual, 3=Skip: ');
        elseif userinput ~= [0 1 2 3]
            warning('Please enter a valid input argument')
            userinput = input('\nIs the new toe-off event OK? 1=YES, 0= NO, 2= manual, 3=Skip: ');
        end

        hold off;
        if size(newHS_EventEstimate,1)>1 & userinput ==1
            which_event = input('Which event estimate(s) do you want? (e.g. 1,2,[1 2]): ');
        else
            which_event =1;
        end

        try
            cursor_info = evalin('base','cursor_info');
        catch
            cursor_info = [];
        end

        if ~isempty(cursor_info)
            %disp(cursor_info)
            remove_HSevent= [remove_HSevent; cursor_info.Position(1,1)];%To remove gait event, right-click on event and export cursor data to workspace. This takes the x-coordinate from the cursor position and saves it in a new array
        end
        if userinput == 1
            newHSLocs = [newHSLocs; newHS_EventEstimate(which_event,:)];
            HeelStrikeLocs = [HeelStrikeLocs; newHS_EventEstimate(which_event,:)];  %Add in fixed gait events
            xx = length(which_event)-1;
            close;
        elseif userinput == 0
            clear events;
            n_events = input('\nEnter number of events to mark: ');
            events = zeros(n_events,2);
            %             ButtonHandle = uicontrol('Style', 'PushButton', ...
            %                 'String', 'Skip', 'Position', [5, 20, 150, 30], ...
            %                 'Callback', @skipevent);

            for i = 1:n_events
                fprintf('Mark heel strike event %d (press enter to skip)\n', i);
                [x, y] = ginput(1);
                if isempty(x)
                    fprintf('Skipped event %d\n', i);
                    continue;  % Skip the rest of the loop and move to the next event
                end
                events(i, :) = [x, y];
                hold on;
                stem(x, 1000, 'm','LineWidth',2);
            end
            hold off;
            pause(1)

            if isempty(n_events) || n_events ==0 %Skip button was pressed, function emptied multiselect variable
                disp('Event has been skipped')
                xx = [];
                close;
            else %"Select multiple events" button was pressed, locations stored in multiselect array
                newHSLocs = [newHSLocs; events(:,1)];
                xx = size(events,1)-1; % xx keeps track of how many events were selected beyond the end of newHSLocs
                HeelStrikeLocs = [HeelStrikeLocs; newHSLocs(end-xx:end)];
                close;
            end
            disp('Finished selecting')
        elseif userinput == 3
            disp('SKIP')
            xx = []; %no event selected
            close;
        elseif userinput ==2 %manually enter
            g = input('Enter new HS locations:');
            if size(g,2) ~=1
                g = g';
            end

            if ~isempty(g)
                newHSLocs = [newHSLocs; g];

                HeelStrikeLocs = [HeelStrikeLocs; newHSLocs(end)]; %Add in fixed gait events
                fprintf('\%i events added',length(newHSLocs))
                xx = length(newHSLocs)-1;
            else
                xx = [];
            end
        end
        close;

        if ~isempty(xx)
            if xx ~=0
                figure; plot(round(badHSwindow(1)):round(badHSwindow(end)),GRF(round(badHSwindow(1)):round(badHSwindow(end))),'Color',[0.7 0.7 0.7]);
                hold on;
                stem(badHSwindow, 1000*ones(size(badHSwindow)),'k');
                stem(newHSLocs(end-xx:end), 900*ones(size(newHSLocs(end-xx:end))),'m','LineWidth',2);
                set(gcf, 'Position',[0 200 1500 650])
                title(['New gait event #', num2str(i)])


                if ~isempty(previously_added_events)
                    stem(previously_added_events, 1000*ones(size(previously_added_events)), ':k')
                    legend('Fz', 'old gait events', 'previously added events','new gait events')
                end
                hold off;
                set(gcf, 'Position',[0 200 1500 650])
            elseif (xx == 0 && ~isempty(newHSLocs))
                figure; plot(round(badHSwindow(1)):round(badHSwindow(end)),GRF(round(badHSwindow(1)):round(badHSwindow(end))));
                hold on;
                stem(badHSwindow, 1000*ones(size(badHSwindow)),'k');
                stem(newHSLocs(end), 900 ,'m','LineWidth',2);
                title(['New gait event #', num2str(i)])
                legend('Fz', 'old gait events','new gait event')
                if ~isempty(previously_added_events)
                    stem(previously_added_events, 1000*ones(size(previously_added_events)), ':k')
                    legend('Fz', 'old gait events', 'previously added events','new gait events')
                end
                hold off
                set(gcf, 'Position',[0 200 1500 650])
            end
        end

        repeated_events = diff(HeelStrikeLocs) == 0; %delete any events that were accidentally repeated
        HeelStrikeLocs(repeated_events) = [];
        disp('Press space bar to continue')
        pause;
        close;
        %         fprintf('\nsize of remove_HSevent %i\n', size(remove_HSevent,2))
        %         disp(remove_HSevent)
    end

    % Remove events manually
    if ~isempty(remove_HSevent)
        for v=1:length(remove_HSevent)
            %find where the heel strike locations equal the locations that were manually selected for removal; create an index
            newHSLocs( newHSLocs== remove_HSevent(v))
            HeelStrikeLocs(HeelStrikeLocs == remove_HSevent(v)) = []; %remove events at this index (cursor locations user selected for manual event removal)
        end
    end
    %% Fix Toe Offs
    %Fix any unmarked/mislabeled toe offs
    disp('Toe Offs')
    newTO_EventEstimate = [];
    disp(['Number of bad toe offs found: ', num2str(length(badTOIndex))]);
    for i = 1:length(badTOIndex) %spans through each gait event outlier
        clear badTOwindow repeated_events
        cursor_info=[];
        previously_added_events=[];
        y = 0;
        if badTOIndex(i) <= 5
            badTOwindow = ToeOffLocs(1:badTOIndex(i)+5);%if first bad event is within the first 5 gait events, have window start at 1
            avgTOdiff = mean(diff(ToeOffLocs)); %if you don't have 5 previous events to take the difference of, just use the mean of all toe off locations
        else
            %try
            badTOwindow = ToeOffLocs(badTOIndex(i)-5:(min([badTOIndex(i)+5,length(ToeOffLocs)])));
            avgTOdiff = mean(diff(ToeOffLocs(badTOIndex(i)-5:badTOIndex(i)))); %average latency between previous 5 toe offs
            %             catch
            %                 badTOwindow = ToeOffLocs(1:badTOIndex(i)+5);%if first bad event is within the first 5 gait events, have window start at 1
            %                 avgTOdiff = mean(diff(ToeOffLocs)); %if you don't have 5 previous events to take the difference of, just use the mean of all toe off locations
            %             end
        end
        %         time_window_ind = find(EEG_time> badTOwindow(1) & EEG_time<badTOwindow(end));
        %         figure; plot(EEG_time(time_window_ind), GRF(time_window_ind,:));
        figure; plot(round(badTOwindow(1)):round(badTOwindow(end)),GRF(round(badTOwindow(1)):round(badTOwindow(end))),'Color',[0.7 0.7 0.7]);
        hold on;
        stem(badTOwindow, 1000*ones(size(badTOwindow)),'k');set(gcf, 'Position',[0 200 1500 650])

        %estimate toe-off based on average of events
        newTO_EventEstimate = ToeOffLocs(badTOIndex(i))+avgTOdiff; %new estimated toe off is the average time away from the previous TO

        stem(newTO_EventEstimate, 900*ones(size(newTO_EventEstimate,1)), ':g','LineWidth',2);%plots estimated toe off, dotted green line
        legend('Fz', 'old gait events','estimated gait event')
        title(['Bad Toe Off #', num2str(i),'/',num2str(length(badTOIndex))]);
        for iv = 1:length(newTOLocs) %mark any new events added to previous windows so you don't accidentally mark a spot more than once
            if newTOLocs(iv) > badTOwindow(1) & newTOLocs(iv) < badTOwindow(end)
                stem(newTOLocs(iv), 1000, ':k')
                previously_added_events = [previously_added_events; newTOLocs(iv)]; %store these, will plot again later
                legend('Fz', 'old gait events', 'previously added events','new gait events')
            end
        end
        disp('To remove a gait event, right-click on event and export cursor data to workspace')
        pause;
        userinput = input('Is the new Toe-off event OK? 1=YES, 0= NO, 2= manual entry, 3=Skip: ');
        if isempty(userinput)
            warning('Please enter a valid input argument')
            userinput = input('\nIs the new toe-off event OK? 1=YES, 0= NO, 2= manual, 3=Skip: ');
        elseif userinput ~= [0 1 2 3]
            warning('Please enter a valid input argument')
            userinput = input('\nIs the new toe-off event OK? 1=YES, 0= NO, 2= manual, 3=Skip: ');
        end

        hold off;
        if size(newTO_EventEstimate,1)>1 & userinput ==1;
            which_event = input('Which event estimate(s) do you want? (e.g. 1,2,[1 2]): ');
        else
            which_event =1;
        end

        try
            cursor_info = evalin('base','cursor_info');
        catch
            cursor_info = [];
        end
        if ~isempty(cursor_info)%Manually remove events in window
            %disp(cursor_info)
            remove_TOevent= [remove_TOevent; cursor_info.Position(1,1)];%To remove gait event, right-click on event and export cursor data to workspace. This takes the x-coordinate from the cursor position and saves it in a new array
        end
        if userinput == 1
            newTOLocs = [newTOLocs; newTO_EventEstimate(which_event,:)];
            ToeOffLocs = [ToeOffLocs;  newTO_EventEstimate(which_event,:)];  %Add in fixed gait events
            xx = length(which_event)-1;
            close;
        elseif userinput == 0
            clear events;
            n_events = input('\nEnter number of events to mark: ');
            events = zeros(n_events,2);
            %             ButtonHandle = uicontrol('Style', 'PushButton', ...
            %                 'String', 'Skip', 'Position', [5, 20, 150, 30], ...
            %                 'Callback', @skipevent);

            for i = 1:n_events
                fprintf('Mark toe-off event %d\n', i);
                [x, y] = ginput(1);
                if isempty(x)
                    fprintf('Skipped event %d (press enter to skip)\n', i);
                    continue;  % Skip the rest of the loop and move to the next event
                end
                events(i, :) = [x, y];
                hold on;
                stem(x, 1000, 'g','LineWidth',2);
            end
            hold off;
            pause(1)

            if isempty(n_events) || n_events ==0 %Skip button was pressed, function emptied multiselect variable
                disp('Event has been skipped')
                xx = [];
                close;
            else %"Select multiple events" button was pressed, locations stored in multiselect array
                newTOLocs = [newTOLocs; events(:,1)];
                xx = size(events,1)-1; % xx keeps track of how many events were selected beyond the end of newHSLocs

                close;
            end

        elseif userinput == 3
            disp('Event has been skipped')
            xx=[];
            close;
        elseif userinput ==2 %manually enter
            g = input('Enter new toe-off locations:');
            if size(g,2) ~=1
                g = g';
            end
            if ~isempty(g)
                newTOLocs = [newTOLocs; g];
                ToeOffLocs = [ToeOffLocs; newTOLocs(end)]; %Add in fixed gait events
                fprintf('\%i events added',length(newTOLocs))
                xx = length(newTOLocs)-1;
            else
                xx = [];
            end

        else
            disp('waiting for input')
            pause;
        end
        close;

        if ~isempty(xx)
            if xx ~=0

                figure; plot(round(badTOwindow(1)):round(badTOwindow(end)),GRF(round(badTOwindow(1)):round(badTOwindow(end))),'Color',[0.7 0.7 0.7]);
                hold on;
                stem(badTOwindow, 1000*ones(size(badTOwindow)),'k');
                stem(newTOLocs(end-xx:end), 900*ones(size(newTOLocs(end-xx:end))),'g','LineWidth',2);
                title(['New gait event #',num2str(i)]); legend('Fz', 'old gait events', 'new gait events')
                if ~isempty(previously_added_events)
                    stem(previously_added_events, 1000*ones(size(previously_added_events)), ':k')
                    legend('Fz', 'old gait events', 'previously added events','new gait events')
                end
                hold off;
                set(gcf, 'Position',[0 200 1500 650])
            elseif (xx==0 & ~isempty(newTOLocs))
                figure; plot(round(badTOwindow(1)):round(badTOwindow(end)),GRF(round(badTOwindow(1)):round(badTOwindow(end))),'Color',[0.7 0.7 0.7]);
                hold on;
                stem(badTOwindow, 1000*ones(size(badTOwindow)),'k');
                stem(newTOLocs(end), 900 ,'g','LineWidth',2);
                set(gcf, 'Position',[0 200 1900 650])
                title(['New gait event #', num2str(i)])
                legend('Fz', 'old gait events', 'new gait event')
                if ~isempty(previously_added_events)
                    stem(previously_added_events, 1000*ones(size(previously_added_events)), ':k')
                    legend('Fz', 'old gait events', 'new gait events','previously added events')
                end
                hold off
                set(gcf, 'Position',[0 200 1500 650])
            end
        end
        repeated_events = find(diff(ToeOffLocs) == 0); %delete any events that were accidentally repeated
        ToeOffLocs(repeated_events) = [];
        disp('Press space bar to continue')
        pause;
        close;
    end
    % Remove events manually
    %Remove events that were selected manually using the cursor
    if ~isempty(remove_TOevent)
        for v = 1:length(remove_TOevent)
            %find where the Toe off locations equal the locations that were manually selected for removal
            newTOLocs(newTOLocs== remove_TOevent(v)) =[];
            ToeOffLocs(ToeOffLocs == remove_TOevent(v)) = []; %remove event at this cursor location index
        end
    end
    if ~isempty(badTOIndex) || ~isempty(badHSIndex)
        cont = questdlg('Do you want to continue to the next step?','Save Gait Events','Yes','No, recheck for bad gait events','Cancel','No, recheck for bad gait events');
    else
        cont = 'Yes'; %continue
    end

end
remove_HSevent = unique(remove_HSevent);
remove_TOevent = unique(remove_TOevent);

% while strcmp(cont, 'No, recheck for bad gait events')
% figure; plot(GRF); hold on;
% stem(HeelStrikeLocs, mean(GRF)*ones(size(ToeOffLocs)), 'c');
% stem(ToeOffLocs, mean(GRF)*ones(size(ToeOffLocs)), 'm');
% legend('forceplate','heel strike','toe-off')
% xlabel('Time (s)'); ylabel('Ground Reaction Force (N)')
% title('Gait Events')
%
% %manually remove events
% disp('To remove a gait event, right-click on event and export cursor data to workspace')
% cursor_info = evalin('base','cursor_info');
% if ~isempty(cursor_info)%Manually remove events in window
%     %disp(cursor_info)
%     remove_TOevent= [remove_TOevent; cursor_info.Position(1,1)];%To remove gait event, right-click on event and export cursor data to workspace. This takes the x-coordinate from the cursor position and saves it in a new array
% end
%
% userinput = input('Do you want to manually remove or add gait events? Press 1=YES, 0=NO: ');
%
% if userinput == 1
%     %manually adjust heel strikes
%     clear g;
%     multiselect = 0; %reset array
%     save('multiselect.mat'); %save array
%     disp('Select Heel Strike event')
%     ButtonHandle = uicontrol('Style', 'PushButton', ...
%         'String', 'Skip', 'Position', [5, 20, 150, 30], ...
%         'Callback', @skipevent);
%     ButtonHandle2 = uicontrol('Style', 'PushButton', ...
%         'String', 'Select multiple events', 'Position', [5, 50, 150, 30], 'Callback', @multiEventSelection);%Creates a "Stop loop" button in the figure that the user can push to escape for loop if they don't want to select a gait event in that window
%     g = ginput(1); %user marks new gait event in figure
%     disp('Press space bar to continue')
%     pause;
%     load('multiselect.mat') %have to load files saved in button functions because variables won't save to workspace
%     multiselect = evalin('base', 'multiselect');
%     if isempty(multiselect) %Skip button was pressed, function emptied multiselect variable
%         disp('Event has been skipped')
%         xx = [];
%         close;
%     elseif (~isempty(multiselect) & multiselect ~= 0) %"Select multiple events" button was pressed, locations stored in multiselect array
%         g = multiselect;
%         disp('Multiple events selected')
%         newHSLocs = [newHSLocs; g(:,1)];
%         %                 HScont = HeelStrikeLocs(badHSIndex(i)+1) - newHSLocs(end)>(2*avgHSdiff);  %Check if there are multiple unmarked gait events in this region (i.e., if the next event happens later than 2x the avg time difference)
%         xx = length(multiselect)-1; % xx keeps track of how many events were selected beyond the end of newHSLocs
%         HeelStrikeLocs = [HeelStrikeLocs; newHSLocs(end-xx:end)];
%         close;
%     elseif (multiselect== 0) %no buttons pressed, only one gait event selected
%         newHSLocs = [newHSLocs; g(end,1)];
%         %                 HScont = HeelStrikeLocs(badHSIndex(i)+1) - newHSLocs(end)>(2*avgHSdiff);  %Check if there are multiple unmarked gait events in this region (i.e., if the next event happens later than 2x the avg time difference)
%         HeelStrikeLocs = [HeelStrikeLocs; newHSLocs(end)]; %Add in fixed gait events
%         disp('One event selected');
%         xx = 0; %only one event selected
%         close;
%     end
%     disp('Finished selecting heel-strikes')
%
%     %manually adjust toe-offs
%     clear g;
%     multiselect = 0; %reset array
%     save('multiselect.mat'); %save array
%     disp('Select Toe Off event')
%     ButtonHandle = uicontrol('Style', 'PushButton', ...
%         'String', 'Skip', 'Position', [5, 20, 150, 30], ...
%         'Callback', @skipevent);
%     ButtonHandle2 = uicontrol('Style', 'PushButton', ...
%         'String', 'Select multiple events', 'Position', [5, 50, 150, 30], 'Callback', @multiEventSelection);%Creates a "Stop loop" button in the figure that the user can push to escape for loop if they don't want to select a gait event in that window
%     g = ginput(1); %user marks new gait event in figure
%     disp('Press space bar to continue')
%     pause;
%     load('multiselect.mat') %have to load files saved in button functions because variables won't save to workspace
%     multiselect = evalin('base', 'multiselect');
%     if isempty(multiselect) %Skip button was pressed, function emptied multiselect variable
%         disp('Event has been skipped')
%         xx = [];
%         close;
%     elseif ~isempty(multiselect) & multiselect ~= 0 %"Select multiple events" button was pressed, locations stored in multiselect array
%         g = multiselect;
%         disp('Multiple events selected')
%         newTOLocs = [newTOLocs; g(:,1)];
%         %                 TOcont = ToeOffLocs(badTOIndex(i)+1) - newTOLocs(end)>(2*avgTOdiff);  %Check if there are multiple unmarked gait events in this region (i.e., if the next event happens later than 2x the avg time difference)
%         xx = length(multiselect)-1;
%         ToeOffLocs = [ToeOffLocs; newTOLocs(end-xx:end)];%Add in fixed gait events
%         close;
%     elseif (multiselect== 0) %no buttons pressed, only one gait event selected
%         newTOLocs = [newTOLocs; g(end,1)];
%         %Add in fixed gait events
%         ToeOffLocs = [ToeOffLocs; newTOLocs(end)];
%         disp('One event selected');
%         xx = 0;
%         close;
%     end
%     disp('Finished selecting toe-offs')
%
% end
%     % Remove events manually
%     %Remove events that were selected manually using the cursor
%         % Remove events manually
%     if ~isempty(remove_HSevent)
%         for v=1:length(remove_HSevent)
%             remove_HSevent_index = find(HeelStrikeLocs == remove_HSevent(v)); %find where the heel strike locations equal the locations that were manually selected for removal; create an index
%             HeelStrikeLocs(remove_HSevent_index) = []; %remove events at this index (cursor locations user selected for manual event removal)
%         end
%     end
%     if ~isempty(remove_TOevent)
%         for v = 1:length(remove_TOevent)
%             remove_TOevent_index = find(ToeOffLocs == remove_TOevent(v)); %find where the Toe off locations equal the locations that were manually selected for removal
%             ToeOffLocs(remove_TOevent_index) = []; %remove event at this cursor location index
%         end
%     end
%
% cont = questdlg('Do you want to continue and save gait events?','Save Gait Events','Yes','No, recheck for bad gait events','Cancel');
% end
%
% remove_HSevent = unique(remove_HSevent);
% remove_TOevent = unique(remove_TOevent);

end



function [HS, TO] = check_gaitevents(HS,TO,GRF)

figure;
plot( GRF,'Color',[0.7 0.7 0.7]); hold on;
stem(HS, max(GRF)*ones(size(HS)), 'c');
stem(TO,max(GRF)*ones(size(TO)), 'm')
legend('forceplate','HS','TO')
ylabel('Ground Reaction Force (N)')
xlabel('Frame')
xlim([HS(1)-100 HS(end)+100])
ylim([-100 1200])
set(gcf, 'Position',[0 200 1500 650])
pause;
userinput = input('\nGood to go? (1|0): ');
while isempty(userinput)
    warning('Please enter a valid input')
    userinput = input('\nGood to go? (1|0): ');
end

if userinput ~= 0 && userinput ~=1
    warning('Please enter a valid input')
    userinput = input('\nGood to go? (1|0): ');
end

while userinput ==0
    cursor_info = [];
    remove_event = [];
    newHSLocs = []; newTOLocs =[];
    cursor_info = evalin('base','cursor_info');
    if ~isempty(cursor_info)%Manually remove events in window
        %disp(cursor_info)
        remove_event= [remove_event; cursor_info.Position(1,1)];%To remove gait event, right-click on event and export cursor data to workspace. This takes the x-coordinate from the cursor position and saves it in a new array
    end

    clear events;
    n_events = input('\nEnter number of HS events to mark: ');
    events = zeros(n_events,2);
    %             ButtonHandle = uicontrol('Style', 'PushButton', ...
    %                 'String', 'Skip', 'Position', [5, 20, 150, 30], ...
    %                 'Callback', @skipevent);
    if ~isempty(n_events)
        for i = 1:n_events
            fprintf('Mark HS event %d\n', i);
            [x, y] = ginput(1);
            if isempty(x)
                fprintf('Skipped event %d (press enter to skip)\n', i);
                continue;  % Skip the rest of the loop and move to the next event
            end
            events(i, :) = [x, y];
            hold on;
            stem(x, max(GRF), ':c','LineWidth',2);
        end
    end
    hold off;
    pause(1)

    if isempty(n_events) || n_events ==0 %Skip button was pressed, function emptied multiselect variable
        disp('Event has been skipped')
        xx = [];
    else %"Select multiple events" button was pressed, locations stored in multiselect array
        newHSLocs = [newHSLocs; events(:,1)];
        xx = size(events,1)-1; % xx keeps track of how many events were selected beyond the end of newHSLocs
        HS = [HS; newHSLocs(end-xx:end)];

    end

    n_events = input('\nEnter number of TO events to mark: ');
    events = zeros(n_events,2);
    %             ButtonHandle = uicontrol('Style', 'PushButton', ...
    %                 'String', 'Skip', 'Position', [5, 20, 150, 30], ...
    %                 'Callback', @skipevent);
    if ~isempty(n_events)
        for i = 1:n_events
            fprintf('Mark toe-off event %d\n', i);
            [x, y] = ginput(1);
            if isempty(x)
                fprintf('Skipped event %d (press enter to skip)\n', i);
                continue;  % Skip the rest of the loop and move to the next event
            end
            events(i, :) = [x, y];
            hold on;
            stem(x, max(GRF), ':m','LineWidth',2);
        end
    end
    hold off;
    pause(1)

    if isempty(n_events) || n_events ==0 %Skip button was pressed, function emptied multiselect variable
        disp('Event has been skipped')
        xx = [];
    else %"Select multiple events" button was pressed, locations stored in multiselect array
        newTOLocs = [newTOLocs; events(:,1)];
        xx = size(events,1)-1; % xx keeps track of how many events were selected beyond the end of newHSLocs
        TO = [TO; newTOLocs(end-xx:end)];
    end
    remove_event = unique(remove_event);
    if ~isempty(remove_event)
        for eventi =1:length(remove_event)
            idx = find(TO== remove_event(eventi));
            if ~isempty(idx)
                TO(idx)=[];
            else
                idx = find(HS== remove_event(eventi));
                HS(idx)=[];
            end
        end
    end

    close;
    figure;
    plot( GRF,'Color',[0.7 0.7 0.7]); hold on;
    stem(HS, max(GRF)*ones(size(HS)), 'c');
    stem(TO,max(GRF)*ones(size(TO)), 'm')
    legend('forceplate','HS','TO')
    ylabel('Ground Reaction Force (N)')
    xlabel('Frame');
    xlim([HS(1)-100 HS(end)+100]);
    ylim([-100 1200])
    set(gcf, 'Position',[0 200 1500 650])
    pause;

    userinput = input('\nGood to go? (1|0): ');
    if userinput ~= 0 && userinput ~=1
        warning('Please enter a valid input')
        userinput = input('\nGood to go? (1|0): ');
    end
end
TO = sort(TO);
HS = sort(HS);
close;
end


