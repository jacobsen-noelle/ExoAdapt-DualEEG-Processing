function plot_gaitevents(EEG,xax)

% twindow  = [30000:60000];
% GE= table2struct(gaitdata.GaitEvents);
% Time = gaitdata.Biomech_Data_Tseries.Time(twindow);
% GRF_L = gaitdata.Biomech_Data_Tseries.Left_GRF(twindow);
% GRF_R = gaitdata.Biomech_Data_Tseries.Right_GRF(twindow);
% RHS = [GE(strcmpi([GE.EventType],'RHS')).Time];
% LHS = [GE(strcmpi([GE.EventType],'LHS')).Time];
% LTO = [GE(strcmpi([GE.EventType],'LTO')).Time];
% RTO = [GE(strcmpi([GE.EventType],'RTO')).Time];
% 
% figure; tiledlayout(2,1);
% ax1= nexttile;
% plot(Time, GRF_L); hold on;
% stem(LHS, mean(GRF_L)*ones(size(LHS)), 'c');
% stem(LTO,mean(GRF_L)*ones(size(LTO)), 'm')
% legend('forceplate','LHS','LTO')
% xlabel('Time (s)'); ylabel('Ground Reaction Force (N)')
% title('Left')
% 
% ax2 = nexttile;
% plot(Time, GRF_R); hold on;
% stem(RHS, mean(GRF_R)*ones(size(RHS)), 'c');
% stem(RTO,mean(GRF_R)*ones(size(RTO)), 'm');
% linkaxes([ax1, ax2],'x')
% legend('forceplate','RHS','RTO')
% xlabel('Time (s)'); ylabel('Ground Reaction Force (N)')
% title('Right')
GRF_L = EEG.data(strcmp({EEG.chanlocs.labels},'Left_GRF'),:);
GRF_R = EEG.data(strcmp({EEG.chanlocs.labels},'Right_GRF'),:);


RHS = [EEG.event(strcmpi({EEG.event.type},'RHS')).latency];
LHS = [EEG.event(strcmpi({EEG.event.type},'LHS')).latency];
LTO = [EEG.event(strcmpi({EEG.event.type},'LTO')).latency];
RTO = [EEG.event(strcmpi({EEG.event.type},'RTO')).latency];
if strcmp(xax,'time')
    RHS = RHS/EEG.srate;
    LHS = LHS/EEG.srate;
    LTO = LTO/EEG.srate;
    RTO = RTO/EEG.srate;
    X = EEG.times/1000;
elseif strcmp(xax,'lat')
    X = 1:EEG.pnts;
end

figure; tiledlayout(2,1);
ax1= nexttile;
plot(X, GRF_L,'Color',[0.7 0.7 0.7]); hold on;
stem(LHS, max(GRF_L)*ones(size(LHS)), 'c');
stem(LTO,max(GRF_L)*ones(size(LTO)), 'm')
legend('forceplate','LHS','LTO')
ylabel('Ground Reaction Force (N)')
title('Left')
if strcmp(xax,'time')
    xlabel('Time (s)');
else
    xlabel('Latency')
end

ax2 = nexttile;
plot(X, GRF_R,'Color',[0.7 0.7 0.7]); hold on;
stem(RHS, max(GRF_R)*ones(size(RHS)), 'c');
stem(RTO,max(GRF_R)*ones(size(RTO)), 'm');
linkaxes([ax1, ax2],'x')
legend('forceplate','RHS','RTO')
ylabel('Ground Reaction Force (N)')
title('Right')
if strcmp(xax,'time')
    xlabel('Time (s)');
    fprintf('Event Type \t Average Time Diff (s)\n')
else
    xlabel('Latency')
    fprintf('Event Type \t Average Latenvcy Diff (s)\n')
end

fprintf('\tLHS\t\t   %.2f\n', mean(diff(LHS)))
fprintf('\tLTO\t\t   %.2f\n', mean(diff(LTO)))
fprintf('\tRHS\t\t   %.2f\n', mean(diff(RHS)))
fprintf('\tRTO\t\t   %.2f\n', mean(diff(RTO)))
end