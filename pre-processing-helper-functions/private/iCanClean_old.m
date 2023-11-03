% iCanClean() - use canonical correlation analysis to remove mutual
%               bad components between two sets of channels
% Usage:
%   >>  out = iCanClean( in1, in2 );
%
% Inputs:
%   in1     - first input of the function
%   in1     - first input of the function
%
% Outputs:
%   out     - output of the function
%
% See also:
%   POP_iCanClean, EEGLAB

% Copyright (C) <2020>  <Ryan Downey>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function cleanEEG = iCanClean( EEG, xChan, yChan, visualizeFinalResults, params)

if nargin < 1
    help iCanClean;
    return;
end

%% initialize EEG structure that will eventually be our clean data and keep track of time
cleanEEG = EEG;
tic;

%% check/set parameters
params = iCC_setParams(params)
fprintf('HERE ARE YOUR PARAMS')
cleanEEG.iCanClean.params = params;
%% define time windows (allows use of moving window)
t = EEG.times/1000; %seconds
tTotal = EEG.times(end)/1000;
tStart = 0:params.stepSize:tTotal-params.windowLength; %array of times defining start of window (moving or static)
if isempty(tStart) tStart = 0; end %helping to avoid a bug momentarily
tEnd = tStart+params.windowLength; %array of end times

%start starts figure (fill in with data later)
if params.plotStatsOn
    myStatsFig = figure();
    myStatsFig_R = subplot(2,1,1); xlabel('Component #'); ylabel('R or Rsq need to check');
    myStatsFig_p = subplot(2,1,2); xlabel('Component #'); ylabel('p val');%set(myStatsFig_p,'yscal','log')
end

%% optionally use "external" data within each window by first computing CCA on the entire dataset
if params.calcCCAonWholeData == true
    X = EEG.data(xChan,:)';
    Y = EEG.data(yChan,:)';

    %calculate noise sources using CCA
    if any(isnan(X(:))) | any(isnan(Y(:)))
        disp('I could not do calib with entire dataset b/c you have NaN')
        A_calib = [];
        B_calib = [];
    else %it's safe to calc

        [A_calib,B_calib,R_calib,U_calib,V_calib,STATS_calib] = canoncorr(X, Y );
        X_MC = (X - mean(X,1));
        U = X_MC*A_calib;
        A_inv_calib = pinv(A_calib(:,:));
        fakeWinv_calib = mrdivide(X_MC(:,:)',U(:,:)');

        clear X Y X_MC Y_MC U V
    end
else
    A_calib = [];
    B_calib = [];
end
%% go thru each window, grab the data, calculate potential noise sources, and clean x and or y with method of choice
% msg_n = 0; %used for trying to display output in one line only
numWindows = length(tStart);
numNoiseCompsRemoved = zeros(1,numWindows);
for window_i = 1:numWindows

    window = find( t>=tStart(window_i) & t<=tEnd(window_i) );
    broadWindow = find( t>=tStart(window_i)-params.extraTime_pre & t<=tEnd(window_i)+params.extraTime_post );
    [lia1, locb1]  = ismember(window,broadWindow);

    %grab data
    X = EEG.data(xChan,broadWindow)';

    if isempty(yChan)
        disp('I guess you want to do auto lagged CCA??');
        lagAmount = 1; %samples
        disp(['Using a lag of ',num2str(lagAmount),' samples.']);
        shiftedBroadWindow = broadWindow + lagAmount;
        Y = EEG.data(xChan,shiftedBroadWindow)';
    else
        Y = EEG.data(yChan,broadWindow)';
    end

    %calculate noise sources using CCA
    if any(isnan(X(:))) | any(isnan(Y(:))) %if missing data (NaN)
        disp('I could not clean this section because you have NaN')
        nT = size(X,1); nX = size(X,2); nY = size(Y,2); nS = min(nX,nY);
        R = zeros(nS,1);
        U = zeros(nT,nS);
        V = U;
    else %we have data available so it's safe to run
        if isempty(A_calib) || isempty(B_calib) %if we don't use to use external data and instead want to calc cca on this smaller window
            [A,B,R,U,V,STATS] = iCC_calcSources(X,Y,params);
        else %we were given external A and B matrices to use
            A = A_calib;
            B = B_calib;
            X_MC = (X - mean(X,1));
            Y_MC = (Y - mean(Y,1));
            U = X_MC*A;
            V = Y_MC*B;
            [tempRHO,tempPVAL] = corr(U,V);
            R = diag(tempRHO);
            STATS.p = diag(tempPVAL);
        end
    end

    % define bad noise sources at source level (cancor stats)
    badComps = find(R.^2>params.rhoSqThres_source & STATS.p < params.pThres_source); % find(R>0.98 & STATS.p < 1E-4) sometimes 1E-3
    numNoiseCompsRemoved(window_i) = length(badComps);

    if params.plotStatsOn
        %plot stats
        hold(myStatsFig_R,'off');stem(myStatsFig_R, R); hold(myStatsFig_R,'on'); stem(myStatsFig_R, badComps,R(badComps),'r'); set(myStatsFig_R,'ylim',[-1 1]);
        hold(myStatsFig_p,'off');stem(myStatsFig_p, R.^2); hold(myStatsFig_p,'on'); stem(myStatsFig_p, badComps,R(badComps).^2,'r'); set(myStatsFig_p,'ylim',[0 1]);
        %     stem(myStatsFig_p, STATS.p);
        drawnow;
    end
    %% define which noise sources to use to clean each of the original channels
    %i.e. do you want to clean using mixtures of X or mixtures of Y?
    if strcmpi(params.cleanXwith,'X') | strcmpi(params.cleanXwith,'U')
        xNoiseSources = U(:,badComps);
    elseif strcmpi(params.cleanXwith,'Y') | strcmpi(params.cleanXwith,'V')
        xNoiseSources = V(:,badComps);
    elseif strcmpi(params.cleanXwith,'XY') | strcmpi(params.cleanXwith,'UV')
        xNoiseSources = (U(:,badComps)+V(:,badComps))./2;
    else
        erorr('you did not define params.cleanXwith properly');
    end

    if strcmpi(params.cleanYwith,'X') | strcmpi(params.cleanYwith,'U')
        yNoiseSources = U(:,badComps);
    elseif strcmpi(params.cleanYwith,'Y') | strcmpi(params.cleanYwith,'V')
        yNoiseSources = V(:,badComps);
    elseif strcmpi(params.cleanXwith,'XY') | strcmpi(params.cleanXwith,'UV')
        yNoiseSources = (U(:,badComps)+V(:,badComps))./2;
    else
        erorr('you did not define params.cleanYwith properly');
    end

    %% clean channels
    %clean x
    %     tempClean = RyanCCA_cleanChansWithNoiseSources(X,xNoiseSources,params.noiseRemovalMethod,params)';
    if isempty(A_calib) || isempty(B_calib)%ryan new 11/30/2020 quick test to see if we want to clean another way
        tempClean = iCC_cleanChansWithNoiseSources(X,xNoiseSources,params.noiseRemovalMethod,params)';
    else %we want to use calib data to clean
        disp('THIS IS EXPERIMENTAL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        X_MC = (X - mean(X,1));
        %         Y_MC = (Y - mean(Y,1));
        U = X_MC*A_calib;

        %         %option 1: reject with psuedo inverse
        %         X_est = X_MC*A_calib(:,badComps)*A_inv_calib(badComps,:);

        %option 2: reject with matrix division (least square?)
        %            fakeWinv_calib = mrdivide(X_MC(:,:)',U(:,:)'); disp('doing local scaling'); %calc locally if line executed. comment out if you want to use what was calc from larger dataset
        X_est = (fakeWinv_calib(:,badComps)*(X_MC*A_calib(:,badComps))')';
        %         X_est = X_MC*A_calib(:,badComps)*fakeWinv_calib(:,badComps)'; %equiv expression to line above? check and see

        tempClean = (X-X_est)';
    end
    cleanEEG.data(xChan,window) = tempClean(:,locb1); clear tempClean;

    %clean y if desired
    if params.cleanYBool && ~isempty(yChan)
        tempClean = iCC_cleanChansWithNoiseSources(Y,yNoiseSources,params.noiseRemovalMethod,params)';
        cleanEEG.data(yChan,window) = tempClean(:,locb1); clear tempClean;
    end

    %update user on progress (assuming moving window)
    %     fprintf(repmat('\b',1,msg_n+1));
    if params.giveCleaningUpdates
        msg = ['iCanClean cleaned [',num2str([tStart(window_i) tEnd(window_i)]),']s. Removed ',num2str(length(badComps)),' bad sources in this window'];
        disp(msg);
    end
end

%% finish keeping track of time
fprintf('\niCC cycle complete\n')
toc
%% update EEG structure on the way out
cleanEEG.iCanClean.numNoiseCompsRemovedPerWindow = numNoiseCompsRemoved;
cleanEEG.iCanClean.numNoiseCompsRemovedOnAvg = mean(numNoiseCompsRemoved);
disp(['Removed ',num2str(cleanEEG.iCanClean.numNoiseCompsRemovedOnAvg),' noise components on average (across all time windows cleaned)']);
%% visualize final results
if visualizeFinalResults
    disp('Visualizing results. This may take a while.');
    %     RyanCCA_visualizeResults(cleanEEG,EEG,xChan(1:4:end),yChan(1:8:end));
    iCC_visualizeResults(cleanEEG,EEG,xChan(1:4:end),[]);
end
end
