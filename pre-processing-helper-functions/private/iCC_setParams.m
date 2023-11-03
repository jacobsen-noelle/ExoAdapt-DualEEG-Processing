% Ryan Downey
function [params] = iCC_setParams(params)

defaultParams.windowLength = 2; %seconds
if isfield(params,'windowLength')
    if ~isempty(params.windowLength) %if user specified a window length, make that the default step size. 
                                    %note if user also specified a step size, then we will use that instead
        defaultParams.stepSize = params.windowLength; %seconds
    else
        defaultParams.stepSize = defaultParams.windowLength; %seconds
    end
else %if user didn't specify a window length then we will use the default window length and default step size that matches the window length
    defaultParams.stepSize = defaultParams.windowLength; %seconds
end
%     %no extra time
%     defaultParams.extraTime_pre = 0;%(2-params.windowLength)/2;
%     defaultParams.extraTime_post = 0;%(2-params.windowLength)/2;
%     
%     %pad to 2 seconds total evenly on right and left
%     defaultParams.extraTime_pre = (1-defaultParams.windowLength)/2;
%     defaultParams.extraTime_post = (1-defaultParams.windowLength)/2;
%     
%     %pad to 2 seconds total on only on left to simulate real-time (only history available, not future data)
%     defaultParams.extraTime_pre = (2-params.windowLength);
%     defaultParams.extraTime_post = 0;

%manual extra time
defaultParams.extraTime_pre = 0;
defaultParams.extraTime_post = 0;

defaultParams.cleanYBool = 1; %setting true may be better for AMICA later

defaultParams.pThres_source = 1; %1E-6 recommended
% defaultParams.rhoThres_source = []; %0.5-.99 recommended  (old)
defaultParams.rhoSqThres_source = .85;
                                       %note EMG parms below may be wrong (refer to R
                                       %instead of Rsq
                                        %1sec_8EMG->0.7; 2sec_8EMG->0.;
                                        %4sec_8EMG->; 8sec_8EMG->0.4;
                                        %16sec_8EMG->0.3; 32sec_8EMG->0.3;
                                        %100sec_8EMG->0.3:
                                        %+300sec_8EMG->0.25
                                        
                                        %rho sq suggested lower limits for Rsq based on
                                        %clean EEG phantom data with noise
                                        %electrodes
                                        %1sec_128N->[0.9,.95]; 2sec_128N->[0.75,0.85];
                                        %4sec_128N->[0.55,0.65]; 8sec_128N->[0.45];
                                        %16sec_128N->[0.3,0.4]; 32sec_128N->[0.2];
                                        %100sec_128N->[0.05,0.1]:
                                        %+300sec_128N->[.05,0.25]
% params.rhoThres_source_min = 0.5;

defaultParams.pThres_ch = 1; %1E-4 recommended
defaultParams.rhoThres_ch = 0; %0-1, as low as 0 probably ok

defaultParams.cleanXwith='X';
defaultParams.cleanYwith='Y';

defaultParams.calcCCAonWholeData = false; %ryan currently editing this feature

defaultParams.ccaCalcMethod = 0; %0 for matlab default, 1 for special version found online that claims to be better

defaultParams.noiseRemovalMethod = 0;

defaultParams.plotStatsOn = 1;
defaultParams.giveCleaningUpdates = 1;

if isempty(params) %if no params given, use all defaults
    params = defaultParams;
else %if some params given, update missing ones with default values as needed
    f = fields(defaultParams);
    for f_i = 1:length(f)
        if ~isfield(params,f(f_i)) %if missing expected default field
%             disp( f{f_i} )
            params.(f{f_i}) = defaultParams.(f{f_i}); %copy over default for missing field
        end
    end
end

%make sure step size isn't bigger than the cleaning window length
if params.stepSize > params.windowLength
params.stepSize = params.windowLength;
end