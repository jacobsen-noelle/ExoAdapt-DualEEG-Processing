function dataOut = AlignDAQ2DAQ(Fs_DAQa, t_DAQa,  t_DAQb, tEvents_DAQa, tEvents_DAQb, data_DAQb)
% Takes data from data acquisition device B (DAQb) and aligns/adds it to 
% data acquisition device A (DAQa), given timing information.
% 
% Fs_DAQa       -   Sampling frequency of DAQa (e.g., EEG)
% t_DAQa        -   time vector as recorded by DAQa
% t_DAQb        -   time vector as recorded by DAQb
% tEvents_DAQa  -   timing of events as recorded by DAQa
% tEvents_DAQb  -   timing of events as recorded by DAQb
% data_DAQb     -   matrix of data recorded from DAQb (e.g., LabVIEW). Each
%                   row vector is a separate signal. For example, a 5 x 200
%                   matrix is 5 signals with 200 samples each
% 
% Note: there needs to be at least two events in order to align the systems
% properly. I.E., tEvents_DAQa and tEvents_DAQb should be vectors of at
% least length = 2.
 
 
 % Author: 
%   Ryan Downey
%   Biomedical Engineering, University of Florida
%

%% Check that number ofevents matches
if length(tEvents_DAQa) ~= length(tEvents_DAQb)
    disp('ERROR: number of events does not agree between DAQa and DAQb');
    disp(['DAQa: ', num2str(length(tEvents_DAQa)),' and ','DAQb: ', num2str(length(tEvents_DAQb))]);
    return;
end

%% Check that the data matrix is channel x time, %Noelle added
if size(data_DAQb,1)>size(data_DAQb,2) 
    data_DAQb = data_DAQb';
end
    
%% Convert DAQb time to DAQa time based on lining up event markers 
t_DAQb_in_DAQa_base = interp1(tEvents_DAQb,tEvents_DAQa,t_DAQb,'linear','extrap'); %linear interpolation done here to avoid breaking causality (sample i+1 should always take place after sample i even  if the timing in general is bad for DAQb). Extrapolation done so we can time warp data outside the event markers
% Plot mapping between DAQa's local time and DAQb's
figure; plot(tEvents_DAQb,tEvents_DAQa,'o',t_DAQb,t_DAQb_in_DAQa_base,':.');
    title('Mapping between DAQa time and DAQb time');
    xlabel('DAQb local time (s)');
    ylabel('DAQa local time (s)');
    legend('Events','Interpolated');


%% Define the time points of DAQa's system that we are going to interpolate DAQb's values into 

% Note: since we are aligning from DAQb to DAQa, we need everything to  
% ultimately be represented in terms of exact sample numbers of DAQa.

% Take the beginning and end time points of DAQb (now represented in DAQa's
% reference frame but as decimal values of time that don't necessarily line 
% up perfectly with a sample of DAQa) and create a vector of time points 
% that are at intervals perfectly matching the sampling points of DAQa

% For example, say the first measurement of DAQb happened at 0.10034 seconds
% in terms of DAQa's time frame (i.e. t_DAQb,t_DAQb_in_DAQa_base(1) = 0.10034 
% and let's say for example that DAQa recorded at 1000 Hz. In this case,
% DAQa's samples are at 0, .001, .002, ..., .009, .100, .101, ... seconds.
% Thus there is no sampling instant of DAQa that exactly matches .10034 seconds 
% so we instead choose .100 (DAQa's sample #11) as the closest possible time point
% that perfectly fits DAQa's sampling. 

t1 = round(t_DAQb_in_DAQa_base(1)*Fs_DAQa)/Fs_DAQa; %start point of DAQb's recording in term of an exact sampling instant of DAQa
t2 = round(t_DAQb_in_DAQa_base(end)*Fs_DAQa)/Fs_DAQa; %end point of DAQb's recording in term of an exact sampling instant of DAQa
step = 1/Fs_DAQa;
t_DAQb_match_DAQa_sampling = [t1:step:t2]'; 

% Note that t_DAQb_match_DAQa_sampling is *NOT* all of DAQb's original time 
% points represented in DAQa's time frame at perfect samples. Rather this 
% variable is simply the first and last time points of DAQb represented in 
% DAQa's time frame at perfect sampling instants with perfect spacing 
% between them.

% Note that t1 and t2 can be outside of the actual range that DAQa recorded 
% (start earlier, end later, or both). That's not important for now. We will
% crop later. Obviously, there needs to be at least some overlap in the 
% "global" time that both systems recorded or it's pointless to time warp
% in the first place. 



%% Do the interpolation
% Note: interp1 throws an error if you send it an independent variable with 
% repeated numbers, so remove them 
%figure; plot(diff(t_DAQb_in_DAQa_base))
rejIndex = diff(t_DAQb_in_DAQa_base) == 0;
t_DAQb_in_DAQa_base = t_DAQb_in_DAQa_base(~rejIndex);
data_DAQb = data_DAQb(:,~rejIndex);
% figure; plot(diff(t_DAQb_in_DAQa_base))

data_interp = [];
data_interp = zeros(size(data_DAQb,1),length(t_DAQb_match_DAQa_sampling));
% disp(size(data_interp));
parfor row_i = 1:size(data_DAQb,1)
%    data_interp(row_i, :) =  interp1(t_DAQb_in_DAQa_base,data_DAQb(row_i,:),t_DAQb_match_DAQa_sampling,'spline');
   data_interp(row_i, :) =  interp1(t_DAQb_in_DAQa_base,data_DAQb(row_i,:),t_DAQb_match_DAQa_sampling,'linear','extrap');
end
% disp(size(data_interp));

% Initialize data output variable
dataOut = zeros(size(data_DAQb,1),length(t_DAQa)); %number of new signals X number of samples originally recorded with DAQa

% Trim edges if need be
index = t_DAQb_match_DAQa_sampling >= t_DAQa(1) & t_DAQb_match_DAQa_sampling <= t_DAQa(end);
t_crop = t_DAQb_match_DAQa_sampling(index); 
data_crop = data_interp(:,index);

if t_DAQb_match_DAQa_sampling(1) < t_DAQa(1)
    display('DAQb started recording before DAQa, cropping excess data at the beginning');
elseif t_DAQb_match_DAQa_sampling(1) > t_DAQa(1)
    display('DAQb started recording after DAQa, filling missing data with zeros');
end

if t_DAQb_match_DAQa_sampling(end) < t_DAQa(end)
    display('DAQb stopped recording before DAQa, filling missing data with zeros');
elseif t_DAQb_match_DAQa_sampling(end) > t_DAQa(end)
    display('DAQb continued to record after DAQa already stopped, cropping excess data at the end');
end



% Update dataOut matrix/vector with properly cropped data
startPoint = find(t_DAQa >= t_crop(1),1); %index of first DAQa time greater than beginning of DAQb recording
endPoint = startPoint+size(data_crop,2)-1;

if length(startPoint:endPoint) ~= length(data_crop) %make sure no issues with final alignment. This check probably not needed anymore. Was used during initial devleopment
    disp('ERROR with lining up DAQa and DAQb');
    disp(startPoint);
    disp(endPoint);
    disp(length(startPoint:endPoint))
    disp(length(data_crop))
    return;
else
     dataOut(:,startPoint:endPoint)=data_crop;
end



