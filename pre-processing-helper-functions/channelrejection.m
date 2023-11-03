
function [clean_data, clean_data_timerej] = channelrejection(original_data,varargin)
 %% Reject bad EEG channels
 % [clean_data, clean_data_timerej] = channelrejection(original_data, Options...)
 % This function uses pop_rejchan (kurtosis and probability), trimOutlier (standard deviation). and
 % clean_artifacts (removes flatline channels, low-frequency drifts, noisy channels) for channel rejection. At least 2/4 of these methods have to suggest
 % removing a channel before it is removed from the original data. 
 
 % Inputs:
 %   orginal_data = Raw continuous EEG recording to clean up (as EEGLAB dataset structure), [channels x time]
 
 %NOTE: The following parameters are the core parameter options should be passed in as Name-Value Pairs.
 %
 %   ChannelRange =[chanStart:chanEnd] channel index that will be cleaned. 
 %                  Default: all channels
 %   KurtosisCriterion: Kurtosis criterion. Default: 5
 %   ProbabilityCriterion: Probability criterion. Default: 5
 %   StandardDeviationCriterion: TrimOutliers criterion. Default: 500
 %
 %######################################################################################################
 %clean_artifacts() parameters (Source: Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
 %                                2012-09-04)
 % 
 %   ChannelCriterion : Minimum channel correlation. If a channel is correlated at less than this
 %                      value to an estimate based on other channels, it is considered abnormal in
 %                      the given time window. This method requires that channel locations are
 %                      available and roughly correct; otherwise a fallback criterion will be used.
 %                      (default: 0.85)
 %
 %   LineNoiseCriterion : If a channel has more line noise relative to its signal than this value, in
 %                        standard deviations based on the total channel population, it is considered
 %                        abnormal. (default: 4)
 %
 %   BurstCriterion : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose
 %                    variance is larger than this threshold relative to the calibration data are
 %                    considered missing data and will be removed. The most aggressive value that can
 %                    be used without losing much EEG is 3. For new users it is recommended to at
 %                    first visually inspect the difference between the original and cleaned data to
 %                    get a sense of the removed content at various levels. A quite conservative
 %                    value is 5. Default: 'off'.
 %
 %   BurstRejection : 'on' or 'off'. If 'on' reject portions of data containing burst instead of 
 %                    correcting them using ASR. Default is 'off'.
 %
 %   WindowCriterion : Criterion for removing time windows that were not repaired completely. This may
 %                     happen if the artifact in a window was composed of too many simultaneous
 %                     uncorrelated sources (for example, extreme movements such as jumps). This is
 %                     the maximum fraction of contaminated channels that are tolerated in the final
 %                     output data for each considered window. Generally a lower value makes the
 %                     criterion more aggressive. Default: 0.25. Reasonable range: 0.05 (very
 %                     aggressive) to 0.3 (very lax).
 %
 %   Highpass : Transition band for the initial high-pass filter in Hz. This is formatted as
 %              [transition-start, transition-end]. Default: [0.25 0.75].
 %
 %   FlatlineCriterion : Maximum tolerated flatline duration. In seconds. If a channel has a longer
 %                       flatline than this, it will be considered abnormal. Default: 5
 %#####################################################################################################
 % Outputs:
 %   clean_data : EEG data after channel rejection
 %   clean_data_timerej : EEG data after channel and time window rejection
 %   
 %
 %
 %                              Noelle Jacobsen, Human Neuromechanics Lab,
 %                              University of Florida
 %                              2-6-2020

 hlp_varargin2struct(varargin,...
       {'kurt_crit','KurtosisCriterion'}, 5,...
       {'prob_crit','ProbabilityCriterion'}, 5 ,...
       {'std_crit','StandardDeviationCriterion'}, 500,...
       {'chan_range','ChannelRange'},1:128,...
       {'chan_crit1','ChannelCriterion'}, 0.85,...
       {'linenoise_crit', 'LineNoiseCriterion'}, 4,...
       {'flatline_crit1','FlatlineCriterion'},5,...
       {'highpass','Highpass'}, [0.25 0.75],...
       {'wind_crit', 'WindowCriterion'}, 0.25,...
       {'burst_crit1','BurstCriterion'}, 'off',...
       {'burst_rej','BurstRejection'}, 'off'); 
  
   %store parameter values in EEG.etc
   ETC = original_data.etc;
   ETC.clean_artifacts.parameters(1).crit = 'Kurtosis';
   ETC.clean_artifacts.parameters(1).value = kurt_crit;
   ETC.clean_artifacts.parameters(2).crit = 'Probability';
   ETC.clean_artifacts.parameters(2).value = prob_crit;
   ETC.clean_artifacts.parameters(3).crit = 'Standard Deviation';
   ETC.clean_artifacts.parameters(3).value = std_crit;
   ETC.clean_artifacts.parameters(4).crit = 'Channel Range';
   ETC.clean_artifacts.parameters(4).value = chan_range;
   ETC.clean_artifacts.parameters(5).crit = 'Channel Correlation';
   ETC.clean_artifacts.parameters(5).value = chan_crit1; 
   ETC.clean_artifacts.parameters(6).crit = 'Line Noise';
   ETC.clean_artifacts.parameters(6).value = linenoise_crit;
   ETC.clean_artifacts.parameters(7).crit = 'Flatline';
   ETC.clean_artifacts.parameters(7).value = flatline_crit1;
   ETC.clean_artifacts.parameters(8).crit = 'Highpass Filter';
   ETC.clean_artifacts.parameters(8).value = [0.25 0.75];
   ETC.clean_artifacts.parameters(9).crit = 'Window';
   ETC.clean_artifacts.parameters(9).value = wind_crit;
   ETC.clean_artifacts.parameters(10).crit = 'Burst';
   ETC.clean_artifacts.parameters(10).value = burst_crit1;
   ETC.clean_artifacts.parameters(11).crit = 'Burst Rejection';
   ETC.clean_artifacts.parameters(11).value = burst_rej;
   
 original_data = eeg_checkset(original_data);
 
 %Channel rejection using clean_artifact
 [clean_artifact_data HP BUR] = clean_artifacts(original_data, 'ChannelCriterion', chan_crit1,'LineNoiseCriterion',linenoise_crit,'FlatlineCriterion',flatline_crit1,'Highpass',[0.25 0.75],...
                        'WindowCriterion', wind_crit,'BurstCriterion', burst_crit1, 'BurstRejection', burst_rej); 
try                     
if size(clean_artifact_data.data,2) == size(BUR.data,2)
   fprintf('No time windows rejected \n');
else
 ETC.clean_artifacts.clean_sample_mask = clean_artifact_data.etc.clean_sample_mask; %store masks 
end
catch ME
    if (contains(ME.identifier, 'does not exist'))
        msg = ['clean_artifact_data.etc.clean_sample_mask does not exist'];
    end
end
 
 chan_range = 1:clean_artifact_data.nbchan; %update channel range to reflect channels removed from clean_artifact
 
 %Channel rejection using kurtosis
 clean_kurtosis_data = pop_rejchan(clean_artifact_data, 'elec',[chan_range],'threshold',kurt_crit,'norm','on','measure','kurt');
 %Channel rejection using probability
 EEG = eeg_checkset(clean_artifact_data);
 clean_proability_data = pop_rejchan(clean_artifact_data, 'elec',[chan_range],'threshold',prob_crit,'norm','on','measure','prob');
 %Channel rejection using TrimOutliers (uses standard deviation) 
 EEG = eeg_checkset(clean_artifact_data);
 disp('Check standard deviation, mark *bad* if channel std>500 \n')
 clean_std_data = trimOutlier(clean_artifact_data, -Inf, std_crit, Inf, 0);
 ETC.channelrejection.clean_drifts_kernel = clean_std_data.etc.clean_drifts_kernel; %store masks 
 if isfield(clean_std_data.etc,'clean_channel_mask')
 ETC.channelrejection.clean_channel_mask = clean_std_data.etc.clean_channel_mask;
 end
 %remove channels found in more than one rejection method
 %keep a 2 versions of channel rejection - one with time rejection
 %(clean_data_timerej) and one without time rejection (clean_data)
    clean_data = BUR;
    clean_data_timerej=clean_artifact_data;
    clean_data_timerej.data=[];
    clean_data_timerej.chanlocs=[];
    clean_data.data = [];
    clean_data.chanlocs= [];
for index1 = chan_range
    a = clean_artifact_data.chanlocs(index1).labels;
    count = 0; %the count indicates how many methods removed that channel. 0 = none removed it, 4 = all removed it.
    if isempty(find(strcmpi(a,{clean_kurtosis_data.chanlocs.labels}))) %the i in strcmpi means case insensitive when searching strings
        count = count + 1;
        fprintf('kurtosis suggested removing %i \n', index1)
    end
    if isempty(find(strcmpi(a,{clean_proability_data.chanlocs.labels})))
        count = count+1;
        fprintf('probability suggested removing %i \n', index1)
    end
    if isempty(find(strcmpi(a,{clean_artifact_data.chanlocs.labels})))
        count = count +2; %always remove channels cleaned by clean_artifact
        fprintf('clean_artifact suggested removing %i \n', index1)
    end
    if isempty(find(strcmpi(a,{clean_std_data.chanlocs.labels})))
        count = count+1;
        fprintf('standard deviation suggested removing %i \n', index1)
    end
    if count < 2 %If less than 2 of the methods removed the indicated channel, keep it
        clean_data.data = [clean_data.data;BUR.data(index1,:)];
        clean_data.chanlocs = [clean_data.chanlocs, BUR.chanlocs(index1)];
        clean_data_timerej.data = [clean_data_timerej.data;clean_artifact_data.data(index1,:)];
        clean_data_timerej.chanlocs = [clean_data_timerej.chanlocs,clean_artifact_data.chanlocs(index1)];
    end
end     
clean_data.etc = ETC;
clean_data_timerej.etc= ETC;
clean_data.nbchan = length(clean_data.chanlocs); %update number of channels after cleaning
clean_data_timerej.nbchan = length(clean_data_timerej.chanlocs);
num_chan_removed = original_data.nbchan - clean_data.nbchan; %find # of channels removed from original set
if num_chan_removed ~= 0 
    fprintf('%i channels were removed', num_chan_removed)
else 
    fprintf('No channels were removed')
end



