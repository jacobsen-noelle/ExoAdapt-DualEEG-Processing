% IC_rejection(EEG, outputFigureFolder, opt)
% Performs IC rejection using intersection of various criterion, outputs store EEG.etc.ic_criteria
% 1) PSD slope
% 2) IC source projection to EMG channels ( if EMG was included in ICA decomp
% 3) Residual variance
% 4) Dipole location using segmented MRI atlas
% 5) ICLabel class probability
%
% Usage: [EEG]= IC_rejection(EEG, outputFigureFolder, opt)
%
% Required inputs:
%   EEG                  - EEG dataset
%   outputFigureFolder   - Folder where you want to store figures from each
%                          rejection criteria
%
% Output:  
%   EEG                  -EEG dataset,criterion calculations stored in 
%                         EEG.etc.ic_criteria, components that are
%                         rejected at each step stored in 
%                         EEG.etc.comp_reject (1= reject, 0= keep)
%
% optional inputs
% opt.RVthresh 				- [0-1] residual variance threshold (default = 0.15)
% opt.PSD_slope_thresh		- (num) threshold for PSD slope (default = 0)
% opt.MRI_folder 		    - folder with subject segmented MRI (.mat containing EEG.subject and "seg" in ACPC coord space),or leave empty to use MNI template
%							   Uses dipole position stored in EEG.dipfit_ACPC.dip.pos, Could modify to use MNI coordsys atlas and dipole positions
% opt.ICLabel_brain_thresh  - [0-1] ICLabel brain probability threshold 

% Authors: 
%   Noelle Jacobsen,University of Florida
%   Makoto's useful EEGLAB code
%   Created 10/17/21
%   Last updated 3/14/23

function [EEG]= IC_rejection(EEG, outputFigureFolder,opt)
myEEGfilename = EEG.filename;
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
if ~isempty(EMG_chans) & size(EEG.icaact,1) == length([EEG_chans EMG_chans])
    EMGproj_flag = 1;
else
    EMGproj_flag = 0;
end

%check for inputs
if isfield(opt,'RVthresh')
	RVthresh=  opt.RVthresh;
else
RVthresh = 0.15;
end

if isfield(opt,'PSD_slope_thresh')
	PSD_slope_thresh = opt.PSD_slope_thresh; 
else
   PSD_slope_thresh = 0;
end

if isfield(opt,'MRI_folder')
	crit4flag = 1;
    MRI_folder = opt.MRI_folder;
	MRIfileList = dir(fullfile(MRI_folder));
	idx = find(strcmp({MRIfileList.name},[EEG.subject,'_seg_i.mat']));
	if isempty(idx)
		error('could not find segmented MRI file')
	else
	segmented_mri = MRIfileList(idx).name;
    end
      %setup fieldtrip
     addpath(opt.fieldtrip_path)   
     rmpath(genpath([opt.eeglab_path,'\plugins\Fieldtrip-lite20230716'])) % CHANGE folder; having more than one path to fieldtrip will confuse matlab, update to your eeglab fieldtrip lite folder
     ft_defaults;
	
else
	crit4flag = 0;
end

if isfield(opt, 'ICLabel_brain_thresh')
	brainLabelthresh = opt.ICLabel_brain_thresh; %IC label percent brain
else
	brainLabelthresh = 0.5; %IC label percent brain
end
	


%% Criteria 1 - PSD slope
%identify components whose PSD slope in specified freq window exceeds
%threshold
%calculate spectras
figure;
freq_factor = 4; %frequency resolution
frames =  size(EEG.icaact,2);
%Calcualate all component spectra
[spectra, freqs] = spectopo(EEG.icaact(:,:), frames, EEG.srate,'freqfac', freq_factor, 'percent', 20, 'icacomps', [], 'nicamaps', 5, 'freqrange',[2 100],'electrodes','off'); %ms of frame gotten from
close;

bad_comp_slope = [];
thresh = PSD_slope_thresh; %threshold for PSD slope criteria
minfreq = min(freqs);% frequency range minimum
maxfreq= 40; %Hz, frequency range max
myfreqs = find(freqs<=maxfreq);
fprintf('=======================================================================\n\t\tCriterion 1: Spectra Slope\n=======================================================================\n');
fprintf('Finding IC spectras with slopes > %i ...',thresh);
fprintf('Linear slope\n')
myfig= figure;%create a new figure so you don't have all ICs on one plot
set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
set(gcf,'Color','w');
xlabel('Frequency (Hz)')
ylabel('Log Power Spectral Density 10*log_10 (uV^2/Hz)')
mytitle= sgtitle({'Component power spectral density',...
    extractBefore(EEG.filename,'.set')});
mytitle.Interpreter = 'none';
for ic = 1:size(spectra,1)
    if ic <=36
        subplot(6,6,ic)
    elseif ic>36
        if (ic == 37 || ic == 73 || ic== 109)
            savethisfig(myfig,strcat(EEG.subject,'_ICspectra_',num2str(ic-1)),outputFigureFolder,'jpg');
            myfig= figure;%create a new figure so you don't have all ICs on one plot
            set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
            set(gcf,'Color','w');
            xlabel('Frequency (Hz)')
            ylabel('Log Power Spectral Density 10*log_10 (uV^2/Hz)')
            mytitle= sgtitle({'Component power spectral density',...
                extractBefore(EEG.filename,'.set')});
            mytitle.Interpreter = 'none';
        end
        if ic>36 && ic<=72
            subplot(6,6,(ic-36))
        elseif ic>72 && ic<=108
            subplot(6,6,(ic-72))
        elseif ic>108
            subplot(6,6,(ic-108))
        end
    end
    
    p = polyfit(freqs(myfreqs),spectra(ic,myfreqs)',1); % linear slope of spectra
    slope(ic) = p(1);
    fprintf('\nIC%i\t%.3f',ic,slope(ic));
    if slope(ic) > thresh
        fprintf(' **BAD**');
        bad_comp_slope = [bad_comp_slope, ic];
        plot(freqs(myfreqs),spectra(ic,myfreqs),'r','LineWidth', 2)
    else
        plot(freqs(myfreqs),spectra(ic,myfreqs),'LineWidth', 2)
    end
    hold on;
    f = polyval(p,freqs(myfreqs));
    plot(freqs(myfreqs),f,'k--')
    hold off;
    title(['IC' int2str(ic)])
    %ylim([min(min(spectra(:,myfreqs))) max(max(spectra(:,myfreqs)))])
    xlim([1 maxfreq]);
    if ic == 1 || ic== 37 || ic == 73 || ic ==109
        legend('PSD','Linear fit');
    end
end
savethisfig(myfig,strcat(EEG.subject,'_ICspectra_',num2str(ic-1)),[outputFigureFolder,'\PSD'],'jpg'); close;
fprintf('\nComponents with more positive spectra slopes:')
disp(bad_comp_slope)
%store all slopes
EEG.etc.ic_criteria.specslope = slope;

%plot components with bad spectras
myfig=figure;
if ~isempty(bad_comp_slope)
    for xx = 1:length(bad_comp_slope)
        i = bad_comp_slope(xx);
        subplot(ceil(sqrt(size(bad_comp_slope,2))),ceil(sqrt(size(bad_comp_slope,2))),xx)
        plot(freqs,spectra(i,:),'r','LineWidth', 2)
        hold on
        axis('tight')
        title(['IC' int2str(i)])
        xlabel('frequency (Hz)')
        ylabel('power (dB)')
        xlim([2 100]);
        %ylim([-55 -20]);
    end % ic plot
end
set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 9]);
set(gcf,'Color','w');
mytitle= sgtitle({'Components with spectra slopes >', num2str(thresh),...
    ' ',extractBefore(EEG.filename,'.set')});
mytitle.Interpreter = 'none';
savethisfig(myfig,strcat(EEG.subject,'_ICspectra_rejcomp'),[outputFigureFolder,'\PSD'],'jpg');
%find components with acceptable spectra slopes
goodIC_slope = setdiff([1:size(EEG.icaact,1)],bad_comp_slope);
EEG.etc.comp_reject.PSDslope= ones(size(EEG.icaact,1),1); 
EEG.etc.comp_reject.PSDslope(goodIC_slope) = 0;

%plot summary psd for good and bad ICs
figure;
[~, ~] = spectopo(EEG.icaact(goodIC_slope,:), frames, EEG.srate,'freqfac', freq_factor, 'percent', 20, 'icacomps', [], 'nicamaps', 5, 'freqrange',[2 100],'electrodes','off'); %ms of frame gotten from
title([EEG.subject,' PSD- Good ICs']);
savethisfig(gcf,strcat(EEG.subject,'_PSD_goodIC'),[outputFigureFolder,'\Summary'],'jpg');close;
if ~isempty(bad_comp_slope)
figure;
[~, ~] = spectopo(EEG.icaact(bad_comp_slope,:), frames, EEG.srate,'freqfac', freq_factor, 'percent', 20, 'icacomps', [], 'nicamaps', 5, 'freqrange',[2 100],'electrodes','off'); %ms of frame gotten from
title([EEG.subject,' PSD- Bad ICs']);
savethisfig(gcf,strcat(EEG.subject,'_PSD_badIC'),[outputFigureFolder,'\Summary'],'jpg');close;
end 
if EMGproj_flag ==1
%% Criteria 2 -IC source projection to EMG channels
fprintf('=======================================================================\n\t\tCriterion 2: IC Projection to EMG channels\n=======================================================================\n');
%list EEG channel labelslocated around the back lower part of head
inferior_chans = {'1-A11','1-A12','1-A13','1-A14','1-A24','1-A25','1-A26','1-A27','1-B8','1-B9','1-B10','1-B11','1-B14','1-D24','1-D31','1-D32'}; %labels of channels and the back bottom of head, nearest neck
inferior_chan_ind=[];
%find non-rejected posterior-inferior EEG channels in this dataset
for index = 1:length(inferior_chans)
    c = inferior_chans{index};
    if ~isempty(find(strcmpi(c,{EEG.chanlocs.labels})))
        x = find(strcmpi(c,{EEG.chanlocs.labels}));
        inferior_chan_ind= [inferior_chan_ind, x];
    end
end


[EEG_chans, EMG_chans] = getchantypes(EEG); %get EEG and EMG channel indices
mychans = [inferior_chan_ind EMG_chans]; %combine EEG and EMG channel index so we can look at projections to both
fprintf('Finding components that project highly onto EMG channels \n Calculating projection (z normalized)\n');
projection = EEG.icawinv; %projection of ICs onto channels [channels x components]
projection_norm2chan = (projection-mean(projection,2))./std(projection,0,2); %z normalization, take mean of all components projected onto each channel and divide by standard deviation
EEG.etc.ic_criteria.projection_norm2chan = projection_norm2chan;
myfig=figure;
%visualize projection of components onto channels
for i = 1:length(mychans)
    subplot(4,6,i)
    if i ==1;
        xlabel('Component')
        ylabel('Normalized Projection')
    end
    stem([1:size(projection_norm2chan,2)],projection_norm2chan(mychans(i),:));
    title(EEG.chanlocs(mychans(i)).labels)
end
set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 9]);
set(gcf,'Color','w');
mytitle= sgtitle({'Component projection to electrodes at back of head and neck ',...
    extractBefore(EEG.filename,'.set')});
mytitle.Interpreter = 'none';
savethisfig(myfig,strcat(EEG.subject,'_ICProjectionByChannel'),outputFigureFolder,'jpg');close;

%identify components that project
bad_proj = [];
rej_thresh = 0.50; %rejection criteria, components that project >50% to EMG channels
mychans = [inferior_chan_ind EMG_chans];
%visualize which channels each component projects to
myfig= figure;
set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);

mytitle= sgtitle({'Component projection to inferior and neck electrodes ',...
    extractBefore(EEG.filename,'.set')});
mytitle.Interpreter = 'none';
fprintf('ICs with >%i%% projection to EMG channels marked "bad"',rej_thresh*100);
fprintf('\tIC #\t EMG Projection(%%)\n\t___________________________\n');
for i = 1:size(projection_norm2chan,2)
    if i <=36
        subplot(6,6,i)
    elseif i>36
        if (i == 37 || i == 73 || i== 109)
            savethisfig(myfig,strcat(EEG.subject,'_ICProjectionByComponent_',num2str(i-1)),outputFigureFolder,'jpg');
            myfig= figure;%create a new figure so you don't have all ICs on one plot
            set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
            set(gcf,'Color','w')
            mytitle= sgtitle({'Component projection to inferior head and neck electrodes',...
                extractBefore(EEG.filename,'.set')});
            mytitle.Interpreter = 'none';
        end
        if i>36 && i<=72
            subplot(6,6,(i-36))
        elseif i>72 && i<=108
            subplot(6,6,(i-72))
        elseif i>108
            subplot(6,6,(i-108))
        end
    end
    stem([1:size(inferior_chan_ind,2)],projection_norm2chan(inferior_chan_ind,i)); hold on;
    start = size(inferior_chan_ind,2);
    stem([start+1:start+size(EMG_chans,2)],projection_norm2chan(EMG_chans,i),'m');
    title(['IC' int2str(i)])
    ylim([min(min(projection_norm2chan(:,:))) max(max(projection_norm2chan(:,:)))])
    %identify components that have high IC weight projection onto EMG
    %channels
    percent_proj = sum(abs(projection_norm2chan(EMG_chans,i)))/sum(abs(projection_norm2chan(:,i)));
    fprintf('\n\t%i\t\t\t%0.2f',i,percent_proj);
    if percent_proj >rej_thresh
        bad_proj = [bad_proj,i];
        fprintf('  **BAD**');
       title(['IC' int2str(i)],'Color','r')
       set(gca,'Color', [1 0.8 0.8]);
    else
       title(['IC' int2str(i)])
    end
    if i == 1 || i== 37 || i == 73 || i ==109
        xlabel('Channel')
        ylabel('ICA projection')
        lgd = legend('Inferior EEG channels','Neck EMG channels');
        lgd.Location = 'northwest';
    end
    hold off
end
lgd = legend('Inferior EEG channels','Neck EMG channels');
ylim([min(min(projection_norm2chan(:,:))) max(max(projection_norm2chan(:,:)))])
savethisfig(gcf,strcat(EEG.subject,'_ICProjectionByComponent_',num2str(i-1)),outputFigureFolder,'jpg');
%plot bad components
pop_viewprops( EEG, 0, [bad_proj], {'freqrange', [2 80]}, {}, 1, '' )

%plot bad component spectras
myfig=figure;
set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
set(gcf,'Color','w');
if ~isempty(bad_proj)
    for xx = 1:2:length(bad_proj)
        i = bad_proj(xx);
        subplot(2,round(length(bad_proj)/2),xx) %% change back to 8x17 for smaller datasets
        stem([1:size(inferior_chan_ind,2)],projection_norm2chan(inferior_chan_ind,i)); hold on;
        start = size(inferior_chan_ind,2);
        stem([start+1:start+size(EMG_chans,2)],projection_norm2chan(EMG_chans,i),'m');
        title(['IC' int2str(i) ' Projection'])
        xlabel('Channel');ylabel('ICA projection')
        xticks([1:length(mychans)]);
        xticklabels({EEG.chanlocs(mychans).labels});xtickangle(90);
        ylim([min(min(projection_norm2chan(:,:))) max(max(projection_norm2chan(:,:)))])
        if xx==1
            legend('Inferior EEG channels','Neck EMG channels');
        end
        subplot(2,round(length(bad_proj)/2),xx+1) %% change back to 8x17 for smaller datasets
        plot(freqs,spectra(i,:),'LineWidth', 2)
        axis('tight')
        title(['IC' int2str(i) ' PSD'])
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        xlim([2 100]);
        %ylim([-55 -20]);
    end % ic plot
end
sgtitle({'Components that project highly to EMG channels(>',num2str(rej_thresh*100),'%)'});
savethisfig(gcf,strcat(EEG.subject,'_ICProjectionByComponent_badICspec'),outputFigureFolder,'jpg');
goodIC_proj = setdiff([1:size(EEG.icaact,1)],bad_proj);
EEG.etc.comp_reject.EMGProjection = ones(size(EEG.icaact,1),1); 
EEG.etc.comp_reject.EMGProjection(goodIC_proj) = 0;
end
%% Criteria 3 - Residual variance
%  =======================================================================
%          Makoto's useful EEGLAB code
%  =======================================================================
% Perform IC rejection using residual variance of the IC scalp maps.
if RVthresh >100
	disp('adjusting RV threshold to max --> 1 (keep all)')
	RVthresh =1; %keep between 0-1
elseif RVthresh >1
	RVthresh = RVthresh/100; %scale between 0-1
end 
RVthresh = abs(RVthresh);
rvList    = [EEG.dipfit.model.rv];
goodRvIdx = find(rvList < RVthresh)'; % < 15% residual variance == good ICs.
badRvIdx = setdiff(1:length(rvList),goodRvIdx); % RV<15% == bad ICs
%  =======================================================================
%  =======================================================================
fprintf('=======================================================================\n\t\tCriterion 3: Residual Variance of IC Scalp Maps\n=======================================================================\n');
fprintf('Perform IC rejection using residual variance of the IC scalp maps\n');
fprintf('ICs with RV<15%%:\n');
disp(goodRvIdx);
EEG.etc.comp_reject.RV = ones(size(EEG.icaact,1),1); 
EEG.etc.comp_reject.RV(goodRvIdx) =0;

if EMGproj_flag==1 
%% plot extended comp properties of components with high EMG projection
% only dipoles with <15% RV were calculated, so only plot those with good
% RV
index = intersect(bad_proj, goodRvIdx);
for comp = 1:length(index)
pop_prop_extended(EEG,0,[index(comp)],NaN,{'freqrange', [2 80]},{},1,'')
savethisfig(gcf,strcat(EEG.subject,'_IC',num2str(index(comp)),'prop'),[outputFigureFolder,'\highEMGProj_CompProp'],'jpg'); close;
end
end
%% Criteria 4- Remove dipoles outside of the head
% % Option 1: ft_sourcedepth()
% % Note: ft_sourcedepth doesn't work for fem head model
% % Perform IC rejection using inside brain criterion.
%  =======================================================================
%          Makoto's useful EEGLAB code
%  =======================================================================
% vol = load(EEG.dipfit.hdmfile); % This returns 'vol'.
% dipoleXyz = zeros(length(goodRvIdx),3);
% for icIdx = 1:goodRvIdx
%     dipoleXyz(icIdx,:) = EEG.dipfit.model(goodRvIdx(icIdx)).posxyz(1,:);
% end
% depth = ft_sourcedepth(dipoleXyz, vol); %A negative depth indicates that the source is inside the source compartment, positive indicates outside.
% depthThreshold = 1;
% insideBrainIdx = find(depth<=depthThreshold);
%   ======================================================================
%   ======================================================================

%% Perform IC rejection using inside brain criterion.
if crit4flag
fprintf('=======================================================================\n\t\tCriterion 4: Inside brain (standard) or skull (custom head model)\n=======================================================================\n');
fprintf('Datasets with standard spherical dipfit already removed dipoles outside of head')
fprintf('Perform IC rejection using inside brain criterion\n');
fprintf('Using segemented MRI as atlas');
fprintf('\nLoading %s ...\n');
disp(segmented_mri)
cd(MRI_folder);
load(segmented_mri); %ACPC coordspace


%find ICs located outside of the skull based on segemented subject mri atlas
functional.coordsys = 'acpc';
%functional.coordsys = 'spm';
opt.queryrange = 3; % not sure what this param is, just copied value from ft_sourceplot param
opt.atlas = seg_i;
outsideSkullIdx = [];
insideSkullIdx = [];
for mydipole= 1:length(goodRvIdx)
    %opt.location= EEG.dipfit.model(goodRvIdx(mydipole)).posxyz'; %xyz position of dipole in MNI space
	opt.location= EEG.dipfit_ACPC.dip(goodRvIdx(mydipole)).pos'; %xyz position of dipole in ACPC space
    label = ft_atlas_lookup(opt.atlas, opt.location, 'coordsys', functional.coordsys, 'queryrange', opt.queryrange);
    EEG.etc.ic_criteria.atlaslabel(goodRvIdx(mydipole)) = {label};
    fprintf('Component %i label: ',goodRvIdx(mydipole));
    disp(label)
    %label output has multiple voters (7 labels total). Use majority vote
    vote = 0;
    if ~isempty(label) %if label is empty, it's outside of the head
        for voter = 1:length(label) % check what each "voter" classified the location as
            if ~contains(label{voter},'scalp') % as long as it's not scalp, count it as brain. Thought skull would be too close of a cutoff given spatial error
                vote = vote +1;
            end
        end
    end
    if vote >3 % if more than 3/7 labels (majority) classify it as skull or inside skull, keep this component
        insideSkullIdx = [insideSkullIdx, goodRvIdx(mydipole)];
        fprintf('\nIC %i Classification: Inside skull',goodRvIdx(mydipole));
    else
        outsideSkullIdx = [outsideSkullIdx,goodRvIdx(mydipole)];
        fprintf('\nIC %i Classification: Inside skull',goodRvIdx(mydipole));
    end
end
EEG.etc.ic_criteria.insideBrainIdx = ones(1,size(EEG.icaact,1));
EEG.etc.ic_criteria.insideBrainIdx(outsideSkullIdx) =0;
EEG.etc.comp_reject.outsideSkull = zeros(size(EEG.icaact,1),1); 
EEG.etc.comp_reject.outsideSkull(outsideSkullIdx) =1;

if ~isempty(outsideSkullIdx)
%plot summary
pop_dipplot( EEG, [outsideSkullIdx] ,'mri','C:\Users\jacobsen.noelle\Desktop\eeglab2022.0\plugins\dipfit\\standard_BEM\\standard_mri.mat','normlen','on');
title([EEG.subject,' Bad ICs']);
view([1,0,0])
savethisfig(gcf,strcat(EEG.subject,'_dipLoc_outsideSkull_saggital'),[outputFigureFolder,'\Summary'],'jpg');
view([0,-1,0])
savethisfig(gcf,strcat(EEG.subject,'_dipLoc_outsideSkull_coronal'),[outputFigureFolder,'\Summary'],'jpg');
view([0,0,1])
savethisfig(gcf,strcat(EEG.subject,'_dipLoc_outsideSkull_top'),[outputFigureFolder,'\Summary'],'jpg');
end

else
EEG.etc.ic_criteria.insideBrainIdx = ones(1,size(EEG.icaact,1));
end
insideSkullIdx = find(EEG.etc.ic_criteria.insideBrainIdx ==1);



%% Criteria 5 - ICLabel class
%Run IC Label on subset of data
%Needs ICLabel toolbox
if brainLabelthresh >100
	disp('adjusting brain label threshold to max --> 1 (keep all)')
	brainLabelthresh =1; %keep between 0-1
elseif brainLabelthresh >1
	brainLabelthresh = brainLabelthresh/100; %scale between 0-1
end 
brainLabelthresh = abs(brainLabelthresh);
fprintf('=======================================================================\n\t\tCriterion 5: ICLabel\n=======================================================================\n');
fprintf('Using IC label on subset of data\n');
tempEEG = pop_select( EEG, 'time',[1 300] );
tempEEG = iclabel(tempEEG);
brainComps = find(tempEEG.etc.ic_classification.ICLabel.classifications(:,1)>brainLabelthresh);
fprintf('Found %i ICs labeled >%i%% brain\n', length(brainComps), brainLabelthresh*100);
%pop_viewprops( tempEEG, 0, [brainComps'], {'freqrange', [2 80]}, {}, 1, '' )
EEG.etc.ic_classification = tempEEG.etc.ic_classification;
EEG.etc.comp_reject.ICLabel = ones(size(EEG.icaact,1),1); 
EEG.etc.comp_reject.ICLabel(brainComps) = 0;


%% Take AND across the five criteria.
fprintf('=======================================================================\n')
fprintf('Finding good ICs based on all four criteria:\n');
% goodIcIdx = intersect(goodIC_proj, goodRvIdx);
% goodIcIdx = intersect(goodIcIdx, goodIC_slope);
goodIcIdx = intersect(goodRvIdx, goodIC_slope);
goodIcIdx = intersect(goodIcIdx, insideSkullIdx);
goodIcIdx = intersect(goodIcIdx, brainComps);
disp(goodIcIdx')
badIcIdx = setdiff(1:size(EEG.icaact,1),goodIcIdx); %final bad component index
% store component labels
EEG.reject.gcompreject(badIcIdx)= 1; %1 = rejection, 0= keep
close all;
end


function savethisfig(fig,name,figfilepath,type)
if ~exist(figfilepath, 'dir') %check
    mkdir(figfilepath)
end
cd(figfilepath)
saveas(fig,name,type);
end

function [EEG_chans, EMG_chans] = getchantypes(EEG)
EEG_chans = find(strcmpi('EEG',{EEG.chanlocs.type}));
EMG_chans = find(strcmpi('EMG',{EEG.chanlocs.type}));
end


