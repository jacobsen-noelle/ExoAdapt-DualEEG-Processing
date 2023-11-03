function [A,B,R,U,V,STATS] = iCC_calcSources(X,Y,params)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% convert to double (for rank computation)
X = double(X); Y = double(Y);

% %% keep track of original data (so we can clean it)
% X_orig = X; Y_orig = Y;

% %% re-reference temporarily (X to av X and Y to av Y)
%update: don't use this because it will cause errors if you send in raw EMG
%and raw noise data (different original references). Basically you won't
%delete the brain activity from the EMG channels because your avg re-ref is
%taking place across all channs (including noise) not just the EMG chans.
%This function knows nothing about what X and Y are
% %Note: re-ref prior to CCA is to avoid accidentally deleting brain activity 
% %For dual-layer EEG to cancel motion artifcats, this shouldn't be a concern.
% %However, for other signals such as neck EMG whose sensors \share a common 
% %recording reference with the EEG channels, not re-referencing prior to CCA
% %can lead to deletion of brain activity. 
% 
% %Note: To prevent accidentally deleting brain activity, we re-reference 
% %both the sets of channels (X and Y) to their individual selves. Note that
% %while full rank re-referencing is better from a rank perspective,
% %we have seen during testing that the rank it retains can be problematic 
% %when it comes to noise cancelling signals. Thus we recommend performing a 
% %'normal' average re-refererencing where a single rank is lost for X and Y.
% %If instead you want to use the full rank re-referencing approach (to try 
% %to delete more noise with fewer channels at the risk of deleting brain 
% %activity), we provide that option commented out below.
% 
% %Final note: since the re-referencing is applied in this function
% %locally, the deletion of sources happens elsewhere, we are able to
% %avoid changing the user's original reference scheme when we clean data.
% tempEEG = struct;
% tempEEG.data = X';%Note: column vectors  expected for X,Y in other parts of code but I have ready made code  for re-referencing  row vectors copied and modified from EEGLAB 
%     chansin = 1:size(tempEEG.data,1);
%     nchansin = length(chansin);
%     refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method,preferred)
%     % refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%full rank method
%     chansout = chansin;
%     tempEEG.data(chansout,:) = refmatrix*tempEEG.data(chansin,:);
% X = tempEEG.data'; clear tempEEG;
% 
% tempEEG.data = Y';%Note: column vectors  expected for X,Y in other parts of code but I have ready made code  for re-referencing  row vectors copied and modified from EEGLAB 
%     chansin = 1:size(tempEEG.data,1);
%     nchansin = length(chansin);
%     refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin; %normal avg ref (rank losing method,preferred)
%     % refmatrix = eye(nchansin)-ones(nchansin)*1/(nchansin+1);%full rank method
%     chansout = chansin;
%     tempEEG.data(chansout,:) = refmatrix*tempEEG.data(chansin,:);
% Y = tempEEG.data'; clear tempEEG;
%         


%% check rank and if  rank deficient, temporarily calulcate PCA scores (we will later clean original data using other functions)


Xrank = rank(X);
% Xrank = 16; %RYAN TEMP OVERIDE TO FORCE PCA TO TEST SOMETHING (DELETE!)
if Xrank < min(size(X)) %if rank deficient
    [~, SCORE] = pca(X); %calc PCA scores as alternative representaiton of data (only temporary)
    if size(SCORE,2) > Xrank %if we found more principal components than we think there is true data rank
        X = SCORE(:,1:Xrank); %take only the most important ones without going over rank
    else %if pca returned equal or less components than the orig number of channels, trust the pca algorithm's rank calculation instead
        X = SCORE;
    end
end

Yrank = rank(Y);
% Yrank =16; %warning('manual rank override');%RYAN TEMP OVERIDE TO FORCE PCA TO TEST SOMETHING (DELETE!)
if Yrank < min(size(Y))
    [~, SCORE] = pca(Y);
    if size(SCORE,2) > Yrank
        Y = SCORE(:,1:Yrank); %may be able to just set to SCORE since PCA automatically accounts for rank (input matrix of 10 columns and you can get back a matrix SCORES of less than 10 col)
    else
        Y = SCORE;
    end
end

%Note PCA SCORES are linear mixtures of the original channel data. CCA will
%later undo these mixtures as if we never converted to PCA. The main reason
%for all of this pca work is to avoid a warning output from matlab about 
%rank deficiency.

%% CCA find sources
if params.ccaCalcMethod == 1
    [theta,U,V,R,STATS] = subspacea(X,Y); A = []; B = []; %new fancy option that claims to be better (found online)
else
  [A,B,R,U,V,STATS] = canoncorr(X, Y );    %default matlab option
%    ncomp = 50; [A,B,U,V,BETA,PCTVAR,MSE,STATS] = plsregress(X,Y,ncomp); STATS.p = zeros(1,ncomp);  [tempRHO,tempPVAL] = corr(U,V); R =  diag(tempRHO)';%RYAN TEMP TEST with partial least squares
end

% % This little bit below is just a sanity check to see how well raw noise
% % ch data cleans EEG data rather than going thru the extra step of
% % calculating U,V. Looks like CCA is better
% U = Y(:,1:4:end);%ryan delete!!!!!!
% V = U;%ryan delete!!!!!!
% A = [];%ryan delete!!!!!!
% B = [];%ryan delete!!!!!!
% STATS.p = zeros(1,size(U,2));  [tempRHO,~] = corr(U,V); R =  diag(tempRHO)';%RYAN TEMP TEST delete

%% PCA option
%         %PCA find sources   
%         [coeff, U, latent, tsquared, explained] = pca(Y);
%         V = U; %for compatibility with code prev. programmed for cca
%         %badComps = 1:size(U,2); 
% %         badComps = 1:5; 
            
%% corr
% [R2,P2]=corr(X,V);
% % figure; plot((R2.^2)','o');
% myTempCorrFig = figure; 
% myTempCorrPlot = plot(1:size(V,2),(min(R2.^2)'),'o'); 
% hold on; 
% % plot(1:size(V,2),(mean(R2.^2)'),'o'); 
% plot(1:size(V,2),(median(R2.^2)'),'o'); 
% %plot(1:size(V,2),(max(R2.^2)'),'o');

% figure; plot(prctile(abs(R2),[90]),'o'); ylim([0 .3]);

end

