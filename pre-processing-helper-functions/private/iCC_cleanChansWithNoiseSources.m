function [chData_clean] = iCC_cleanChansWithNoiseSources(chData_raw,noiseSources,noiseRemovalMethod, params)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% disp(['using removal method ',noiseRemovalMethod]);
X = chData_raw; %column vectors
V = noiseSources;
broadWindow = 1:size(X,1);
chData_clean = chData_raw;

%% Mean center the data
X_MC = detrend(X,'constant');%_MC = mean centered
V_MC = detrend(V,'constant');


%% define bad noise sources at source level (cancor stats)
% badComps = find(R>params.rhoThres_source & STATS.p < params.pThres_source); % find(R>0.98 & STATS.p < 1E-4) sometimes 1E-3
badComps = 1:size(V_MC,2); %temporary as we update code

%% clean data by finding noise projection onto dirty channels
%newest of new ready to speed away
chanNoiseReconstruct = zeros(size(X,2),length(broadWindow));
  
%    figure; imagesc(R.^2);

    
                                %using V
%                                         [tempRho, tempP ] = corr(X_MC(:,:),V_MC(:,badComps));
%                                         badCompsSubsetInd = find(max(abs(tempRho)) >= params.rhoThres_ch & min(tempP) <= params.pThres_ch);
                                        badCompsSubsetInd = badComps;
                                        fakeWinv = mrdivide(X_MC(:,:)',V_MC(:,badComps(badCompsSubsetInd))'); %possibly try replacing with V???
                                        chanNoiseReconstruct = chanNoiseReconstruct + fakeWinv*V_MC(:,badComps(badCompsSubsetInd))';
    
                                X_clean = X(:,:)-chanNoiseReconstruct';
                                 chData_clean(:,:) = X_clean;





% %previously newest and best and slowest
% for xCh_i = 1:size(X,2)
% %                             %Option 1: Remove all badcomps releveant to particular
% %                             %channel at once. Note this is very similar to other code in that it can
% %                             %overfit accidentally if mixing too many badcomponents onto a
% %                             %single channel. The advantage here is that it checks for
% %                             %subsets of U/V to remove from X rather than just grabbing the
% %                             %first however many badComps in the original sorted order).
% %                             %Essentially it allows a little better flexibility for
% %                             %removing different noise components from each channel rather
% %                             %than it being common. Thus it should be able to clean each
% %                             %channel a little better given the same number of bad
% %                             %components you are trying to remove.
% %     
% %                             chanNoiseReconstruct = zeros(1,length(broadWindow));
% %     
% %     
% %                                 %using V
% %                                         [tempRho, tempP ] = corr(X_MC(:,xCh_i),V_MC(:,badComps));
% %                                         badCompsSubsetInd = find(abs(tempRho) >= params.rhoThres_ch & tempP <= params.pThres_ch);
% %                                         %note line below was from when we were building up
% %                                         %chanNoiseReconstruct on an iterative basis. It's not needed here
% %                                         %but also doesn't really hurt so I'm leaving it
% %                                         fakeWinv = mrdivide(X_MC(:,xCh_i)',V_MC(:,badComps(badCompsSubsetInd))'); %possibly try replacing with V???
% %                                         chanNoiseReconstruct = chanNoiseReconstruct + fakeWinv*V_MC(:,badComps(badCompsSubsetInd))';
% %     
% %     
% %                                 X_clean = X(:,xCh_i)-chanNoiseReconstruct';
% % %                                 [lia1, locb1]  = ismember(window,broadWindow);
% % %                                 tempEEG.data(Xchans(xCh_i),window) = X_clean(locb1,:)';
% %                     chData_clean(:,xCh_i) = X_clean;
% %                             %End option 1
%     
%     
% %     %Option 2: Remove bad components channel by channel and one by
% %     %one. This should avoid the overfitting problem because we are
% %     %only projecting one bad component onto our cortical EEG data
% %     %at any given time (not allowing for mixtures of multiple bad
% %     %components at the same time). By doing multiple passes we
% %     %avoid overfitting and removing data we shouldn't. On the other
% %     %hand, it may increase the potential to accidentally add some random noise
% %     %that wasn't there originally.
% %     chanNoiseReconstruct = zeros(1,length(broadWindow));
% %     %             V_MC = U_MC; %ryan look here to fix
% %     %             [tempRho, tempP ] = corr(X_MC(:,xCh_i),U_MC(:,badComps));
% %     [tempRho, tempP ] = corr(X_MC(:,xCh_i),V_MC(:,badComps));
% %     %             tempcrap = find(abs(tempRho) > rhoThres_ch  & tempP < pThres_ch );
% %     %             display(length(tempcrap));
% %     for badComp_i = 1:length(badComps)
% %         %         tempCorr = corr(X_MC(:,xCh_i),U_MC(:,badComps(badComp_i)));
% %         %                 [tempRho, tempP ] = corr(X_MC(:,xCh_i),V_MC(:,badComps(badComp_i)));
% %         if ( abs(tempRho(badComp_i)) > params.rhoThres_ch  && tempP(badComp_i) < params.pThres_ch )
% %             %                 %U
% %             % %                 fakeWinv = mrdivide(X_MC(:,xCh_i)',U_MC(:,badComps(badComp_i))'); %possibly try replacing with V???
% %             %                     fakeWinv = mrdivide( (X_MC(:,xCh_i)'-chanNoiseReconstruct) ,U_MC(:,badComps(badComp_i))');
% %             %                     %note we are now building up chanNoiseReconstruct on an iterative basis
% %             %                     chanNoiseReconstruct = chanNoiseReconstruct + fakeWinv*U_MC(:,badComps(badComp_i))';%possibly try replacing with V???
% %             
% %             %V
% %             fakeWinv = mrdivide(X_MC(:,xCh_i)',V_MC(:,badComps(badComp_i))'); %possibly try replacing with V???
% %             %fakeWinv = mrdivide( (X_MC(:,xCh_i)'-chanNoiseReconstruct) ,V_MC(:,badComps(badComp_i))'); %possibly try replacing with V???
% %             %note we are now building up chanNoiseReconstruct on an iterative basis
% %             chanNoiseReconstruct = chanNoiseReconstruct + fakeWinv*V_MC(:,badComps(badComp_i))';%possibly try replacing with V???
% %             
% %             
% %         else
% %             %do nothing
% %         end
% %     end
% %     %     X_clean = X(:,xCh_i)-chanNoiseReconstruct'-mean(chanNoiseReconstruct,2)'; %best so far
% %     %     X_clean = X_clean - mean(X_clean); %remove the mean again cause why not
% %     %     tempEEG.data(EEG_chans(xCh_i),window) = X_clean';
% %     %             X_clean = X(:,xCh_i)-chanNoiseReconstruct'; %best so far
% %     %             tempEEG.data(EEG_chans(xCh_i),window) = X_clean';
% %     
% %     X_clean = X(:,xCh_i)-chanNoiseReconstruct';
% % %     [lia1, locb1]  = ismember(window,broadWindow);
% % % %     tempEEG.data(Xchans(xCh_i),window) = X_clean(locb1,:)';
% % % 
% % % chData_clean(window,xCh_i) = X_clean(locb1,:);
% % chData_clean(:,xCh_i) = X_clean;
% %     %End option 2
%     
% end


               
%         display(['Num bad comps = ',num2str(length(badComps))]);
        
        % %old overfitting prevention
        % badComps =1:5;
        % if length(badComps > length(window)/100)
        %     badComps = 1:1:length(window)/100;
        % end
        

        
        %         %new better method? (V may be better for evrything but EOG. V probably better for noise to prevent deletion of brain activity)
        %         %     fakeWinv = mrdivide(X(:,:)',U(:,badComps)'); %possibly try replacing with V???
        %         %     chanNoiseReconstruct = fakeWinv*U(:,badComps)';%possibly try replacing with V???
        
        %         fakeWinv = mrdivide(X_MC(:,:)',V_MC(:,badComps)'); %CCA mean cetners X and Y when it calculates U and V so we should calculate the fake inverse weights by looking at U or V's projection to X_MC rather than their projections to X)
        %         %     fakeWinv = mrdivide(X(:,:)',V(:,badComps)'); %possibly try replacing with V???
        %         chanNoiseReconstruct = fakeWinv*V_MC(:,badComps)';%possibly try replacing with V???
        %         X_clean = X(:,:)-chanNoiseReconstruct';
        % %         X_clean = X(:,:)-chanNoiseReconstruct'-mean(chanNoiseReconstruct,2)'; %best so far
        % %         X_clean = X_clean - mean(X_clean);
        %
        %         [lia1, locb1]  = ismember(window,broadWindow);
        %         tempEEG.data(Xchans,window) = X_clean(locb1,:)';
        %
        %         %clean noise for funsies
        %         fakeWinv = mrdivide(Y_MC(:,:)',V_MC(:,badComps)');
        %         chanNoiseReconstruct = fakeWinv*V_MC(:,badComps)';
        %         Y_clean = Y(:,:)-chanNoiseReconstruct';
        %         [lia1, locb1]  = ismember(window,broadWindow);
        %         tempEEG.data(Ychans,window) = Y_clean(locb1,:)';
        
        % %untested option
        % chNoiseReconstruct_test = (V_MC(:,badComps)*pinv(A(:,badComps)))';
        % V_rej = V_MC;
        % % V_rej(:,~badComps) = deal(0);
        % chNoiseReconstruct_test2 = (V_rej*pinv(A))';


end

