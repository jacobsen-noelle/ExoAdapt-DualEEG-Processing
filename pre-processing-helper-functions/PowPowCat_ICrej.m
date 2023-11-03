% PowPowCat for IC rejection
% Reject ICs that have moderate cross-freq coupling in low (<8Hz) and high
% (>30Hz) frequency windows
% Takes the median correlation coefficient value in both high and low
% frequency windows, not including the identity coeffient (which always
% ==1)
%
% % Usage:
% >> [badPCC_IC] = PowPowCat_ICrej(EEG, varargin);
%
% Required inputs:
%   EEG            - EEG dataset structures with PowPowCat precomputed
%
% Optional inputs:
%   plotstuff      - plot extended component properties of "bad" ICs and
%                    spectral covariance (default = 0 (off))
%   outputFigureFolder - main folder path to store figures, subfolder will
%                       be created. (default = EEG.filepath)
% Output:
%   badPPC_IC      - Independent components identified as having moderate
%                    correlation in low and high freq. windows (e.g. eye
%                    blinks and muscle artifacts, respectively)
%
% Author: Noelle Jacobsen, University of Florida, 10/17/2021

function badPPC_IC = PowPowCat_ICrej(EEG,plotstuff, savePath)
% hlp_varargin2struct(varargin,...
%     {'plotstuff','plotstuff'}, 0,...
%     {'savePath','savePath'}, EEG.filepath );

covMatrix = EEG.etc.PowPowCAT.covMatrix;
freqs =  EEG.etc.PowPowCAT.freqs;
T = 0.3; %correlation coefficient threshold, 0.5-0.7= moderate correlation
badPPC_IC =[];
clear lowfreq_coupling highfreq_coupling
fprintf('\tIC#\t\tLowFreq Cov.\t\tHighFreq Cov.\n');
fprintf('________________________________________________');
for IC = 1:size(covMatrix,3)
    lowfreq_coupling(IC)= {covMatrix(find(freqs<8), find(freqs<8),IC)}; %below alpha
    highfreq_coupling(IC) = {covMatrix(find(freqs>30), find(freqs>30),IC)}; %above beta

    CC = lowfreq_coupling{1, IC};
    CC(logical(eye(size(CC)))) =[]; %remove diagnoal identity (coeff always equals one with itself)
    CC=reshape(CC,size(lowfreq_coupling{1, IC},1)-1,size(lowfreq_coupling{1, IC},1));
    lowfreq_coupling_avg = median(median(CC));

    CC = highfreq_coupling{1, IC};
    CC(logical(eye(size(CC)))) =[]; %remove diagnoal identity (coeff always equals one with itself)
    CC=reshape(CC,size(highfreq_coupling{1, IC},1)-1,size(highfreq_coupling{1, IC},1));
    highfreq_coupling_avg = median(median(CC));

    fprintf('\n\t%i\t\t%.2f\t\t\t\t%.2f',IC,lowfreq_coupling_avg,highfreq_coupling_avg)
    %identify ICs with corr. coeff. above thresh in low and high freq windows
    if lowfreq_coupling_avg > T || highfreq_coupling_avg >T
        badPPC_IC = [badPPC_IC,IC];
        fprintf('\t**BAD**');
    end
end



if plotstuff ==1

    for badIC = 1:length(badPPC_IC)
        fprintf('\nPlotting extended component properties of bad ICs');
        % plot extended comp properties
        pop_prop_extended(EEG,0,[badPPC_IC(badIC)],NaN,{'freqrange', [2 80]},{},1,'')
        savethisfig(gcf,strcat(EEG.subject,'_IC',num2str(badPPC_IC(badIC)),'prop'),[savePath,'\PowPowCAT'],'jpg'); close;
    end
    myfig=figure;
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
    set(gcf,'Color','w');
    mytitle= sgtitle([EEG.subject,' Spectral Covariance ']);
    for IC = 1:size(covMatrix,3)
        if IC <=36
            subplot(6,6,IC)
        elseif (IC == 37 || IC == 73 || IC== 109)
            savethisfig(myfig,strcat(EEG.subject,'_SpectralCovariance_',num2str(IC-1)),[savePath,'\PowPowCAT'],'jpg');
            myfig= figure;%create a new figure so you don't have all ICs on one plot
            set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
            set(gcf,'Color','w')
            mytitle= sgtitle({'Component projection to inferior head and neck electrodes',...
                extractBefore(EEG.filename,'.set')});
            mytitle.Interpreter = 'none';
        end

        imagesc(covMatrix(:,:,IC), [-0.8 0.8]);
        customColorMap = colormap(jet);
        colormap(customColorMap)
        currentAxesPosition = get(gca, 'position');
        colorbarHandle = colorbar;
        set(get(colorbarHandle, 'title'), 'String', 'Corr. Coef','fontsize',8)
        set(colorbarHandle, 'fontsize',8);
        axis xy
        axis square
        tickLabels = round(freqs(10:10:length(freqs))*10)/10;
        tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
        set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', tickLabels,...
            'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
            'fontsize', 5)
        xlabel('Frequency (Hz)', 'fontsize', 8)
        ylabel('Frequency (Hz)', 'fontsize', 8)
        if ismember(IC,badPPC_IC)
            title(['IC' int2str(IC)],'Color','r','fontsize',10)
        else
            title(['IC' int2str(IC)],'fontsize',10)
        end

    end
    savethisfig(myfig,strcat(EEG.subject,'_SpectralCovariance_',num2str(IC-1)),[savePath,'\PowPowCAT'],'jpg');
    fprintf('\nFigures stored here: %s\n',[savePath,'\PowPowCAT']);

    %plot bad ICs spectral covariance
    %subplot dimensions
    if length(badPPC_IC)>5
        M = 2;
        N = 5;
    else
        M = 1;
        N = length(badPPC_IC);
    end
    myfig = figure;
    set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[0 0 15 10]);
    set(gcf,'Color','w');
    mytitle= sgtitle([EEG.subject,' Spectral Covariance -Bad ICs']);
    for IC = 1:length(badPPC_IC)
        if length(badPPC_IC)>1
            subplot(M,N,IC)
        end
        imagesc(covMatrix(:,:,badPPC_IC(IC)), [-0.8 0.8]);
        customColorMap = colormap(jet);
        colormap(customColorMap)
        currentAxesPosition = get(gca, 'position');
        colorbarHandle = colorbar;
        set(get(colorbarHandle, 'title'), 'String', 'Corr. Coef','fontsize',8)
        set(colorbarHandle, 'fontsize',12);
        axis xy
        axis square
        tickLabels = round(freqs(10:10:length(freqs))*10)/10;
        tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
        set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', tickLabels,...
            'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
            'fontsize', 8)
        xlabel('Frequency (Hz)', 'fontsize',12)
        ylabel('Frequency (Hz)', 'fontsize', 12)
        title(['**IC' int2str(IC),'**'],'Color','r','fontsize',16)
    end
    savethisfig(myfig,strcat(EEG.subject,'_BadSpectralCovariance_',num2str(IC-1)),[savePath,'\Summary'],'jpg');
    fprintf('Summary of ICs with bad spectral covariance stored here: %s\n',[savePath,'\Summary']);
end
end

