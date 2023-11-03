%Load STUDY that's been clustered and has spec data
%written using EEGlab v14
function all_cond = plot_smooth_spec_from_study(STUDY,ALLEEG,savePath,clusters_to_plot,one_sub_per_cl,myplotParams,plotTopo)
addpath('C:\Users\jacobsen.noelle\Documents\GitHub\EEG_Processing\helper-functions\Plotting')
design = myplotParams.design;
plot_freq_sections =1;
plot_mov_bandxcond = 0;
notes = strcat(savePath,'\Notes.txt');
cd(savePath)
fopen(notes,'wt');
fprintf(notes,'%f',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond_labels = myplotParams.labels;
% ['on'|'off]
statsMethod ='perm';    % ['param'|'perm'|'bootstrap']
Alpha = 0.05;           % [NaN|alpha], Significance threshold (0<alpha<<1)
mcorrect = 'cluster'; %['fdr'|'holms'|'bonferoni'|'none'] correction for multiple
groupstats = 'off';
mode = 'fieldtrip';
plotstuff = 1;


%% load spec data from study
c =0;
CL_PSD = struct([]);
for CL = clusters_to_plot %cluster number one at a time
    c = c+1;
    cd(STUDY.filepath)

    %set stastical method
    condstats ='off';
    switch mode
        case 'eeglab'
            STUDY = pop_statparams(STUDY, 'condstats', condstats,'method',statsMethod,'alpha',Alpha,'mcorrect',mcorrect,'singletrials','off');
        case  'fieldtrip'
            STUDY = pop_statparams(STUDY, 'condstats', condstats,...
                'method',statsMethod,...
                'singletrials','off','mode',mode,'fieldtripalpha',Alpha,...
                'fieldtripmethod','montecarlo','fieldtripmcorrect',mcorrect);
    end
    [STUDY specdata specfreqs pgroup pcond pinter] = std_specplot(STUDY,ALLEEG,'clusters', CL);
    close;

    %% normalize using subject mean psd
    tmp_specdata = specdata;
    df = [];
    for cond = 1:size(tmp_specdata,1)
        %concatenate subject data across conditions
        if size(tmp_specdata{cond,1},2) < size(df,2) %if subject isn't in all conditions (e.g. A2 no errors) , then fill with NaN so specdata can be concatenated
            tmp_specdata{cond,1}(:,end+1:size(df,2)) = NaN;
        elseif size(tmp_specdata{cond,1},2) > size(df,2) && cond ~=1
            df(:,end+1:size(tmp_specdata{cond,1},2),:) = NaN;
        end
        df(:,:,cond) = tmp_specdata{cond,1};
    end
    %%
    sub_avg = nanmean(df,[1,3]); %subject mean across all conditions and freqs
    specdata_norm = cell(size(specdata));
    for condi = 1:size(tmp_specdata,1)
        psd = tmp_specdata{condi,1};
        psd_norm = psd-repmat(sub_avg,[size(specfreqs,1),1]);
        specdata_norm{condi,1} = psd_norm;
    end

    clear cond_spec cond_spec_norm
    width = 5; %size of smoothing window
    cond = [];

    for condi=1:size(specdata_norm,1)
        nancolumn = isnan(specdata_norm{condi,1}(1,:));
        specdata_norm{condi,1}(:,nancolumn) = [];
        specdata{condi,1} = filtfilt(1/width*(ones(width,1)),1,double(specdata_norm{condi,1})); %smooth data
        specdata{condi,1}(:,nancolumn) = NaN;
        cond(condi).name = {STUDY.design(design).variable(1).value{1, condi}};
    end

    %% weighted average across components in each cluster
    if one_sub_per_cl ~= 0
        design_subs = STUDY.design(design).cases.value;
        [~, rm_setidx] = setdiff({STUDY.datasetinfo.subject},design_subs);
        idx = ~ismember(STUDY.cluster(CL).sets, rm_setidx);
        CL_sets = STUDY.cluster(CL).sets(idx);
        unique_clus_subs = unique(CL_sets); %set index

        for i=1:size(specdata,1)
            cond(i).name = {STUDY.design(design).variable(1).value{1, i}};
            % cond(i).ersp = double(mean(allersp{i, 1},3))
            %conslidate subjects that appear more than once in a cluster
            for uc = 1:length(unique_clus_subs)
                x = find(CL_sets == unique_clus_subs(uc));
                CL_cond_spec = specdata{i, 1} ;%all spec data in cluster
                index = uc;
                if strcmp(cond(i).name,'SB2 late') %S19 SB2 late data was reject during preprocessing
                    s19_setnum = find(strcmp({STUDY.datasetinfo.subject},'S19'));
                    s19_i = find(CL_sets == s19_setnum );
                    temp = CL_sets;
                    if ~isempty(s19_i)
                        %find new index b/c S19 spec data is missing for this condition
                        temp([s19_i]) = [];
                        x = find(temp == unique_clus_subs(uc));
                        if uc>s19_i
                            index = index-1;
                        end
                    end
                    clear s19_i s19_setnum
                end
                if ~isempty(x)
                    sub_CL_cond_spec= CL_cond_spec(:,x);%all spec data in cluster belonging to a subject
                    if size(x,2)>1 %if subject appears more than once in cluster
                        disp('Subject IC ERSPs consolidated to')
                        if one_sub_per_cl ==1
                            cond(i).specdata(:,index) = nanmean(sub_CL_cond_spec,2);
                            disp('one subject/cluster by averaging')
                        elseif one_sub_per_cl ==2
                            cond(i).specdata(:,index) = sub_CL_cond_spec(:,1);
                            disp('one subject/cluster using lowest IC number')
                        end
                    else
                        cond(i).specdata(:,index) = sub_CL_cond_spec;
                    end
                end
            end

            if ~strcmp(cond(i).name,'SB2 late')
                if (size(cond(i).specdata,2) ~= length(unique_clus_subs))
                    error('Dimensions of spec data do not match number of unique subjects in cluster')
                end
            elseif (strcmp(cond(i).name,'SB2 late')&&(size(cond(i).specdata,2) ~= size(unique(temp),2)))
                error('Dimensions of spec data do not match number of unique subjects in cluster')
            end
            cond(i).avg_psd = double(nanmean(cond(i).specdata(:,:),2));
            specdata{i,1}= cond(i).specdata;
        end
    else
        for i=1:size(specdata,1)
            cond(i).avg_psd = double(nanmean(specdata{i, 1},2));
            cond(i).specdata = specdata{i,1};
        end
    end

    %% reorder conditions
    mylegend = [STUDY.design(design).variable(1).value];
    changePlotOrder =1;
    if changePlotOrder ==1
        %plot order, write condition names in order you want them plotted
        StringNames = myplotParams.labels;
        mylegend = myplotParams.legend;
        order= [];
        reorderedData_avgpsd =[];
        reorderedData_allpsd = {};
        for s = 1:length(StringNames)
            try
                order(1,s) = find(strcmp(StringNames{1,s},[STUDY.design(design).variable(1).value]));

            catch
                disp('Variable name not found in STUDY.design.variable')
            end
            reorderedData_avgpsd(s,:)= [cond(order(1,s)).avg_psd]';
            reorderedData_allpsd{s,1} = cond(order(1,s)).specdata;
        end
        specdata_all = reorderedData_allpsd;
        specdata_avg = reorderedData_avgpsd;

    else
        specdata_avg = [cond.avg_psd]';
        specdata_all = {cond.specdata}.';
    end

    %% Stats
    % set parameters
    % -------------
    condstats = 'on' ;
    subject = '';
    comps = [];
    %set stastical method
    switch mode
        case 'eeglab'
            STUDY = pop_statparams(STUDY, 'condstats', condstats,'method',statOpt.statsMethod,'alpha',Alpha,'mcorrect',mcorrect,'singletrials','off');
        case  'fieldtrip'
            STUDY = pop_statparams(STUDY, 'condstats', condstats,...
                'method',statsMethod,...
                'singletrials','off','mode',mode,'fieldtripalpha',Alpha,...
                'fieldtripmethod','montecarlo','fieldtripmcorrect',mcorrect,'fieldtripnaccu',10000);
    end
    if length(comps) == 1
        stats.condstats = 'off'; stats.groupstats = 'off';
        disp('Statistics cannot be computed for single component');
    end

    stats = STUDY.etc.statistics;
    stats.fieldtrip.channelneighbor = struct([]); % assumes one channel or 1 component
    if isempty(STUDY.design(STUDY.currentdesign).variable)
        stats.paired = { };
    else
        stats.paired = { STUDY.design(STUDY.currentdesign).variable(:).pairing };
    end


    [pcond, pgroup, pinter, statscond, statsgroup, statsinter, pval, F] = std_stat_clusterpval(specdata_all, stats);
    cluster_perm_test(1).name = 'all_cond';
    cluster_perm_test(1).freqrange = STUDY.etc.specparams.freqrange;
    cluster_perm_test(1).pval = pval{1,1};
    cluster_perm_test(1).pcond = pcond;

    %% if there is a sig effect of condition on psd, then test specific
    sigcondeffect = any(pval{1,1}< Alpha);

    if sigcondeffect
        %hypotheses pairs--- pre-adapt vs adapt.
        %         ref_cond = 'A2-correct';
        %         test_cond = 'A2-error';
        ref_cond = myplotParams.refCond;
        test_cond = myplotParams.testCond;
%         ref_cond = 'SB2 early';
%         test_cond = 'SB1 early';
        refi = find(strcmpi(StringNames,ref_cond));
        testi = find(strcmpi(StringNames,test_cond));
        conditions2compare = [refi testi];

        %test specific hypotheses
        if CL == 7 || CL ==10 || CL ==3 || CL ==13 %ACC
            % theta
            stats_freq_range = [4 7];
            freqi = find(specfreqs>=stats_freq_range(1) & specfreqs<=stats_freq_range(2));
            specdata_all_tmp = {};
            for condi = 1:2
                specdata_all_tmp{condi,1} = specdata_all{conditions2compare(condi),1}(freqi,:); %selected frq range
            end
            [pcond, pgroup, pinter, statscond, statsgroup, statsinter, pval] = std_stat_clusterpval(specdata_all_tmp, stats);

            pmask = find(pcond{1,1}==0); % 0 = n.s.
            pval = pval{1,1};
            pval(pmask) =1;%apply mask?

            if any(pval<0.05) %if any pvalues <alpha, continue to next step
                %determine effect size for cluster with lowest p-val
                [effect] = calc_clust_effectsize(specdata_all_tmp,specfreqs(freqi), pval,0);

                cluster_perm_test(end+1).name = 'theta';
                cluster_perm_test(end).freqrange = stats_freq_range;
                cluster_perm_test(end).pval = pval;
                cluster_perm_test(end).pcond = pcond{1,1};
                cluster_perm_test(end).effect = effect;
            end

        end

        if CL == 6 || CL ==8 || CL ==14 %PPC, SMI
            % alpha
            stats_freq_range = [8 12];
            freqi = find(specfreqs>=stats_freq_range(1) & specfreqs<=stats_freq_range(2));
            specdata_all_tmp = {};
            for condi = 1:2
                specdata_all_tmp{condi,1} = specdata_all{conditions2compare(condi),1}(freqi,:); %selected frq range
            end
            [pcond, pgroup, pinter, statscond, statsgroup, statsinter, pval] = std_stat_clusterpval(specdata_all_tmp, stats);
            pmask = find(pcond{1,1}==0); % 0 = n.s.
            pval = pval{1,1};
            pval(pmask) =1;%apply mask?
            if any(pval<0.05) %if any pvalues <alpha, continue to next step
                [effect] = calc_clust_effectsize(specdata_all_tmp,specfreqs(freqi), pval,0);
                cluster_perm_test(end+1).name = 'alpha';
                cluster_perm_test(end).freqrange = stats_freq_range;
                cluster_perm_test(end).pval = pval;
                cluster_perm_test(end).pcond = pcond{1,1};
                cluster_perm_test(end).effect = effect;
            end
            clear pmask pval pcond effect
            stats_freq_range = [13 30];
            freqi = find(specfreqs>=stats_freq_range(1) & specfreqs<=stats_freq_range(2));
            specdata_all_tmp = {};
            for condi = 1:2
                specdata_all_tmp{condi,1} = specdata_all{conditions2compare(condi),1}(freqi,:); %selected frq range
            end
            [pcond, pgroup, pinter, statscond, statsgroup, statsinter, pval] = std_stat_clusterpval(specdata_all_tmp, stats);
            pmask = find(pcond{1,1}==0); % 0 = n.s.
            pval = pval{1,1};
            pval(pmask) =1;%apply mask?
            if any(pval<0.05) %if any pvalues <alpha, continue to next step
                [effect] = calc_clust_effectsize(specdata_all_tmp,specfreqs(freqi), pval,0);
                cluster_perm_test(end+1).name = 'beta';
                cluster_perm_test(end).freqrange = stats_freq_range;
                cluster_perm_test(end).pval = pval;
                cluster_perm_test(end).pcond = pcond{1,1};
                cluster_perm_test(end).effect = effect;
                clear pmask pval pcond effect
            end
        end
    end
    %% plot PSD
    % eeglab plot with significant stats bar on bottom
    % mask array is a binary mask created with std_specplot when stats
    % parameters are set.
    if plotstuff ==1

        try
            label = cellstr(STUDY.cluster(CL).label);
            label = label{1,1};
            mytitle= [strcat(myplotParams.title,num2str(CL),'-',label)];
        catch
            mytitle= [strcat(myplotParams.title,num2str(CL))];
        end
        %colors = [{lightOrange},{red},{lightRed},{purple},{darkPurple}, {lightPurple},];
        colors = myplotParams.colors;

        plotopt = {'highlightmode','bottom','plotmean','off','ylim',[], 'xlabel','Frequency (Hz)','ylabel',...
            'Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)','legend',mylegend};
        myspecplot = figure('InvertHardcopy', 'off', 'PaperType', 'a2', 'PaperOrientation', 'landscape');
        
        %shaded error bars- TEMP
        myrange= [];
        for xx = 1:size(specdata_all,1)
        shadedErrorBar(specfreqs, specdata_all{xx,1}', {@mean,@std}, 'lineprops',  {'-','Color',colors{1,xx},'LineWidth',2},'patchSaturation',0.1)
        myrange(xx,1) = min(specdata_avg(xx,:)- std(specdata_all{xx,1},0,2)');
        myrange(xx,2) = max(specdata_avg(xx,:)+ std(specdata_all{xx,1},0,2)');
        hold on;
        end

        plotcurve_colors( specfreqs,specdata_avg, 'colors', colors, 'maskarray', cluster_perm_test(1).pcond{1,1}', plotopt{1:end}, 'title', mytitle);

        %set figure settings
        set(gcf,'Color','w')
        set(gcf,'PaperUnits','Inches','Units','Inches','PaperPosition',[1 0 8 4.2],'Position', [1 0 8 4.2]);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        pos = get(gca,'Position');
        xticklabel = get(gca,'XTickLabel');
        yticklabel = get(gca,'YTickLabel');
        ytick = get(gca,'yTick');
        ylim([min(myrange(:,1)) max(myrange(:,2))])
        %set(gca,'Visible','off')
        % axes('Position',pos,'XAxisLocation','bottom','YAxisLocation','left',...
        %  'Color','none','XTickLabel',xticklabel,'YTickLabel',yticklabel,...
        %  'XColor','k','YColor','k')
        %axes('Position',pos,'XAxisLocation','bottom','YAxisLocation','left')
        ylabel({'Log Power'; 'Spectral Density ';'10*log_{10}(\muV^{2}/Hz)'});
        ylh = get(gca,'ylabel');ylp = get(ylh, 'Position');
        %         set(ylh,'FontSize',12)
        set(ylh, 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',12);
        set(ylh,'Position',[ylp(1)-3 ylp(2) ylp(3)])
        xlh = get(gca,'xlabel'); set(xlh,'FontSize',12);

        % add dashed vertical lines for freq bands of interest
        if plot_freq_sections ==1
            hold on;
            yl = get(gca,'YLim');
            if STUDY.etc.specparams.freqrange(1)>=4
                plot([4 4],[yl(1) yl(2)],':k')
                text(5,yl(2)-1,'\Theta','FontSize',12)
                plot([8 8],[yl(1) yl(2)],':k')
            end
            text(9,yl(2)-1,'\alpha','FontSize',12)
            plot([13 13],[yl(1) yl(2)],':k')
            text(14,yl(2)-1,'\beta','FontSize',12)
            if STUDY.etc.specparams.freqrange(2) >30
                plot([30 30],[yl(1) yl(2)],':k')
                text(32.5,yl(2)-1,'\gamma','FontSize',12)
            end
            legend off
            legend(mylegend)
        end

        lgd = findobj(gcf, 'Type', 'Legend'); lgd.FontSize = 7;
        %         set(lgd,'Position',[pos(1)-150 pos(2)+90 pos(3) pos(4)]);
        %         set(lgd, 'EdgeColor','none','Orientation','horizontal')
        set(lgd, 'EdgeColor','none','Orientation','horizontal')
        pos = get(lgd, 'Position');
        set(lgd,'Position',[0 0.95 pos(3) pos(4)]);
        box off


        %save figure
        savethisfig(myspecplot,['PSD_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(STUDY.etc.specparams.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg'],'jpg')
        savethisfig(myspecplot,['PSD_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(STUDY.etc.specparams.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.fig'],[savePath,'\PSD\',myplotParams.figname,'\fig'],'fig')
        savethisfig(myspecplot,['PSD_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(STUDY.etc.specparams.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.svg'],[savePath,'\PSD\',myplotParams.figname,'\svg'],'svg')
        if ~exist([savePath,'\PSD\',myplotParams.figname,'\dpdf\'], 'dir') %check
            mkdir([savePath,'\PSD\',myplotParams.figname,'\dpdf\'])
        end
        print([savePath,'\PSD\',myplotParams.figname,'\dpdf\PSD_',myplotParams.figname,num2str(CL),'_',num2str(round(max(specfreqs))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.dpdf'], '-dpdf', '-painters'); % Makoto's print method. On Linux.
        close;

        %% plot condition A vs B stats
        if sigcondeffect
            conditionIndex = [refi testi];
            for perm_testi = 2:size(cluster_perm_test,2) %skip the first all condition cluster-test-- just do results from individual freq band tests
                plot_psd_comparison(STUDY, CL, specfreqs,specdata_avg,statsMethod,cluster_perm_test(1,perm_testi),conditionIndex, myplotParams, mylegend,savePath)

            end
        end

     
    end

    all_cond(c).cluster = CL;
    all_cond(c).label = STUDY.cluster(CL).label;
    all_cond(c).data.allpsd = specdata_all;
    all_cond(c).data.avgpsd = specdata_avg;
    all_cond(c).data.freqs = specfreqs;
    all_cond(c).cluster_perm_test = cluster_perm_test; clear cluster_perm_test


    %%  plot moving avg band power across conditions
    if plot_mov_bandxcond ==1
        cd(STUDY.filepath)
        clust_spec = struct([]);
        x = 1;
        for seti = 1:size(STUDY.cluster(CL).comps,2)
            setnum = STUDY.cluster(CL).sets(seti);
            [STUDY, specdata, freqs, yvals, events, fileparams] = std_readdata(STUDY, ALLEEG, 'clusters', CL, 'freqrange', [2:100], ...
                'component', seti,'subject',STUDY.datasetinfo(setnum).subject, 'singletrials', 'on', 'design', myplotParams.design, 'datatype', 'spec');
            if ~isempty(specdata{1,1}) %specdata will be empty if this subject isn't included in the design
                %normalize using subject mean psd of all conditions and all freqs as baseline
                specdata_norm = cell(size(specdata));
                all_spec = [];
                for cond = 1:size(specdata,1)
                    if ~isempty(specdata{cond,1}) %some subjects don't have SLA errors during A2
                        all_spec(:,cond) = mean(specdata{cond,1},2);
                    else
                        all_spec(:,cond) = NaN;
                    end
                end
                avg_spec = nanmean(all_spec,[1,2]); %take average across conditions and frequency bands

                for cond = 1:size(specdata,1)
                    cond_spec = specdata{cond,1};
                    specdata_norm{cond,1} = specdata{cond,1}-avg_spec;
                end

                %filling missing strides in each cond with average
                for i = 1:size(specdata,1)
                    nanflag = 0;
                    if contains(STUDY.design(myplotParams.design).variable(1).value{1,i}, 'perturbation')
                        if size(specdata{i,1},2) < 10 %10 strides in perturbation cond
                            specdata{i,1}(:,end:10) = repmat(mean(specdata{i,1},2),[1,10-size(specdata{i,1},2)+1]);
                        end
                    elseif size(specdata{i,1},2) < 30 %30 strides in all other conditions
                        if isempty(specdata{i,1})
                            specdata{i,1} = repmat(nan,[length(freqs),30]);
                            nanflag =1;
                        else
                            specdata{i,1}(:,end:30) = repmat(mean(specdata{i,1},2),[1,30-size(specdata{i,1},2)+1]);
                        end
                    elseif size(specdata{i,1},2) > 30
                        specdata{i,1} = specdata{i,1}(:,1:30);
                    end

                    if ~nanflag
                        clust_spec(i).spec(:,:,x) = specdata{i,1};
                        %              clust_spec(i).spec_norm(:,:,sub) = specdata{i,1}-avg_spec;
                        %              %subtract mean spec across all conditions (subject's grand
                        %              mean) % i think subject mean is removed if it's turned on in
                        %              STUDY.etc.specparams
                        clust_spec(i).spec_smooth(:,:,x) = filtfilt(1/10*(ones(10,1)),1,double(clust_spec(i).spec(:,:,x))); %smooth data
                    else
                        clust_spec(i).spec(:,:,x) = specdata{i,1};
                        clust_spec(i).spec_smooth(:,:,x) =specdata{i,1};
                    end
                    clust_spec(i).condition = STUDY.design(myplotParams.design).variable(1).value{1,i};

                end
                x = x+1;
                clear all_spec avg_spec
            end
        end
        cd(savePath)

        %% reorder conditions
        mylegend = [STUDY.design(myplotParams.design).variable(1).value];
        changePlotOrder =1;
        if changePlotOrder ==1
            %plot order, write condition names in order you want them plotted
            StringNames = myplotParams.labels;
            mylegend = myplotParams.legend;
            order= [];
            reorderedData_allpsd = struct([]);
            for s = 1:length(StringNames)
                try
                    order(1,s) = find(strcmp(StringNames{1,s},[STUDY.design(myplotParams.design).variable(1).value]));

                catch
                    disp('Variable name not found in STUDY.design.variable')
                end

                reorderedData_allpsd(s).spec = double(clust_spec(order(1,s)).spec);
                reorderedData_allpsd(s).spec_smooth = clust_spec(order(1,s)).spec_smooth;
                reorderedData_allpsd(s).condition =  clust_spec(order(1,s)).condition;
            end
            clust_spec = reorderedData_allpsd;
        end

        % get mean spec across all conditions
        allspec = [];
        for i = 1:size(clust_spec,2)
            allspec = [allspec, clust_spec(i).spec_smooth];
        end
        %%
        allspec_mean = squeeze(nanmean(allspec,2)); %avg spec across all conditions


        [allspec_mean bad_sub_ind] = rmoutliers(allspec_mean',"mean"); %If input is a matrix, then rmoutliers detects outliers in each column of A separately and removes the entire row; method "mean":efines an outlier as an element of A more than three standard deviations from the mean.
        allspec_mean = allspec_mean';

        %temp, uncommment later
        %     myspecplot = figure;
        %     plot(freqs,allspec_mean);
        %     title ({[label,'-CL ',num2str(CL)]},{'Condition average PSD by subject'},'interpreter','none')
        %     ylabel('Power 10*log_1_0 (uV^2/Hz)');
        %     xlabel('Frequency (Hz)');
        %     set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[1 1 6.5 4.5],'Color','w');
        %     set(gca,'Fontsize',14,'FontName','Arial');
        %     xlim([2 100])
        %     savethisfig(myspecplot,['CondAvgPSD_bySub',myplotParams.figname,num2str(CL),'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg\'],'jpg')
        %     close;


        diary on
        fprintf('\n%i outlier subjects (>3std from mean) removed',size(find(bad_sub_ind==1),1))
        diary off
        %find individual center frequency in each band
        pk_alpha = []; pk_beta =[]; pk_theta = []; pk_gamma= [];
        for i=1:size(allspec_mean,2) %loop through each subject
            data = allspec_mean(:,i);
            beta =  data(freqs>=13 & freqs<30,:,:);
            theta = data(freqs>=4 & freqs<7,:,:);
            alpha = data(freqs>=8 & freqs<13,:,:);
            gamma = data(freqs>=30 & freqs<128,:,:);

            [pks,locs] = findpeaks(theta,freqs(freqs>=4 & freqs<7),'SortStr','descend');
            if ~isempty(locs)
                pk_theta(i) = locs(1); % peak alpha freq is first loc in list
            else
                pk_theta(i) = NaN; %replace NaN with group avg center freq later
            end

            [pks,locs] = findpeaks(alpha,freqs(freqs>=8 & freqs<13),'SortStr','descend');
            if ~isempty(locs)
                pk_alpha(i) = locs(1); % peak freq is first loc in list
            else
                pk_alpha(i) = NaN;
            end

            [pks,locs] = findpeaks(beta,freqs(freqs>=13 & freqs<30),'SortStr','descend');
            if ~isempty(locs)
                pk_beta(i) = locs(1); % peak freq is first loc in list
            else
                pk_beta(i) = NaN;
            end

            [pks,locs] = findpeaks(gamma,freqs(freqs>=30 & freqs<128),'SortStr','descend');
            if ~isempty(locs)
                pk_gamma(i) = locs(1); % peak freq is first loc in list
            else
                pk_gamma(i) = NaN;
            end
        end

        try
            label = cellstr(STUDY.cluster(CL).label);
            label = label{1,1};
            mytitle= [strcat(myplotParams.title,num2str(CL),'-',label)];
        catch
            mytitle= [strcat(myplotParams.title,num2str(CL))];
        end
        diary on
        fprintf('\n\nCL-%i-%s\n',CL,label)
        fprintf('Freq. band\t\tMean Center Freq. (Hz)\t Std Center Freq. (Hz)\t\n =========================================================\n')
        fprintf('\t%s\t\t\t %0.2f\t\t\t\t\t\t%0.2f\n','theta',nanmean(pk_theta),nanstd(pk_theta))
        fprintf('\t%s\t\t\t %0.2f\t\t\t\t\t\t%0.2f\n','alpha',nanmean(pk_alpha),nanstd(pk_alpha))
        fprintf('\t%s\t\t\t %0.2f\t\t\t\t\t\t%0.2f\n','beta',nanmean(pk_beta),nanstd(pk_beta))
        fprintf('\t%s\t\t\t %0.2f\t\t\t\t\t\t%0.2f\n','gamma',nanmean(pk_gamma),nanstd(pk_gamma))
        diary off

        %replace subjects w/o peak freq with group avg center freq for each
        %band
        pk_theta(isnan(pk_theta)) = nanmean(pk_theta);
        pk_alpha(isnan(pk_alpha)) = nanmean(pk_alpha);
        pk_beta(isnan(pk_beta)) = nanmean(pk_beta);
        pk_gamma(isnan(pk_gamma)) = nanmean(pk_gamma);

        %consolidate subjects PSD  by band
        window = 2; %(Hz) size of window centered at center peak freq
        alpha = []; beta =[]; theta = []; gamma= []; condition = []; cond_num  =[];

        for i=1:size(clust_spec,2)
            data = clust_spec(i).spec_smooth(:,:,~bad_sub_ind);
            for seti = 1:size(data,3)
                B(:,seti) = mean(data(freqs >=(pk_beta(1,seti)-(window/2)) & freqs <= (pk_beta(1,seti)+(window/2)),:,seti),1)'; %grab spec data in window centered at pk freq
                T(:,seti) = mean(data(freqs >=(pk_theta(1,seti)-(window/2)) & freqs <= (pk_theta(1,seti)+(window/2)),:,seti),1)';
                A(:,seti) = mean(data(freqs >=(pk_alpha(1,seti)-(window/2)) & freqs <= (pk_alpha(1,seti)+(window/2)),:,seti),1)';
                G(:,seti) = mean(data(freqs >=(pk_gamma(1,seti)-(window/2)) & freqs <= (pk_gamma(1,seti)+(window/2)),:,seti),1)';
            end
            beta = [beta;B];
            theta = [theta;T];
            alpha = [alpha;A];
            gamma = [gamma;G];
            M = size(data,2);
            condition = [condition; repmat({myplotParams.labels{1,i}},[M,1])];
            cond_num = [cond_num; repmat(i,[M,1])];
            clear B T A G
            %         cluster = [cluster; repmat(opt.clusters(CL),[M,1])];
        end

        all_cond(c).data.subPSDbands.theta = theta;
        all_cond(c).data.subPSDbands.alpha = alpha;
        all_cond(c).data.subPSDbands.beta = beta;
        all_cond(c).data.subPSDbands.gamma = gamma;
        all_cond(c).data.subPSDbands.condition = condition;
        all_cond(c).data.subPSDbands.cond_num = cond_num;

        % % or use commonly defined frequency bands
        %     alpha = []; beta =[]; theta = []; gamma= []; condition = []; cond_num  =[];
        %     for i=1:size(clust_spec,2)
        %         data = clust_spec(i).spec;
        %         beta = [beta, data(freqs>=13 & freqs<30,:,:)];
        %         theta = [theta, data(freqs>=4 & freqs<7,:,:)];
        %         alpha = [alpha, data(freqs>=8 & freqs<13,:,:)];
        %         gamma = [gamma, data(freqs>=30 & freqs<128,:,:)];
        %         M = size(data,2);
        %         condition = [condition; repmat({myplotParams.labels{1,i}},[M,1])];
        %         cond_num = [cond_num; repmat(i,[M,1])];
        % %         cluster = [cluster; repmat(opt.clusters(CL),[M,1])];
        %     end

        %% avg PSD power in each freq band

        %theta band
        data_movmean = movmean(theta,10,1); %moving average across k strides
        myspecplot = figure('Position',[80 800 1000 400],'Units','pixels','InvertHardcopy', 'off', 'PaperType', 'a2', 'PaperOrientation', 'portrait');
        set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[1 2 6.5 2],'PaperPositionMode','auto');
        subplot(1,4,1);
        sgtitle({'Power Spectral Density at Peak Frequency',strcat('CL-', num2str(CL),'-',label)},'Interpreter','none')
        plotPSD_across_cond(data_movmean', clust_spec,cond_num,myplotParams)
        title('Theta')


        %alpha band
        data_movmean = movmean(alpha,10,1); %moving average across k strides
        subplot(1,4,2);
        plotPSD_across_cond(data_movmean', clust_spec,cond_num,myplotParams)
        title('Alpha')


        % beta band
        % plot beta moving mean, 10 strides
        data_movmean = movmean(beta,10,1); %moving average across k strides
        subplot(1,4,3);
        plotPSD_across_cond(data_movmean',clust_spec,cond_num,myplotParams)
        title('Beta')


        data_movmean = movmean(gamma,10,1); %moving average across k strides
        subplot(1,4,4);
        plotPSD_across_cond(data_movmean',clust_spec,cond_num,myplotParams)
        title('Gamma')

        % Give common xlabel, ylabel and title to your figure
        han=axes(myspecplot ,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        Y = ylabel(han,'Power 10*log_1_0 (uV^2/Hz)');
        Y.Position(1) = Y.Position(1)-0.02;
        X = xlabel(han,'Strides');
        X.Position(2) = X.Position(2)-0.01;
        set(gcf,'Color','w')
        set(X,'Fontsize',12)
        set(Y,'Fontsize',12)

        %% save figure
        savethisfig(myspecplot,['Band_10strideMovingAvgPSD_xCondition',myplotParams.figname,num2str(CL),'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg\'],'jpg')
        savethisfig(myspecplot,['Band_10strideMovingAvgPSD_xCondition',myplotParams.figname,num2str(CL),'.fig'],[savePath,'\PSD\',myplotParams.figname,'\fig\'],'fig')
        print([savePath,'\PSD\',myplotParams.figname,'\dpdf\Band_10strideMovingAvgPSD_xCondition',myplotParams.figname,num2str(CL),'.dpdf'], '-dpdf', '-painters'); % Makoto's print method. On Linux.

        %%
        subplot(1,4,1);
        plot(movmean(theta,10,1))
        ylim([min(movmean(theta,10,1),[],'all'),max(movmean(theta,10,1),[],'all')])

        subplot(1,4,2);
        plot(movmean(alpha,10,1))
        ylim([min(movmean(theta,10,1),[],'all'),max(movmean(alpha,10,1),[],'all')])

        subplot(1,4,3);
        plot(movmean(beta,10,1))
        ylim([min(movmean(beta,10,1),[],'all'),max(movmean(beta,10,1),[],'all')])

        subplot(1,4,4);
        plot(movmean(gamma,10,1))
        ylim([min(movmean(gamma,10,1),[],'all'),max(movmean(gamma,10,1),[],'all')])

        savethisfig(myspecplot,['Band_10strideMovingAvgPSD_xCondition_bySub',myplotParams.figname,num2str(CL),'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg\'],'jpg')
        savethisfig(myspecplot,['Band_10strideMovingAvgPSD_xCondition_bySub',myplotParams.figname,num2str(CL),'.fig'],[savePath,'\PSD\',myplotParams.figname,'\fig\'],'fig')
        print([savePath,'\PSD\',myplotParams.figname,'\dpdf\Band_10strideMovingAvgPSD_xCondition_bySub',myplotParams.figname,num2str(CL),'.dpdf'], '-dpdf', '-painters'); % Makoto's print method. On Linux
        close;

        %temp, uncommment later
            %% plot all condition psd by subject (subplots)
            myspecplot = figure('Position',[80 800 1000 1000],'Units','pixels');
            plotPSD_across_cond_bySub (clust_spec,freqs,CL,myplotParams)
            sgtitle(['Power Spectral Density-',num2str(CL),'-',label],'Interpreter','none')
            set(gcf,'PaperUnits','inches','Units','Inches','PaperPosition',[1 1 8 8]);
            han=axes(myspecplot ,'visible','off');
            han.Title.Visible='on';
            han.XLabel.Visible='on';
            han.YLabel.Visible='on';
            Y = ylabel(han,'Power 10*log_1_0 (uV^2/Hz)');
            X = xlabel(han,'Frequency (Hz)');
            set(X,'Fontsize',16)
            set(Y,'Fontsize',16)
        
            %%
            savethisfig(myspecplot,['PSD_',myplotParams.figname,num2str(CL),'_bySubject.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg\'],'jpg')
            savethisfig(myspecplot,['PSD_',myplotParams.figname,num2str(CL),'_bySubject.fig'],[savePath,'\PSD\',myplotParams.figname,'\fig\'],'fig')
            %     savethisfig(myspecplot,['Band_10strideMovingAvgPSD_xCondition',myplotParams.figname,num2str(CL),'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\svg\'],'svg')
            close;
    end
%         CL_PSD(c).clust_num = CL;
%         CL_PSD(c).clust_spec = clust_spec;
%         CL_PSD(c).freqs = freqs;

end



cd([savePath,'\PSD\',myplotParams.figname])
save('CL_PSD.mat','CL_PSD')

%make summary stats and export to csv
% calcualte summary stats
stats = psd_stats(all_cond, mylegend);
%export stats
writetable(struct2table(stats),[savePath,'\PSD\',myplotParams.figname,'\summaryStatistics.csv']);
% all_cond.summary_stats = stats;

%%
disp('Done making spectral plots')
fclose('all')
diary ON
fprintf('%s\nSTUDY Name:%s\nPath: %s',date,STUDY.filename,STUDY.filepath);
fprintf('Number of subjects= %i\n',length(unique({STUDY.datasetinfo.subject})));
fprintf('\nBaseline normalization using average frequency power across all frequencies and all conditions');
fprintf('\nSpectra Paramters:\n\t-Spectra smoothing width = %i\n\t-Stats method = %s\n\t-alpha (p-value thresh)= %.2f\n\t-correction for multiple comparisons= %s\n',width,statsMethod, Alpha, mcorrect);
fprintf('\nBand power calculated using subject peak freq. within each predefined band. Subjects w/o peaks used the group mean as their center freq for each band. Outliers removed. Shaded bars showing standard error of the mean')
diary OFF

end




function plotPSD_across_cond(data,clust_spec,cond_num,myplotParams)
% plot condition filler colors first
y = min(data,[],'all'):0.1:max(data,[],'all');
y2 = [y, fliplr(y)];
for i = 1:size(clust_spec,2)
    mycond = find(cond_num == i);
    x = repmat(mycond(1),[size(y)]);
    x2 = repmat(mycond(end),[size(y)]);
    inBetween = [x, fliplr(x2)];
    fill(inBetween,y2,myplotParams.colors{1,i},'FaceAlpha',0.2,'EdgeColor','none');
    hold on;
end
xlim([0 90])

%     for i = 1:size(data,2)
%         plot(data(:,i));
%     end
% data = squeeze(data)';
% eBar = std(data_movmean,0,2); %std across subjects
eBar = nanstd(data,0,1)/sqrt(size(data,1)); %std across subjects  , SEM = std(data)/sqrt(length(data));
shadedErrorBar([1:size(data,2)], nanmean(data,1) ,eBar, 'lineprops',{'-','Color',[0 0 0],'LineWidth',2},'transparent',true)

%shadedErrorBar([1:size(data,2)], data ,{@mean @std}, 'lineprops',{'-','Color',[0 0 0],'LineWidth',2},'transparent',true)
ylim([min(nanmean(data,1),[],'all')-2,max(nanmean(data,1),[],'all')+2]) %#ok<NOPRT>
xticks([0 100 170])
set(gca,'Fontsize',10,'FontName','Arial');
end

function plotPSD_across_cond_bySub (clust_spec,freqs,CL,myplotParams)
global STUDY

num_subs = size(clust_spec(1).spec,3);
for sub = 1:num_subs
    subplot(5,6,STUDY.cluster(CL).sets(sub))
    for cond = 1:size(clust_spec,2)
        data = mean(clust_spec(cond).spec_smooth(:,:,sub),2);
        plot(freqs,data,'Color',myplotParams.colors{1,cond});
        hold on;
    end


    if sub ==1
        %         legend(mylegend)
        myRange(1) = min(data);
        myRange(2) = max(data);
    else
        currMin = min(data);
        currMax = max(data);

        if currMin < myRange(1)
            myRange(1) = currMin;
        end
        if currMax > myRange(2)
            myRange(2) = currMax;
        end
    end
    xlim([4 30])
    ylim([myRange(1) myRange(2)+5])

    set(gca,'FontName','Arial');
    title(['Subject ',num2str(STUDY.cluster(CL).sets(sub))]);
    set(gca,'Fontsize',10,'FontName','Arial');
end
set(gcf,'Color','w')

end

%plot psd comparision between two conditions
function plot_psd_comparison(STUDY, CL, specfreqs,specdata_avg,statsMethod,cluster_perm_test_i,conditionIndex, myplotParams, mylegend,savePath)
freqi = find(specfreqs>=cluster_perm_test_i.freqrange(1) & specfreqs<=cluster_perm_test_i.freqrange(2));
try
    label = cellstr(STUDY.cluster(CL).label);
    label = label{1,1};
    mytitle= [strcat(myplotParams.title,num2str(CL),'-',label)];
catch
    mytitle= [strcat(myplotParams.title,num2str(CL))];
end
%colors = [{lightOrange},{red},{lightRed},{purple},{darkPurple}, {lightPurple},];
colors = {myplotParams.colors{1,conditionIndex}};

plotopt = {'highlightmode','bottom','plotmean','off','ylim',[], 'xlabel','Frequency (Hz)','ylabel',...
    'Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)','legend',{mylegend{1,conditionIndex}}};
myspecplot = figure('InvertHardcopy', 'off', 'PaperType', 'a2', 'PaperOrientation', 'landscape');
plotcurve_colors( specfreqs(freqi),specdata_avg(conditionIndex,freqi), 'colors', colors, 'maskarray', cluster_perm_test_i.pcond', plotopt{1:end}, 'title', mytitle);

%set figure settings
set(gcf,'Color','w')
set(gcf,'PaperUnits','Inches','Units','Inches','PaperPosition',[1 0 8 4.2],'Position', [1 0 8 4.2]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
pos = get(gca,'Position');
xticklabel = get(gca,'XTickLabel');
yticklabel = get(gca,'YTickLabel');
ytick = get(gca,'yTick');
%set(gca,'Visible','off')
% axes('Position',pos,'XAxisLocation','bottom','YAxisLocation','left',...
%  'Color','none','XTickLabel',xticklabel,'YTickLabel',yticklabel,...
%  'XColor','k','YColor','k')
%axes('Position',pos,'XAxisLocation','bottom','YAxisLocation','left')

ylh = get(gca,'ylabel');ylp = get(ylh, 'Position');
%         set(ylh,'FontSize',12)
set(ylh, 'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',12);
% set(ylh,'Position',[ylp(1)-3 ylp(2) ylp(3)])
xlh = get(gca,'xlabel'); set(xlh,'FontSize',12);


legend({mylegend{1,conditionIndex}})

lgd = findobj(gcf, 'Type', 'Legend'); lgd.FontSize = 7;
%         set(lgd,'Position',[pos(1)-150 pos(2)+90 pos(3) pos(4)]);
%         set(lgd, 'EdgeColor','none','Orientation','horizontal')
set(lgd, 'EdgeColor','none','Orientation','horizontal')
pos = get(lgd, 'Position');
set(lgd,'Position',[0 0.95 pos(3) pos(4)]);
box off


Alpha = STUDY.etc.statistics.fieldtrip.alpha;
mcorrect = STUDY.etc.statistics.fieldtrip.mcorrect;
mode = STUDY.etc.statistics.mode;
%save figure
savethisfig(myspecplot,['PSD_',cluster_perm_test_i.name,'comparison_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(cluster_perm_test_i.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.jpg'],[savePath,'\PSD\',myplotParams.figname,'\jpg'],'jpg')
savethisfig(myspecplot,['PSD_',cluster_perm_test_i.name,'comparison_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(cluster_perm_test_i.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.fig'],[savePath,'\PSD\',myplotParams.figname,'\fig'],'fig')
savethisfig(myspecplot,['PSD_',cluster_perm_test_i.name,'comparison_',myplotParams.figname,'_CL',num2str(CL),'_',num2str(round(max(cluster_perm_test_i.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.svg'],[savePath,'\PSD\',myplotParams.figname,'\svg'],'svg')
if ~exist([savePath,'\PSD\',myplotParams.figname,'\dpdf\'], 'dir') %check
    mkdir([savePath,'\PSD\',myplotParams.figname,'\dpdf\'])
end
print([savePath,'\PSD\',myplotParams.figname,'\dpdf\PSD_',cluster_perm_test_i.name,'comparison_',myplotParams.figname,num2str(CL),'_',num2str(round(max(cluster_perm_test_i.freqrange))),'Hz_',statsMethod,'_alpha',num2str(Alpha),'_',mcorrect,'_',mode,'.dpdf'], '-dpdf', '-painters'); % Makoto's print method. On Linux.
close;
end


function stats = psd_stats(allcond, mylegend)
stats = struct([]);
cond = mylegend;
for CLi=1:size(allcond,2)
    allpsd = allcond(CLi).data.allpsd;
    if size(cond,2)~= size(allpsd,1)
        error('Number of conditions doesn''t match size of psd data');
    end
    freqs = allcond(CLi).data.freqs;
    for condi = size(allpsd,1)
        data = allpsd{condi,1}; %[freq x sub]
        theta = mean(data(freqs>=4 & freqs<7,:),1); %avg data in freq band, [1 x subs]
        alpha = mean(data(freqs>=8 & freqs<13,:),1);
        beta = mean(data(freqs>=13 & freqs<30,:),1);
        gamma = mean(data(freqs>=30 & freqs<128,:),1);

        freqBand = {'theta','alpha','beta','gamma'};
        for f = 1:length(freqBand)
            tmpdata = eval(freqBand{1,f});
            stats(end+1).CL = allcond(CLi).cluster;
            stats(end).CL_label = allcond(CLi).label;
            stats(end).Condition = mylegend{1,condi};
            stats(end).Freq = freqBand{1,f};
            stats(end).Mean = mean(tmpdata);
            stats(end).SD = std(tmpdata);
        end
    end
end
end

%determine effect size for cluster with lowest p-val
%modified from Arnald Delorme: https://github.com/Donders-Institute/infant-cluster-effectsize/blob/main/do_group_analysis.m
function [effect] = calc_clust_effectsize(specdata_all_tmp,freq,pval,method)
%determine cluster with lowest p-val
p= unique(pval);
p = p(p<0.05);

for clusti = 1:length(p)
    effectWindow = find(pval==p(clusti));
    %calculate pairwise difference in psd between conditions
    %for each participant
    cond1psd = specdata_all_tmp{1,1}(effectWindow,:); %[freq x sub]
    cond2psd = specdata_all_tmp{2,1}(effectWindow,:);
    c = cond2psd-cond1psd;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Option 1: Calculate Cohen's d for the average difference
    % in the respective cluster
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if method == 1
        psd_diff = nanmean(c,1); %avg across freq band
        %calculate Cohen's d
        effect(clusti).method = 'avg difference in cluster';
        effect(clusti).SD = std(psd_diff);
        effect(clusti).MEAN = mean(psd_diff);
        effect(clusti).COHENS_D = mean(psd_diff)/std(psd_diff);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Option 2: Determine at maximum effect size and at which channel/time it
    % is maximal (upper bound)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if method ==2
        % Determine maximum effect sizavg across subjectse and at which frequency Cohen's d is maximal
        cd = abs(nanmean(c,2)./std(c,0,2));%%abs so doesn't matter if clusters is positive or negative, avg across subjects
        maxcd= max(cd); %
        maxeffectfreq = freq(cd == maxdiff);
        effect(clusti).method = 'max effect size';
        effect(clusti).COHENS_D = maxcd;
        effect(clusti).maxeffectfreq=  freq(effectWindow(cd == maxcd));
    end


    if method ==0
        psd_diff = nanmean(c,1); %avg across freq band
        cd = abs(nanmean(c,2)./std(c,0,2));%%abs so doesn't matter if clusters is positive or negative, avg across subjects
        maxcd= max(cd); %
        maxeffectfreq = freq(effectWindow(cd == maxcd));
        rowi = find(cd==maxcd);

        %calc 95% CI just for freq with max cd
        e = meanEffectSize(abs(c(rowi,:)),Effect="cohen",ConfidenceIntervalType="bootstrap", ...
            BootstrapOptions=statset(UseParallel=true,type='norm'),NumBootstraps=3000); %idk why this function isn't working after matlab was reinstalled



        effect(clusti).SD = std(psd_diff);
        effect(clusti).MEAN = mean(psd_diff);
        effect(clusti).COHENS_D_avg = mean(psd_diff)/std(psd_diff);
        effect(clusti).COHENS_D_max = e.Effect;
        effect(clusti).CI95 = [e.ConfidenceIntervals];
        effect(clusti).maxeffectfreq =  [maxeffectfreq];
        effect(clusti).window = [freq(effectWindow(1)), freq(effectWindow(end))];

    end
end
end

