function STUDY = diplotfig(STUDY,ALLEEG,clusters_to_plot, colors2use,filesavepath,centrDipsTogether,one_sub_per_cl)

%DIPLOTFIG() plots centroids and ICs for selected clusters
% assumes STUDY is already loaded and ALLEEG is populated
% clusters_to_plot = array of the cluster numbers you want to plot
% colors2use = cell with [R G B] values you want to use to for each cluster
% colors2use{1} = [1 0 0];
% colors2use{2} = [0 1 0];
% colors2use{3} = [0 0 1];
% etc.
% filesavepath = where the figures will be saved

%modification 7/28/20, NJacobsen
%changed export_fig function to saveas for all png files
%modification 11/4/21
%cropped images, option to plot one dipole per subject

if ~isdir(filesavepath)
    mkdir(filesavepath)
end

cd(filesavepath);
if centrDipsTogether==1
    centroidPlotSize=50; %60; %
else
    centroidPlotSize=0;
end

k1=1;
k2=1;
compsall=[]; % structure with all of the component posxyz, momxyz, rv
centsall=[]; % structure with all of the cluster centroid posxyz, momxyz, rv
dipsizes = [];

ct = 1;
for i=1:length(clusters_to_plot)
    k2=k2+1;
    if one_sub_per_cl ==0
        for j=1:length(STUDY.cluster(clusters_to_plot(i)).comps)

            cls_set_i = STUDY.datasetinfo(STUDY.cluster(clusters_to_plot(i)).sets(1,j)).index;
            if ~isfield(ALLEEG(cls_set_i), 'dipfit')
                warndlg2(['No dipole information available in dataset ' ALLEEG(cls_set_i).filename ' , abort plotting'], 'Aborting plot dipoles');
                return;
            end
            comp = STUDY.cluster(clusters_to_plot(i)).comps(j);
            cluster_dip_models(1,j).posxyz = ALLEEG(cls_set_i).dipfit.model(comp).posxyz;
            cluster_dip_models(1,j).momxyz = ALLEEG(cls_set_i).dipfit.model(comp).momxyz;
            cluster_dip_models(1,j).rv = ALLEEG(cls_set_i).dipfit.model(comp).rv;

            colorsall{ct} = colors2use{i};
            colorsc{k1} = colors2use{i};
            colorsc_1{k1} = [0 0 1];
            %         coltemp = colors2use{i}/2; %+[.5 .5 .5];
            %         colors1cls{j} = coltemp; %/(max(coltemp));
            %         colorsall{ct} = coltemp; %/(max(coltemp));
            colors1cls{j} = colors2use{i};
            ct= ct+1; k1=k1+1;

        end
    else
        unique_clus_subs = unique(STUDY.cluster(clusters_to_plot(i)).sets);
        for j = 1:length(unique_clus_subs)
            cls_set_i = find(STUDY.cluster(clusters_to_plot(i)).sets == unique_clus_subs(j));
            comp = STUDY.cluster(clusters_to_plot(i)).comps(cls_set_i);
           
                if size(cls_set_i,2)>1 %if subject appears more than once in cluster

                    if one_sub_per_cl ==1 %average values
                    for c = 1:length(comp)
                        posxyz(c,1) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).posxyz(1);
                        posxyz(c,2) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).posxyz(2);
                        posxyz(c,3) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).posxyz(3);
                        momxyz(c,1) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).momxyz(1);
                        momxyz(c,2) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).momxyz(2);
                        momxyz(c,3) = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).momxyz(3);
                        rv(c)= ALLEEG(unique_clus_subs(j)).dipfit.model(comp(c)).rv;
                    end
                    cluster_dip_models(1,j).posxyz = mean(posxyz,1);
                    cluster_dip_models(1,j).momxyz = mean(momxyz,1);
                    cluster_dip_models(1,j).rv = mean(rv);
                    elseif one_sub_per_cl ==2 %select lowest IC #
                        cluster_dip_models(1,j).posxyz  = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(1)).posxyz;
                        cluster_dip_models(1,j).momxyz = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(1)).momxyz;
                        cluster_dip_models(1,j).rv = ALLEEG(unique_clus_subs(j)).dipfit.model(comp(1)).rv;
                    end
                else
                    cls_set_i = STUDY.cluster(clusters_to_plot(i)).sets(cls_set_i);
                    cluster_dip_models(1,j).posxyz = ALLEEG(cls_set_i).dipfit.model(comp).posxyz;
                    cluster_dip_models(1,j).momxyz = ALLEEG(cls_set_i).dipfit.model(comp).momxyz;
                    cluster_dip_models(1,j).rv = ALLEEG(cls_set_i).dipfit.model(comp).rv;
                end
                colorsall{ct} = colors2use{i};
                colorsc{k1} = colors2use{i};
                colorsc_1{k1} = [0 0 1];
                %         coltemp = colors2use{i}/2; %+[.5 .5 .5];
                %         colors1cls{j} = coltemp; %/(max(coltemp));
                %         colorsall{ct} = coltemp; %/(max(coltemp));
                colors1cls{j} = colors2use{i};
                ct= ct+1; k1=k1+1;
                clear posxyz momxyz rv
                end
    end


    compsall=[compsall,cluster_dip_models(1,:)];
    centsall = [centsall, computecentroid(cluster_dip_models)];
    dipsizes = [dipsizes 25*ones(size(1,length(cluster_dip_models(1,:)))) 40];
    %     coltemp = colors2use{i}/2+[.5 .5 .5];
    %     colors1cls{end+1} = coltemp/(max(coltemp));
    colors1cls{end+1} = colors2use{i}/2;

    %     plots individual cluster and its components
    dipplot([cluster_dip_models(1,:), computecentroid(cluster_dip_models)],'spheres','on','dipolelength',0,...
        'dipolesize',[20*ones(1,length(cluster_dip_models(1,:))) centroidPlotSize],...
        'mri',ALLEEG(1).dipfit.mrifile,'meshdata',ALLEEG(1).dipfit.hdmfile,...
        'coordformat',ALLEEG(1).dipfit.coordformat,'color',colors1cls);
    %     dipplot([compsall centsall],'spheres','on','dipolelength',0,'dipolesize',[20*ones(size(compsall)) centroidPlotSize*ones(size(centsall))],...
    %     'mri','C:\Program Files\MATLAB\R2018a\toolbox\eeglab14_1_2b\plugins\dipfit3.3\standard_BEM\standard_mri.mat','color',colors1cls); %
    %     set(findobj('tag','img'), 'facealpha', 0.6);
    %     set(findobj('facealpha',1),'facelighting','phong');

    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents.fig'])
    view([1,0,0])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_side.jpg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_side.svg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_side.png'])
    I = imread([filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_side.jpg']);
    D1  = imcrop(I,[193.5 62.5 612 510]); %crop image
    imwrite(D1,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_side.jpg']);

    view([0,-1,0])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_back.jpg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_back.svg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_back.png'])
    I = imread([filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_back.jpg']);
    D2  = imcrop(I,[244.5 79.5 509 492]); %crop image
    imwrite(D2,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_back.jpg']);
    view([0,0,1])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_top.jpg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_top.svg'])
    saveas(gcf,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_top.png'])
    I = imread([filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_top.jpg']);
    D3  = imcrop(I,[243.5 11.5 509 610]); %crop image
    imwrite(D3,[filesavepath filesep 'CLS_'  num2str(clusters_to_plot(i)) '_comps_cents_top.jpg']);

    STUDY.cluster(clusters_to_plot(i)).centroid = computecentroid(cluster_dip_models);
    clear cluster_dip_models colors1cls
end

STUDY.etc.centroids =  centsall;
%%
% colors2use
for i=1:length(clusters_to_plot)
    colorsall{ct-1+i} = colors2use{i}/2;
end

% plots all clusters and their components, color coded
dipplot([compsall centsall],'spheres','on','dipolelength',0,'dipolesize',[20*ones(size(compsall)) centroidPlotSize*ones(size(centsall))],...
    'mri',ALLEEG(1).dipfit.mrifile,'meshdata',ALLEEG(1).dipfit.hdmfile,...
    'coordformat',ALLEEG(1).dipfit.coordformat,'color',colorsall); %

set(findobj('tag','img'), 'facealpha', 0.6);
set(findobj('facealpha',1),'facelighting','phong');
set(gcf,'Color','k')
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps.fig'])
view([1,0,0]) %view([40 50]); %
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_side.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_side.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_side.png'])
view([0,-1,0])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_back.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_back.png'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_back.svg'])
view([0,0,1])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_top.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_top.png'])
saveas(gcf,[filesavepath filesep STUDY.name '_allclust_comps_top.svg'])

% plots all components, not color coded, all blue, no centroids
dipplot(compsall,'spheres','on','dipolelength',0,'dipolesize',20,...
    'mri',ALLEEG(1).dipfit.mrifile,'meshdata',ALLEEG(1).dipfit.hdmfile,...
    'coordformat',ALLEEG(1).dipfit.coordformat,'color',colorsc_1);
% dipplot(compsall,'spheres','on','dipolelength',0,'dipolesize',20,...
%     'mri','C:\Program Files\MATLAB\R2018a\toolbox\eeglab14_1_2b\plugins\dipfit3.3\standard_BEM\standard_mri.mat');
set(findobj('tag','img'), 'facealpha', 0.6);
set(findobj('facealpha',1),'facelighting','phong');
set(gcf,'Color','k')
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped.fig'])
view([1,0,0])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_side.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_side.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_side.png'])
view([0,-1,0])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_back.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_back.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_back.png'])
view([0,0,1])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_top.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_top.png'])
saveas(gcf,[filesavepath filesep STUDY.name '_allcomps_ungrouped_top.svg'])

% plots only cluster centroids, no components
dipplot(centsall,'spheres','on','dipolelength',0,'dipolesize',50,... %1,... %
    'mri',ALLEEG(1).dipfit.mrifile,'meshdata',ALLEEG(1).dipfit.hdmfile,...
    'coordformat',ALLEEG(1).dipfit.coordformat,'color',colors2use); %,'projlines','on');
set(findobj('tag','img'), 'facealpha', 0.6);
set(findobj('facealpha',1),'facelighting','phong');
set(gcf,'Color','k')
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only.fig'])
view([1,0,0]) %view([40 50]); %
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_side.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_side.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_side.png'])
view([0,-1,0])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_back.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_back.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_back.png'])
view([0,0,1])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_top.jpg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_top.svg'])
saveas(gcf,[filesavepath filesep STUDY.name '_centroids_only_top.png'])

%%
function dipole = computecentroid(alldipoles)

max_r = 0;
len = length(alldipoles);
dipole.posxyz = [ 0 0 0 ];
dipole.momxyz = [ 0 0 0 ];
dipole.rv = 0;
ndip = 0;
count = 0;
warningon = 1;
for k = 1:len
    if size(alldipoles(k).posxyz,1) == 2
        if all(alldipoles(k).posxyz(2,:) == [ 0 0 0 ])
            alldipoles(k).posxyz(2,:) = [];
            alldipoles(k).momxyz(2,:) = [];
        end;
    end;
    if ~isempty(alldipoles(k).posxyz)
        dipole.posxyz = dipole.posxyz + mean(alldipoles(k).posxyz,1);
        dipole.momxyz = dipole.momxyz + mean(alldipoles(k).momxyz,1);
        dipole.rv     = dipole.rv     + alldipoles(k).rv;
        count = count+1;
    elseif warningon
        disp('Some components do not have dipole information');
        warningon = 0;
    end;
end
dipole.posxyz = dipole.posxyz/count;
dipole.momxyz = dipole.momxyz/count;
dipole.rv     = dipole.rv/count;
if isfield(alldipoles, 'maxr')
    dipole.maxr = alldipoles(1).max_r;
end;