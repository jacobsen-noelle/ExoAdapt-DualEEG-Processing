%Custom Head Model from MRI
% Creates a custom head model using subject specific MRI
% Expected folder hierarchy: ...\Data\SXX\MRI                           
% Required inputs:
%   savePath             - output folder where new EEG dataset is saved after
%                          dipole fitting
%   eeglab_path          - path to eeglab
%   electrodelocs_path   - input path to electrode locations electrodes file
%                          format in folder should be as follows: SXX_chanlocs.txt
%   subjectList          - array of subject numbers to process [] 
%                          currently does not handle string list, but could be updated
%   EEGchanfile          - (automatically searching for SXX_chanlocs.txt) name of channel location .txt file, same as output
%                          text file created from getchanlocs().
%                          Must contain fiducial locations if you want to
%                          coregister headmodel/mri and electrode locs based 
%                          on fiduaicls {'nas','lhj','rhj'};
%   MRIfolder            - input folder with MRI file (.nii extension)
%   MRIfilename          - file with raw or MNI normalized MRI (.nii)
%   inputMRItype         - options:'raw','norm'; indicates whether the MRI file
%                          is raw or has already been normalized using SPM or
%                          other software
%   outputCoordSys       - options:'MNI','CTF' (recommended); output coordinate system of
%                          headmodel. 'MNI' means spatially normalized to MNI coordinate
%                          system, with origin at anterior commissure.
%                          'CTF' means subject coordinate space. If you use
%                          CTF, you should normalized your dipfit results
%                          to MNI space later for group analyses.
% Outputs:
% mri_segmented
% mri_seg_i
% fidLocsMRI             - location of fiducials (nas,lhj,rhj) marked on
%                          MRI. Units are in 'mm', to ensure compatibility
%                          with EEG electrodes
% headmodel_fem          - head model (.mat) created with finite element
%                          analysis (FEM)
% mesh_fem               - mesh (.mat)
%
% Optional inputs:
%  plotStuff             - turn plotting of figures ON/OFF
%
%
% Author:
% Noelle Jacobsen, University of Florida
% Created 8/10/2021
% Last updated 7/20/23
% TO DO: Add MRI nomalization at the end of script so dipole moments can be recomputed on HPG. Dont worry about this unless you are using dipole moments (e.g. if you use them as a feature to cluster)
% =========================================================================
% =========================================================================
%parameters
subjectList = [10 15 17 20:31]; % [already did S01,exclude exo adapt S12,S18, S19]
inputMRItype = 'raw';   
plotStuff =1; %turn plotting of figures ON/OFF
outputCoordSys = 'CTF'; % ['MNI'|'CTF'] mni coordinate space or local ** feature not actively used rn
%set paths
mydir = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation';
% fieldtrip_path = 'C:\Users\jacobsen.noelle\Desktop\fieldtrip-20210614'; %Fieldtrip folder
fieldtrip_path = 'R:\Ferris-Lab\share\MindInMotion\fieldtrip-master'; %path to full fieldtrip
electrodelocs_path =[mydir,'\Data\raw_data\electrode_locations']; % add your path to electrode locations .txt files from getchanlocs
eeglab_path ='C:\Users\jacobsen.noelle\Desktop\eeglab2022.0'; %CHANGE TO YOUR PATH to EEGlab-- only purpose is to remove Fieldtrip lite plugin
MRI_dir = [mydir,'\Data\raw_data\']; 
seg_MRI_dir = ''; %leave empty if you don't already have segmented MRI
savePath = [mydir,'\Data\raw_data']; %folder where subject head model files will be saved
codeRepo = 'R:\Ferris-Lab\jacobsen.noelle\Matlab_Scripts\ExoAdapt-EEG-Processing\';
addpath(genpath(codeRepo))

% or load subject list this way
folderList = dir('\\exasmb.rc.ufl.edu\blue\dferris\jacobsen.noelle\CustomHeadModel\Data\2023-09-22-Exo');
folderList = folderList(3:end);
subjectList = [];
for folderi = 1:size(folderList ,1)
    subjectList(folderi) = str2num(extractAfter(folderList(folderi).name,'S'));
end
subjectList = unique(subjectList);
%% load directory contents
folderList = dir(MRI_dir);

for subi = 1:length(subjectList)
    subject = subjectList(subi);
    if subject <10
        subject = ['0',num2str(subject)];
    else
        subject = num2str(subject);
    end
    subject = ['S',subject];
    subject = string(subject);
    fprintf('\n\nStarting %s\n',subject);

    folderi = find(contains({folderList.name},subject) & [folderList.isdir] ==1); %find subject folder in list,there should only be one
    if length(folderi) >1; error('too many folders found for this subject');end

    
    %folder
    MRIfolder = strcat(MRI_dir,'\',folderList(folderi).name,'\MRI');
    switch inputMRItype
        case {'raw'}
            MRIfilename = strcat('ExoAdapt_',subject,'_T1W.nii');
        case {'norm'}
            MRIfilename = strcat('wSBP_',subject,'_MRI_*.nii');
    end
    
    outputFolder = strcat(savePath,'\',folderList(folderi).name,'\Head Model');
    
    %setup fieldtrip
    addpath(fieldtrip_path)
    rmpath(genpath([eeglab_path,'\plugins\Fieldtrip-lite20230716'])) % CHANGE folder; having more than one path to fieldtrip will confuse matlab, update to your eeglab fieldtrip lite folder
    ft_defaults;
    %% find mri file
    cd (MRIfolder)
    mrifile = dir(fullfile(MRIfolder, MRIfilename)); %find MRI file in directory
    if isempty(mrifile)|| size(mrifile,1)>2
        
        error('Error selecting MRI file')
    else
        fprintf('Loading MRI file: %s',mrifile.name);
    end
    
    if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
        mkdir(outputFolder);
    end
    %% 1. Load MRI
    %Read MRI
    mri = ft_read_mri(mrifile.name);

%     %visualize mri (optional)
%    cfg=[];
%    ft_sourceplot(cfg,mri);

%     mri = ft_determine_coordsys(mri); %not necessary as we'll change
%     mrirs = ft_volumereslice(cfg,mri);
%     ft_sourceplot(cfg,mrirs);

%% 2. Realign to acpc coordinate system
% #####################################################################
% This code aims to identify fiducials ACPC and CTF on T1 image
% Originally written by Vanessa Cruz  
% Editted by Chang Liu - 2022-02-17
% Modified by Noelle Jacobsen - 2023-02-14

% check 
check_acpc = isfile(fullfile(outputFolder,'acpc_fiducials.mat')); 
check_ctf = isfile(fullfile(outputFolder,'ctf_fiducials.mat'));

if check_acpc == 1
    prompt = 'ACPC fiducials already exists, do you want to continue? Y/N : ';
    str = input(prompt,'s');
    if isempty(str) || str == 'N'
        return
    end
end
if check_ctf == 1
    prompt = 'CTF fiducials already exists, do you want to continue? Y/N : '; %default N
    str = input(prompt,'s');
    if isempty(str) || str == 'N'
        return
    end
end  
fprintf('\nMark anterior and posterior commissure\n')
%define anterior commissure and posterior commmissure (ACPC)
cfg              = [];
cfg.coordsys     = 'acpc'; 
cfg.method       = 'interactive';
cfg.viewmode    = 'ortho';
[mri_acpc] = ft_volumerealign(cfg, mri); % have figure selected when pressing keys, not typing in command line
acpc_fiducials   = mri_acpc.cfg.fiducial;

save(fullfile(outputFolder,strcat(subject,'_mri_acpc.mat')),'mri_acpc');
save(fullfile(outputFolder,strcat(subject,'_acpc_fiducials.mat')),'acpc_fiducials');

%realign the T1 raw image to acpc coordinate
% Save mri_acpc_rs.nii
cfg             = [];
cfg.filename    = char(fullfile(outputFolder,strcat(subject,'_MRI_acpc')));
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
cfg.datatype    = 'double';
ft_volumewrite(cfg, mri_acpc);

%% reslice data and make the scan isotropic
% While ft_volumerealign does not change the anatomical data, 
% but it adjusts the transformation matrix of the data, 
% ft_volumereslice will change the anatomical data, 
% i.e. it will arrange data in field anatomy according to the coordinate system
cfg             = [];
%         cfg.dim         = [300 300 300]; % default is [256 256 256]. increase the dim so that the head won't get cut-off at the boundary
cfg.resolution  = [1];
mri_acpc_rs     = ft_volumereslice(cfg,mri_acpc); %can rearrange anatomical data
mri_acpc_rs.coordsys    = 'acpc';

% ---- plot 
cfg     = [];
ft_sourceplot(cfg,mri_acpc_rs);

% Save mri_acpc_rs.nii
cfg             = [];
cfg.filename    = char(fullfile(outputFolder,strcat(subject,'_MRI_acpc_rs')));
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
cfg.datatype    = 'double';
ft_volumewrite(cfg, mri_acpc_rs);

% Save data as mat
save(fullfile(outputFolder,strcat(subject,'_mri_acpc_rs.mat')),'mri_acpc_rs')

% plot and save acpc fiducial points
% TO DO: Add figure to put crosshairs the fiducials just marked
% plotting fiducials 
ac = mri_acpc.cfg.fiducial.ac;
pc = mri_acpc.cfg.fiducial.pc;
xz = mri_acpc.cfg.fiducial.xzpoint;

vox2head = mri_acpc.transform; %same transform as mr_realign.transformprig
ac_mri = ft_warp_apply(vox2head, ac);
pc_mri = ft_warp_apply(vox2head, pc);
xz_mri = ft_warp_apply(vox2head, xz);

cfg = [];
cfg.location = ac_mri;
ft_sourceplot(cfg, mri_acpc); %
title('Anterior Commissure location')
saveas(gcf,fullfile(outputFolder,strcat(subject,'_fiducial_ac.jpeg')))

cfg.location = pc_mri;
ft_sourceplot(cfg, mri_acpc); %
title('Posterior Commissure location')
saveas(gcf,fullfile(outputFolder,strcat(subject,'_fiducial_pc.jpeg')))

cfg.location = xz_mri;
ft_sourceplot(cfg, mri_acpc); %
title('XZ Point')
saveas(gcf,fullfile(outputFolder,strcat(subject,'_fiducial_xz.jpeg')))

close all;

%normalize MRI to MNI space-- will be used later after dipole fitting
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);

% template location
%cfg             = [];
%cfg.filename    = fullfile(strcat(outputFolder,'\',subject,'_mri_norm'));
%cfg.filetype    = 'nifti';
%cfg.parameter   = 'anatomy';
%ft_volumewrite(cfg, mri_norm);
cd(outputFolder)
save(strcat(subject,'_mri_norm.mat'),'mri_norm')
%#####################################################################

%% 3. Locate fiducials in MRI
    %use ft_volumerealign to mark fiducials (nas, lhj, rhj), but don't actually
    %use the output MRI
    %use same fiducials marked in EEG system (using 2lhj/rhj, instead of
    %lpa/rpa)
fprintf('\nMark fiducials\n')
cfg             = [];
cfg.coordsys    = 'ctf';
cfg.method      = 'interactive';
cfg.viewmode    = 'ortho';
[mri_ctf] = ft_volumerealign(cfg, mri_acpc_rs); 
ctf_fiducials = mri_ctf.cfg.fiducial;

   %grab fiducial locations you just marked using ft_volume realign
    nas = mri_ctf.cfg.fiducial.nas;
    lhj = mri_ctf.cfg.fiducial.lpa;
    rhj = mri_ctf.cfg.fiducial.rpa;

    vox2head = mri_ctf.transform; %same transform as mr_realign.transformprig
    lhj = ft_warp_apply(vox2head, lhj,'homogenous');
    rhj = ft_warp_apply(vox2head, rhj,'homogenous');
    nas = ft_warp_apply(vox2head, nas,'homogenous');


    %verify fiducials are marked properly (look at cross hair location
    cfg = [];
    cfg.location = nas;
    nas_fig = figure;
    ft_sourceplot(cfg, mri_ctf); %
    title('Nasion location'); fig = gcf; fig.OuterPosition = [2 300 500 500];
    cfg.location = lhj;
    lhj_fig = figure;
    ft_sourceplot(cfg,mri_ctf); %
    title('LHJ location');fig = gcf; fig.OuterPosition = [502 300 500 500];
    cfg.location = rhj;
    rhj_fig = figure;
    ft_sourceplot(cfg, mri_ctf); %
    title('RHJ location'); fig = gcf; fig.OuterPosition = [1002 300 500 500];
    
    disp('Check fiducial locations')
    pause; 
    
    prompt = 'No, remark';
    while ~strcmp(prompt,'Yes')
        %asks if fiducial locations are okay
        prompt = questdlg('Are the fiducials in the correct location? ', 'Verify fiducial location','Yes', 'No, remark','Plot in template MRI to verify','Yes');
        %response
        switch prompt
            case 'Yes'
                disp([' Everything looks good!']);
            case 'No, remark'
                disp(['Remark fiducial locations']);
                cfg = [];
                cfg.method = 'interactive';
                cfg.coordsys = 'ctf'; %spm doesn't work bc won't let me mark nas, need to mark in this coodsys
                mri_ctf= ft_volumerealign(cfg, mri_acpc); %you won't actually use this mri, just getting fiducial locations
                ctf_fiducials = mri_ctf.cfg.fiducial;
                %grab fiducial locations you just marked using ft_volume realign
                nas = mri_ctf.cfg.fiducial.nas;
                lhj = mri_ctf.cfg.fiducial.lpa;
                rhj = mri_ctf.cfg.fiducial.rpa;
                %transform to coordinate system identified in above
                %case/switch as vox2head variable

                vox2head = mri_ctf.transform; %same transform as mr_realign.transformprig
                nas = ft_warp_apply(vox2head, nas, 'homogenous');
                lhj = ft_warp_apply(vox2head, lhj, 'homogenous');
                rhj = ft_warp_apply(vox2head, rhj, 'homogenous');
                %verify fiducials are marked properly (look at cross hair location
                cfg = []; cfg.location = nas;
                nas_fig = figure;
                ft_sourceplot(cfg, mri_ctf); title('Nasion location');fig = gcf; fig.OuterPosition = [2 300 500 500];
                cfg.location = lhj;
                lhj_fig = figure;
                ft_sourceplot(cfg, mri_ctf); title('LHJ location'); fig = gcf; fig.OuterPosition = [502 300 500 500];
                cfg.location = rhj;
                rhj_fig = figure;
                ft_sourceplot(cfg, mri_ctf); title('RHJ location'); fig = gcf; fig.OuterPosition = [1002 300 500 500];
             end
    end

   saveas(nas_fig,fullfile(outputFolder,strcat(subject,'_fiducial_nas.jpeg')))
   saveas(rhj_fig,fullfile(outputFolder,strcat(subject,'_fiducial_rhj.jpeg')))
   saveas(lhj_fig,fullfile(outputFolder,strcat(subject,'_fiducial_lhj.jpeg')))
    close all;
   
    save(fullfile(outputFolder,strcat(subject,'_mri_ctf.mat')),'mri_ctf');
    save(fullfile(outputFolder,strcat(subject,'_ctf_fiducials.mat')),'ctf_fiducials');

end




%=========================================================================
%% this part can be looped through after all fiducials are marked, as it doesn't require user
%interface
badsubs={}; fin =[]; count = 0;
for subi = 1:length(subjectList)
    mytic = tic;
    subject = subjectList(subi);
    if subject<10
        subject = ['0',num2str(subject)];
    else
        subject = num2str(subject);
    end
    subject = ['S',subject];
    subject = string(subject);
    
    try
    folderi = find(contains({folderList.name},subject) & [folderList.isdir] ==1); %find subject folder in list,there should only be one
    if length(folderi) >1; error('too many folders found for this subject');end

    outputFolder = strcat(savePath,'\',folderList(folderi).name,'\Head Model');
    cd(outputFolder)
    load(strcat(subject,'_mri_acpc_rs.mat'))

    %% 4. Segment the MRI
    if exist([strcat(subject,'_seg_i.mat')],'file')~=2 || exist([strcat(subject,'_mri_segmented.mat')],'file') ~=2
        % Using Finite Element Method (FEM)
        cfg           = [];
        cfg.output         = {'scalp','skull','csf','gray','white'};
        cfg.brainsmooth    = 1;
        cfg.scalpthreshold = 0.11;
        cfg.skullthreshold = 0.15;
        cfg.brainthreshold = 0.15;
        mri_segmented_5_compartment = ft_volumesegment(cfg, mri_acpc_rs);
        seg_i = ft_datatype_segmentation(mri_segmented_5_compartment,'segmentationstyle','indexed');
        
        if plotStuff == 1
        cfg              = [];
        cfg.funparameter = 'tissue'; %%NJ changed from 'seg' to 'tissue' bc seg wasn't a field
        cfg.funcolormap  = gray(5); % distinct color per tissue
        cfg.location     = 'center';
        cfg.atlas        = seg_i;    % the segmentation can also be used as atlas
        ft_sourceplot(cfg, seg_i);
        end

        cd(outputFolder)
        save([strcat(subject,'_seg_i')], 'seg_i');
        
        
        % copy the anatomy into the segmented mri
        mri_segmented = mri_segmented_5_compartment;
        mri_segmented.anatomy = mri_acpc_rs.anatomy;

        if plotStuff ==1
        cfg = [];
        cfg.funparameter = 'gray';
        ft_sourceplot(cfg, mri_segmented);
        %print -dpng natmeg_dip_segmented_brain.png
        cfg.funparameter = 'white';
        ft_sourceplot(cfg, mri_segmented);
        
        cfg.funparameter = 'skull';
        ft_sourceplot(cfg, mri_segmented);
        
        
        cfg.funparameter = 'scalp';
        ft_sourceplot(cfg, mri_segmented);
        end
        
        cd(outputFolder)
        save([strcat(subject,'_mri_segmented')], 'mri_segmented');
    else
        cd(outputFolder)
        seg_mri = [strcat(subject,'_seg_i.mat')];
        subject_mri_segmented = [strcat(subject,'_mri_segmented.mat')];
        load(seg_mri)
        load(subject_mri_segmented)
    end


   saveas(gcf,fullfile(outputFolder,strcat(subject,'_segmentation.jpg')))
   close;
    %% 5B. Create the mesh
    % NOTE: you might have to copy simbio folder and add it to
    % ...\fieldtrip-20190920\external if simbio is not being recognized and you
    % get this error  "Error using ft_hastoolbox (line 496)
    % the SIMBIO toolbox is not installed"
    %6/15/21 I think fieldtrip is expected to be in a folder higher than
    %matlab, not within EEGlab plugins
    %
    plotStuff = 0; 
    if exist(strcat(subject,'_mesh_fem.mat'),'file') ~=2
    cfg.shift  = 0.3;
    cfg.method = 'hexahedral';
    tic;
    mesh_fem = ft_prepare_mesh(cfg,mri_segmented);
    toc;
    cd(outputFolder)
    save([strcat(subject,'_mesh_fem')], 'mesh_fem');
    else
        load(strcat(subject,'_mesh_fem.mat'))
    end
    %visualize mesh and electrodes
    if plotStuff ==1
%     figure; ft_plot_mesh(mesh_fem, 'surfaceonly', 'yes'); % NJ added
    figure;
    hold on
    ft_plot_mesh(mesh_fem,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'face alpha',0.7)
    camlight;
    end
    % % load elec; %load the electrodes
    % ft_plot_sens(elec,'style', '*g');
    
    %% 6 Align Electrodes
    % Co-register MRI to  based on fiducials (nasion and left and right pre-auricular points)
    %load electrode locations
    cd(electrodelocs_path)
    EEGchanfile = strcat('PHD',num2str(subjectList(subi)+14),'_',subject,'_chanlocs.txt');
    chanlocs = readtable(EEGchanfile);% Same output text file from getchalocs.
    chanlocs.Properties.VariableNames = {'labels','X','Y','Z'};
    elec.chanpos(:,1) = [chanlocs.X];
    elec.chanpos(:,2) = [chanlocs.Y];
    elec.chanpos(:,3) = [chanlocs.Z];
    elec.elecpos = elec.chanpos;
    elec.label(:,1) = [chanlocs.labels]';
    
    %load fiducial locations
    cd(outputFolder)
    load(strcat(subject,'_ctf_fiducials.mat'));

    %grab fiducial locations w/o transform
    nas = ctf_fiducials.nas;
    lhj = ctf_fiducials.lpa;
    rhj = ctf_fiducials.rpa;

    %transform fids to to headmodel coordinate space
    vox2head = mri_segmented.transform; %same transform as mr_realign.transformprig
    lhj = ft_warp_apply(vox2head, lhj,'homogenous');
    rhj = ft_warp_apply(vox2head, rhj,'homogenous');
    nas = ft_warp_apply(vox2head, nas,'homogenous');

    % create a structure similar to a template set of electrodes
    fid = [];
    fid.chanpos       = [nas; lhj; rhj];% head coordinates of fiducials
    fid.elecpos       = fid.chanpos;
    fid.label         = {'nas','lhj','rhj'};    % use the same labels as those in elec
    fid.unit          = 'mm';

    % alignment with fiducial method
    cfg               = [];
    cfg.method        = 'fiducial';
    cfg.elec          = elec;                  % the electrodes we want to align
    cfg.template      = fid;                   % the template we want to align to
    cfg.fiducial      = {'nas', 'lhj', 'rhj'};  % labels of fiducials in fid and in elec
    elec_aligned      = ft_electroderealign(cfg);
    
    %plot
    figure;
    ft_plot_mesh(mesh_fem,'surfaceonly','yes','vertexcolor','none','edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    camlight
    ft_plot_sens(elec_aligned,'style','.k');
    ft_plot_sens(fid,'style','*m');
    lgd = legend('Head model', 'EEG electrodes','Fiducials in MRI');
    lgd.FontSize = 18;
    view(90,0)  
    %save fig
    saveas(gcf,fullfile(outputFolder,strcat(subject,'_alignment_cor.jpeg')))
    view(0,90)
    saveas(gcf,fullfile(outputFolder,strcat(subject,'_alignment_ax.jpeg')))
    view(0,0)
    saveas(gcf,fullfile(outputFolder,strcat(subject,'_alignment_sag.jpeg')))
    close;


    %save realigned electrode locations and locations as txt file w/o fidciuals
    cd(outputFolder)
    save([strcat(subject,'_elec_aligned')], 'elec_aligned')

    %don't think I need this following part,using elec_aligned works:
    %===================================================
%     chanlocs_realigned_nofid = table(elec_aligned.label(1:end-3,:) ,elec_aligned.chanpos(1:end-3,1),elec_aligned.chanpos(1:end-3,2), elec_aligned.chanpos(1:end-3,3));
%     writetable(chanlocs_realigned_nofid,'chanlocs_realigned_nofid.txt','Delimiter',' ','WriteVariableNames',false)
%     chanlocs_realigned_nofid = readlocs('chanlocs_realigned_nofid.txt','format',{'labels','X','Y','Z'}); %load channel locations from txt file
%     save chanlocs_realigned_nofid chanlocs_realigned_nofid
    %==================================================
    %other electrode alignment methods:
    % % alignment with headshape method
    % cfg               = [];
    % cfg.method        = 'headshape';
    % cfg.headshape     = mesh_fem;    % the template we want to align to
    % cfg.elec          = elec;        % the electrodes we want to align
    % % cfg.template      = fid;
    % cfg.fiducial      = {'nas', 'lhj', 'rhj'};  % labels of fiducials in fid and in elec
    % elec_aligned      = ft_electroderealign(cfg);
    %%interactive electrode alignement
    % cfg          = [];
    % cfg.method   = 'interactive';
    % cfg.elec     = elec_aligned;
    % cfg.headshape = headmodel_fem;
    % elec_aligned = ft_electroderealign(cfg);
    
    %% 7. Create the headmodel (1-2 mins)
    cd(strcat(fieldtrip_path,'\external\simbio'))% path to simbio folder (some simbio functions weren't working unless I did this)
    cfg = [];
    cfg.method = 'simbio';
    %Search for matching tissue conductivity values (recommend)
    %conductivity order should match tissue label order
    % Create reference index: skin = 0.43, skull=0.01, CSF = 1.79, GM = 0.33, WM = 0.14
    conductivity_ref(1).tissue = 'scalp';
    conductivity_ref(2).tissue = 'skull';
    conductivity_ref(3).tissue = 'csf';
    conductivity_ref(4).tissue = 'gray';
    conductivity_ref(5).tissue = 'white';
    conductivity_ref(1).standard_cond = 0.43;
    conductivity_ref(2).standard_cond = 0.01;
    conductivity_ref(3).standard_cond = 1.79;
    conductivity_ref(4).standard_cond = 0.33;
    conductivity_ref(5).standard_cond = 0.14;
    % for each tissue label, find it's matching conductivity value in the reference index
    for i = 1:length(mesh_fem.tissuelabel)
        idx = find(strcmp({conductivity_ref.tissue},mesh_fem.tissuelabel{i}));
        cfg.conductivity(i) = conductivity_ref(idx).standard_cond;
    end
    %%or hardcode conductivity values (don't recommend, order might change)
    %cfg.conductivity = [1.79 0.33 0.43 0.01 0.14 ]; % same as tissuelabel in vol_simbio, {'scalp','skull','csf','gray','white'};
    cfg.tissuelabel = mesh_fem.tissuelabel; %NJ added
    headmodel_fem = ft_prepare_headmodel(cfg, mesh_fem);
    cd (outputFolder)
    save([strcat(subject,'_headmodel_fem')], 'headmodel_fem');


    %% 
    %disp('ADD MRI NORM HERE SO DIPOLE MOMENTS CAN BE RECOMPUTED ON HPG ') %%to do 
    %%
    close all;
    clearvars -except subi badsubs subject fieldtrip_path plotStuff electrodelocs_path subjectList savePath seg_MRI_dir folderList fin count mytic

    catch ME
         warning('Skipping this subject due to issues: ')
         fprintf('%s\n',ME.identifier)
         badsubs{end+1} = subject;
    end
      %estimate time remaining
        fprintf('\nFinished file %i/%i\n', subi,length(subjectList));
        t_remaining(mytic,fin,count,length(subjectList))
        close all;
    
    end
    
    %% Visualize the headmodel and the electrodes (it might take time and memory)
    % csf: 1, gm: 2, scalp: 3, skull: 4, wm: 5
    % doesn't work for me, get error "Undefined function 'mesh2edge' for input
    % arguments of type 'struct'"
    % ts = 3;
    % figure;
    % mesh2 =[];
    % mesh2.hex = headmodel_fem.hex(headmodel_fem.tissue==ts,:); %mesh2.hex(1:size(mesh2.hex),:);
    % mesh2.pos = headmodel_fem.pos;
    % mesh2.tissue = headmodel_fem.tissue(headmodel_fem.tissue==ts,:); %mesh.tissue(1:size(mesh2.hex),:);
    %
    % mesh_ed = mesh2edge(mesh2);
    % patch('Faces',mesh_ed.poly,...
    %   'Vertices',mesh_ed.pos,...
    %   'FaceAlpha',.5,...
    %   'LineStyle','none',...
    %   'FaceColor',[1 1 1],...
    %   'FaceLighting','gouraud');
    %
    % xlabel('coronal');
    % ylabel('sagital');
    % zlabel('axial')
    % camlight;
    % axis on;
    %
  %  figure;
    %cfg = [];
   % ft_plot_headmodel(headmodel_fem);
    
    % % Steps 8 & 9 are done in DIPFIT, so don't need to compute here---
    % % recommend doing here though because dipfit is dumb and recomputes things
    % % it doesn't need to. It would take about 30days/subject at the rate dipfit
    % % goes currently FYI. Had to modify functions within pop_multifit to allow
    % % for other inputs like sourcemodel and electrode file .mat files
    % %% 8. Create the sourcemodel
    % cfg = [];
    % cfg.resolution = 7.5;
    % cfg.threshold = 0.1;
    % cfg.smooth = 5;
    % cfg.headmodel = headmodel_fem;
    % cfg.inwardshift = 1; %shifts dipoles away from surfaces
    % sourcemodel = ft_prepare_sourcemodel(cfg);
    % % Visualize the sourcemodel
    % figure, ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:))
    % hold on, ft_plot_mesh(mesh_fem,'surfaceonly','yes','vertexcolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.1,'edgealpha',0.1)
    % ft_plot_sens(elec_aligned, 'style', '.g');
    % hold off
    % save sourcemodel sourcemodel
    %
    % %% Test atlas - doesn't work
    % % atlas = ft_read_atlas(strcat(fieldtrip_path,'/template/atlas/afni/TTatlas+tlrc.HEAD'));
    % % atlas = ft_convert_units(atlas,'mm');
    % %
    % % cfg = [];
    % % cfg.atlas      = atlas;
    % % cfg.roi        = atlas.brick1label  ;  % here you can also specify a single label, i.e. single ROI
    % % cfg.inputcoord = 'mni';
    % % mask           = ft_volumelookup(cfg, sourcemodel);
    %
    % %% 9 Compute the leadfield
    % % Note- this step takes a while! When I ran 128 channels, ft_prepare_leadfield took ~8.5 hrs
    % %compute the transfer matrix. When using parfor, 2 workers, took ~5 hrs. Makes my computer freeze for a couple minutes, but always comes around
    % %:)
    % tic;
    % [headmodel_fem_tr, elec] = ft_prepare_vol_sens(headmodel_fem, elec_aligned);
    % toc;
    % cd(outputFolder); save('headmodel_fem_tr', 'headmodel_fem_tr', '-v7.3')
    %
    % cfg = [];
    % cfg.sourcemodel = sourcemodel; %NJ changed cfg.grid to cfg.sourcemodel
    % cfg.headmodel= headmodel_fem_tr;
    % cfg.elec = elec_aligned;
    % cfg.reducerank = 3;%modify the leadfields by reducing the rank (i.e. remove the weakest orientation),default = 3 for EEG
    % leadfield = ft_prepare_leadfield(cfg);
    % cd(outputFolder);save('leadfield','leadfield', '-v7.3')
    %
    %
    %#####################################################################################
    %				Other optional steps
    %#####################################################################################
    %     %% 2. Normalize MRI (if not done already) -- alternatively normalize
%     dipole coordinates after making head model. THis solves the issues of
%     some MRIs gettting "cropped" after normalization. 
%     switch inputMRItype
%         case {'raw'}
%             % 2. Normalize original MRI using fieldtrip (raw_mri = 'SBP_S18_MRI_WIP_T1W_3D_MPRAGE_SENSE_1_1.nii';
%             %test template comparison mri_temp_file = 'C:\Users\jacobsen.noelle\Desktop\fieldtrip-20210614\template\anatomy\single_subj_T1_1mm.nii';
%             %mri_template = ft_read_mri(mri_temp_file);
%             disp('normalizing MRI')
%             cfg = [];
%             mri_norm = ft_volumenormalise(cfg,mri);
%             mri_norm = ft_determine_coordsys(mri_norm); %verify coordinate system is labeled correctly
%             if plotstuff ==1%plot
%                 cfg = [];
%                 ft_sourceplot(cfg, mri_norm); %when you click on anterior commissure, coordinates should be [0 0 0] (or close)
%             end
%             cd(outputFolder)
%             save([strcat(subject,'mri_norm')], 'mri_norm');
%             cfg = [];
%             cfg.parameter = 'anatomy';
%             cfg.filename = [subject,'_mri_norm'];
%             cfg.filetype = 'nifti';
%             cfg.datatype = 'double';
%             ft_volumewrite(cfg, mri_norm)
%         case {'norm'}
%             mri_norm = mri;
%         otherwise
%             disp('Invalid input MRI type')
%     end
    % %% 3. Realign MRI- don't do it this way unless you want output to be in
    % CTF coordinate system... will have issues with group analysis later
    % unless you can figure out how to transform dipole positions to the
    % template MNI space
    % %realign to ctf coordinate system by marking fiducials (nas, lhj, rhj)
    % %use same fiducials marked in EEG system (using lhj/rhj, instead of
    % %lpa/rpa)
    % cfg = [];
    % cfg.method = 'interactive';
    % cfg.coordsys = 'ctf';
    % mri_realigned = ft_volumerealign(cfg, mri_resliced);
    % %plot
    % cfg = [];
    % ft_sourceplot(cfg, mri_realigned);
    % cd(outputFolder)
    % save mri_realigned mri_realigned
    % % save the resliced anatomy in a FreeSurfer compatible format
    % cfg             = [];
    % cfg.filename    =strcat(extractBefore(mrifile.name,'.nii'), '_mri_realigned');
    % cfg.datatype = 'double'; %datatype of mri
    % cfg.filetype    = 'nifti_spm'; %save .nii using spm12
    % cfg.parameter   = 'anatomy'; %anatomical mri
    % ft_volumewrite(cfg, mri_realigned);

    %% 
    
% %     %transform fiducial coordinates to new coordinate system
% %     switch outputCoordSys
% %         case {'MNI'}
% %             mri = mri_norm;
% %             %use original transformation from original MRI to MNI to bring fid to
% %             %normalize coord space
% %             vox2head = mri_norm.transform; %same transform as mri_realign.transformprig
% %             nas = ft_warp_apply(vox2head, nas, 'homogenous');
% %             lhj = ft_warp_apply(vox2head, lhj, 'homogenous');
% %             rhj = ft_warp_apply(vox2head, rhj, 'homogenous');
% %         case {'CTF'}
% %             mri = mri_realigned;
% %             % save the realigned anatomy in a FreeSurfer compatible format
% %             cfg             = [];
% %             cfg.filename    =strcat(extractBefore(mrifile.name,'.nii'), '_mri_realigned');
% %             cfg.datatype = 'double'; %datatype of mri
% %             cfg.filetype    = 'nifti_spm'; %save .nii using spm12
% %             cfg.parameter   = 'anatomy'; %anatomical mri
% %             ft_volumewrite(cfg, mri_realigned);
% %             % transform fiducial locations to CTF coordinate space
% %             vox2head = mri_realigned.transform; %same transform as mr_realign.transformprig
% %             nas = ft_warp_apply(vox2head, nas, 'homogenous');
% %             lhj = ft_warp_apply(vox2head, lhj, 'homogenous');
% %             rhj = ft_warp_apply(vox2head, rhj, 'homogenous');
% %         otherwise
% %             disp('Invalid output coordinate system')
% %     end
%     
%     %verify fiducials are marked properly (look at cross hair location
%     cfg = [];
%     cfg.location = nas;
%     ft_sourceplot(cfg, mri_ctf); %
%     title('Nasion location'); fig = gcf; fig.OuterPosition = [2 300 500 500];
%     cfg.location = lhj;
%     ft_sourceplot(cfg,mri_ctf); %
%     title('LHJ location');fig = gcf; fig.OuterPosition = [502 300 500 500];
%     cfg.location = rhj;
%     ft_sourceplot(cfg, mri_ctf); %
%     title('RHJ location'); fig = gcf; fig.OuterPosition = [1002 300 500 500];
%     
%     prompt = 'No, remark';
%     while ~strcmp(prompt,'Yes')
%         %asks if fiducial locations are okay
%         prompt = questdlg('Are the fiducials in the correct location? ', 'Verify fiducial location','Yes', 'No, remark','Plot in template MRI to verify','Yes');
%         %response
%         switch prompt
%             case 'Yes'
%                 disp([' Everything looks good!']);
%             case 'No, remark'
%                 disp(['Remark fiducial locations']);
%                 cfg = [];
%                 cfg.method = 'interactive';
%                 cfg.coordsys = 'ctf'; %spm doesn't work bc won't let me mark nas, need to mark in this coodsys
%                 mri_realigned = ft_volumerealign(cfg, mri_norm); %you won't actually use this mri, just getting fiducial locations
%                 %grab fiducial locations you just marked using ft_volume realign
%                 nas = mri_realigned.cfg.fiducial.nas;
%                 lhj = mri_realigned.cfg.fiducial.lpa;
%                 rhj = mri_realigned.cfg.fiducial.rpa;
%                 %transform to coordinate system identified in above
%                 %case/switch as vox2head variable
%                 nas = ft_warp_apply(vox2head, nas, 'homogenous');
%                 lhj = ft_warp_apply(vox2head, lhj, 'homogenous');
%                 rhj = ft_warp_apply(vox2head, rhj, 'homogenous');
%                 %verify fiducials are marked properly (look at cross hair location
%                 cfg = []; cfg.location = nas;
%                 ft_sourceplot(cfg, mri); title('Nasion location');fig = gcf; fig.OuterPosition = [2 300 500 500];
%                 cfg.location = lhj;
%                 ft_sourceplot(cfg, mri); title('LHJ location'); fig = gcf; fig.OuterPosition = [502 300 500 500];
%                 cfg.location = rhj;
%                 ft_sourceplot(cfg, mri); title('RHJ location'); fig = gcf; fig.OuterPosition = [1002 300 500 500];
%             case 'Plot in template MRI to verify'
%                 %verify fiducials are marked properly (look at cross hair
%                 %location in mri template file
%                 disp(['Plotting on MNI template. Fiducial locations should be approximately the same']);
%                 mri_temp_file = 'C:\Users\jacobsen.noelle\Desktop\fieldtrip-20210614\template\anatomy\single_subj_T1_1mm.nii';
%                 disp(mri_temp_file)
%                 mri_template = ft_read_mri(mri_temp_file);
%                 cfg = []; cfg.location = nas;
%                 ft_sourceplot(cfg, mri_template); title('Nasion location'); fig = gcf; fig.OuterPosition = [2 300 500 500];
%                 cfg.location = lhj;
%                 ft_sourceplot(cfg, mri_template); title('LHJ location');fig.OuterPosition = [502 300 500 500];
%                 cfg.location = rhj;
%                 ft_sourceplot(cfg, mri_template); title('RHJ location'); fig = gcf; fig.OuterPosition = [1002 300 500 500];
%                 prompt = questdlg('Are the fiducials in the correct location? ','Verify fiducial location', 'Yes', 'No, remark','Plot in template MRI to verify','Yes');
%         end
%     end
%     close all;
%     % create a structure similar to a template set of electrodes
%     fid = [];
%     fid.chanpos       = [nas; lhj; rhj];% head coordinates of fiducials
%     fid.elecpos       = fid.chanpos;
%     fid.label         = {'nas','lhj','rhj'};    % use the same labels as those in elec
%     fid.unit          = 'mm';
%     cd(outputFolder)
%     save([strcat(subject,'_fidLocsMRI')], 'fid');