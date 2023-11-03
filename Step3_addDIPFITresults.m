%load dipfit results from HPG
%NOTE: dipole moments are currently inaccurate after MNI normalization-
%shouldn't affect anything other than figs w/ moments or if you cluster
%with orientation as a feature
mydir = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation';
eeglabpath = 'C:\Users\jacobsen.noelle\Desktop\eeglab2022.0';
codeRepo = 'R:\Ferris-Lab\jacobsen.noelle\Matlab_Scripts\ExoAdapt-EEG-Processing\';
addpath(genpath(codeRepo))
savePath = [mydir,'\Data\processed_data']; %subfolders will be created and added
%blueDrivePath = uigetdir(,'Choose a blue drive folder to parent folder containing subject subdir and dipfit data');% blue drive path to Data folder with ppt subfolders
blueDrivePath= '\\exasmb.rc.ufl.edu\blue\dferris\jacobsen.noelle\CustomHeadModel\Data\2023-09-22-Exo';
TITLE = 'Choose a folder with EEG files';
% EEGinputFolder = uigetdir(pwd,TITLE);
EEGinputFolder = [savePath,'\2023-09-22-Step2-AMICA']; %local folder with EEG sets to add dipfit to, likely your AMICA results folder


downsample_flag = 1; % [0|1] , downsample data
Fs = []; % downsampling value
subjects= [31]; %numeric range or leave empty to proceess all subjects in folder
rawdataFolder = [mydir,'\Data\raw_data']; %folder on R drive with subfolders containing subject data and MRIs
outputFolder = strcat(savePath,'\\',datestr(now, 'yyyy-mm-dd'),'_Step3-DIPFIT');
if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder);
end

%% Load files based on subject number
folderList = {};
if ~isempty(subjects)
    count = 1;
    for subi = subjects
        subject_num= num2str(subi);
        if (size(subject_num,2)==1)
            subject_num= strcat('0',subject_num);
        end
        folderList{count} = strcat(blueDrivePath,'\S', subject_num);
        count = count+1;
    end
else %load all folders in blueDrivePath
    mydirFolders = dir(fullfile(blueDrivePath));
    count = 1;
    for folderi = 1:size(mydirFolders,1)
        if ~contains(mydirFolders(folderi).name,'.') && mydirFolders(folderi).isdir == 1
            folderList{count} = strcat(blueDrivePath,'\', mydirFolders(folderi).name);
            count = count+1;
        end
    end
end



%% Load the custom MRI results
fin =[]; count = 0;
for folderi = 1:size(folderList,2)
    mytic =   tic;
    folder = folderList{folderi};
    subject = extractAfter(folderList{folderi},strcat(blueDrivePath,'\'));
    subject_num = extractAfter(subject,'S'); subject_num = str2num(subject_num);
    hmFolder = [rawdataFolder,'\PHD',num2str(subject_num+14),'_',subject,'\Head Model']; %head model folder
 
    %MRI path
    MRI_path = fullfile([rawdataFolder,'\PHD',num2str(subject_num+14),'_',subject,'\Head Model\']);
    mri_rs_filepath = fullfile( hmFolder,strcat(subject,'_MRI_acpc_rs.nii'));
    mri_rs_matfilepath = fullfile([ hmFolder,subject,'_mri_acpc_rs.mat']);

    %% loa eeg dataset
    if ~exist('ALLCOM')
        addpath(eeglabpath);
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    end
    EEGdatasets = dir(fullfile(EEGinputFolder));
    EEGfileIdx = find(contains({EEGdatasets.name},subject) & contains({EEGdatasets.name},'.set'));

    try
        EEG = pop_loadset('filename',EEGdatasets(EEGfileIdx).name,'filepath',EEGinputFolder);
    catch
        pause(60) %pause in case loading error is just temp network connection isssue
        EEG = pop_loadset('filename',EEGdatasets(EEGfileIdx).name,'filepath',EEGinputFolder);
    end
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    %load dipfit results
    cd(folder)
    clear dipfit_fem
    EEG.dipfit = [];
    fprintf('loading custom head model dipfit results: folder %s...\n',folderList{folderi});
    try
        disp('loading dipfit result info...'),load('dipfit_fem.mat');
    catch
        error('Could not find dipfit_fem.mat in folder')
    end
    % disp('loading headmodel info...'),load(headmodel_filepath)

    %%
    for i=1:length(dipfit_fem.component)
        EEG.dipfit.model(i).posxyz = dipfit_fem.dip(i).pos;
        EEG.dipfit.model(i).momxyz = reshape(dipfit_fem.dip(i).mom, 3, length(dipfit_fem.dip(i).mom)/3)';
        if ~isempty(dipfit_fem.dip(i).rv)
            EEG.dipfit.model(i).rv     = dipfit_fem.dip(i).rv;
        else
            EEG.dipfit.model(i).rv     = NaN;
        end
        EEG.dipfit.model(i).diffmap = dipfit_fem.Vmodel(i) - dipfit_fem.Vdata(i);
        EEG.dipfit.model(i).sourcepot = dipfit_fem.Vmodel(i);
        EEG.dipfit.model(i).datapot   = dipfit_fem.Vdata(i);

        dippos(i,:) = dipfit_fem.dip(i).pos;
    end
    %         [dipole.inside] = ft_inside_headmodel(dippos, headmodel_fem_tr);
    %         for i=1:length(dipfit_fem.component)
    %             EEG.dipfit_fem.model(i).inside = dipole.inside(i);
    %         end
    IC_RV = vertcat(EEG.dipfit.model.rv);


    %% Normalize MRI
     normalizeMRI = 1;
    if normalizeMRI %
        cd(hmFolder)
        %use mri_norm in head model folder if it exists
        disp("TEMP")
        temp = exist([subject,'_mri_norm.mat']);
        if temp ==2 %TEMP-- add back underscore before MRI
            mri_norm = ft_read_mri([subject,'_mri_norm.mat']);
        else
            load(mri_rs_matfilepath)
            cfg = [];
            cfg.nonlinear = 'yes';
            cfg.spmmethod = 'old'; % SPM new method takes forever, Mailing list argues that the spmmethod new works better
            %             cfg.write = 'yes';
            %             cfg.name = fullfile('R:\Ferris-Lab\liu.chang1\DATA\Headmodel_MOBI\normalize',[subjStr,'_norm']);
            mri_norm = ft_volumenormalise(cfg,mri_acpc_rs);

            % template location
            % R:\Ferris-Lab\liu.chang1\Code\spm12\toolbox\OldNorm\T1.nii
            cfg             = [];
            cfg.filename    = fullfile([hmFolder,subject,'_mri_norm']);
            cfg.filetype    = 'nifti';
            cfg.parameter   = 'anatomy';
            ft_volumewrite(cfg, mri_norm);
            cd(hmFolder)
            save([subject,'_mri_norm.mat'],'mri_norm')
            %ft_sourceplot(cfg,mri_norm)
        end
    end
    %% Convert the dipole location to MNI space
    dipfit_fem_pos = reshape([dipfit_fem.dip(:).pos],3,[])';
    dipfit_fem_mnipos = ft_warp_apply(mri_norm.params,ft_warp_apply(mri_norm.initial,dipfit_fem_pos), 'individual2sn');
    dipfit_fem_mni_voxinds = round(ft_warp_apply(pinv(mri_norm.transform), dipfit_fem_mnipos ));
    for i=1:length(dipfit_fem.component)
        dipfit_fem.dip(i).mnipos = dipfit_fem_mnipos(i,:);
        dipfit_fem.dip(i).mni_voxinds = dipfit_fem_mni_voxinds(i,:);
        EEG.dipfit.model(i).posxyz = dipfit_fem_mnipos(i,:);
        EEG.dipfit.model(i).voxinds = dipfit_fem_mni_voxinds(i,:);
    end
    %save other headmodel info
    EEG.dipfit.hdmfile = strcat(eeglabpath,'\plugins\dipfit\standard_BEM\standard_vol.mat');
    EEG.dipfit.mrifile = strcat(eeglabpath,'\plugins\dipfit\standard_BEM\standard_mri.mat');
    EEG.dipfit.chanfile = strcat(eeglabpath,'\plugins\dipfit\standard_BEM\elec\standard_1005.elc');
    EEG.dipfit.chansel = 1:length(EEG.chanlocs);
    EEG.dipfit.coordformat = 'MNI';
    EEG.dipfit.coord_transform = [0,0,0,0,0,0,1,1,1];

    EEG.dipfit_ACPC = dipfit_fem;
    EEG.dipfit.comments = {'Custom FEM head model--> dipole fitting--> warp coord to MNI headspace. Original results in ACPC space found in EEG.dipfit_ACPC. MOMENTS ARE INACCURATE-- were not recalcuated after normalizing to MNI'};
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG,EEG, CURRENTSET);
    EEG = eeg_checkset(EEG);
    %     eeglab redraw
    EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_CHMdipfit.set'),'filepath',outputFolder);

    %feedback
    %% estimate time remaining
     %estimate time remaining
     fprintf('\nFinished file %i/%i\n', folderi,size(folderList,2));
     t_remaining(mytic,fin,count,size(folderList,2))
    close all; clear mri_norm dipfit_fem
end


%             dipfit_fem_mni = dipfit_fem;


%   % nonlinear fit
%             cfg = [];
%             cfg.numdipoles    =  1;
%             cfg.headmodel     = headmodel_fem_tr;
%             cfg.sourcemodel   = leadfield;
%             cfg.elec          = elec_aligned;
%             cfg.dipfit.metric = 'rv';
%             cfg.nonlinear     = 'yes';
%             cfg.component     = 1:size(comp.topo,2);
%             % cfg.component     = [4 11:15];
%             % for each component scan the whole brain with dipoles using FIELDTRIPs
%             % dipolefitting function
%             dipfit_fem        = ft_dipolefitting(cfg,comp);
%fix dipole moments -- excerpt from ft_dipolefitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the model potential distribution and the residual variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%