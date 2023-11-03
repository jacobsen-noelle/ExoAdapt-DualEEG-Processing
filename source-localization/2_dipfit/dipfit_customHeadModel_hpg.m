%run_dipfit_customHeadmodel_hpg()
% Runs DIPFIT Autofit on Hipergator after headmodel settings have been saved
% EEG channels must be labeled as 'EEG' under EEG.chanlocs.type
% Results are output in a subfolder of the input folder called dipfit_DONE
% Dataset is stored as current filename with '_dipfit' appended
% Usage:
% >> [EEG] = run_dipfit_hpg(Folder,'key', 'val', ...);
%
% Required inputs:
%  mydir         - BATCH Folder containing subject subfolders with outputs from compute_headmodel_hpg.m
%                    presaved. New EEG dataset with dipole fitting results
%                    will be save in a new subfolder called
%                    Folder/dipfit_DONE
%
%  RV     		- [integer] for residual variance rejection threshold.
%               Recommended threshold = 15 (15% RV). 100 = fit all
%               dipole locations
% AMICApath		- path to BATCH folder containing subject subfolders with AMICA outputs- ** must match subject folder names in mydir
%					or leave empty if you already have comp.mat saved in mydir/SXX
%
%
function run_dipfit_hpg(mydir,AMICApath,RV)
RV = str2num(RV);
fprintf('Input folder: %s\n Custom head model\n RV=%i\n', mydir, RV');

folderList = dir(fullfile(mydir));

for folderi = 1:size(folderList,1)
    if ~contains(folderList(folderi).name,'.') && folderList(folderi).isdir == 1
        model_folder = strcat(mydir,'/',folderList(folderi).name);
        cd(model_folder)
        fileList = dir(fullfile(model_folder));
        compf_idx = find(contains({fileList.name},'comp') & contains({fileList.name},'.mat'));

        if size(compf_idx ==1)
            load(fileList(compf_idx).name); %load *comp.mat containing ICA data: comp = eeglab2fieldtrip(EEG, 'componentanalysis');
        else %load AMICA output .set file
            if ~isempty(AMICApath)
                disp('Loading comp data from AMICA output folder')
                comp_folder = strcat(AMICApath,'/',model_folder);
                if ~exist(comp_folder, 'dir') %check to see if AMICA folder exists
                    error('Cannot find corresponding AMICA folder with component data: %s', comp_folder)
                else
                    cd(comp_folder)
                    setfile = dir('*.set');
                    if size(setfile)==1
                        EEG = pop_loadset(dat);
                  		 comp = eeglab2fieldtrip(EEG, 'componentanalysis');
                    else
                        error('too many .set files in %s',comp_folder);
                    end
                end
              else
                    error('unable to load component data because AMICA filepath is empty')
              end
          end

            cd(model_folder)
            disp('loading head model...')
            load('headmodel_fem_tr.mat')
            disp('loading channel locations...')
            load('elec_aligned.mat')
            disp('loading leadfield...')
            load('leadfield.mat')

            if isfield(elec_aligned,'tra')
                elec_aligned    = rmfield(elec_aligned,'tra');% Remove this field to force average referencing of leadfield matrix
            end
            %% ========================================================================= %%
            %% Chang Liu's code-- use ft_dipolefitting instead of eeglab functions
            %% ========================================================================= %%

            %% Setup parallel processing
            mycluster = parcluster('local')
            thePool = parpool(16)
            %% Dipole fitting
            % coarse fit
            cfg = [];
            cfg.numdipoles    =  1;
            cfg.headmodel     = headmodel_fem_tr;
            cfg.sourcemodel   = leadfield;
            cfg.elec          = elec_aligned;
            cfg.dipfit.metric = 'rv';
            cfg.nonlinear     = 'no';
            cfg.component     = 1:size(comp.topo,2);
            % cfg.component     = 2;
            % for each component scan the whole brain with dipoles using FIELDTRIPs
            % dipolefitting function
            dipfit_fem_grid        = ft_dipolefitting(cfg,comp);
            save(fullfile(model_folder,'dipfit_fem_grid.mat'),'dipfit_fem_grid','-v7.3');

            % nonlinear fit
            cfg = [];
            cfg.numdipoles    =  1;
            cfg.headmodel     = headmodel_fem_tr;
            cfg.sourcemodel   = leadfield;
            cfg.elec          = elec_aligned;
            cfg.dipfit.metric = 'rv';
            cfg.nonlinear     = 'yes';
            cfg.component     = 1:size(comp.topo,2);
            % cfg.component     = [4 11:15];
            % for each component scan the whole brain with dipoles using FIELDTRIPs
            % dipolefitting function
            dipfit_fem        = ft_dipolefitting(cfg,comp);
            save(fullfile(model_folder,'dipfit_fem.mat'),'dipfit_fem','-v7.3');

            try
                p = gcp('nocreate');
                disp(p.NumWorkers)
            end
            delete(gcp('nocreate'))
            %% ========================================================================= %%
        end
end

% %maindir = '/blue/dferris/';
% %mydir = '/blue/dferris/jacobsen.noelle/DIPFIT/';
% %inputdir = strcat(mydir,folder);
% %addpath('/blue/dferris/jacobsen.noelle/eeglab2021.0'); %add path to eeglab functions
% %addpath('/blue/dferris/jacobsen.noelle/DIPFIT');

% %create output folder as subfolder of input folder
% %outputFolder = fullfile(mydir,'dipfit_DONE');
% outputFolder = mydir;
% if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
% mkdir(outputFolder);
% end

% %Load EEG datasets
% fileList = dir(fullfile(pwd, '*.set'));

% %load dataset and fit dipoles
% fprintf('Calculating Dipoles...\n');
% for i = 1:length(fileList)
% fileName = fileList(i).name;
% %[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;%lets try just using pop_loadset instead of starting up EEGLab

% EEG = pop_loadset('filename',fileName, 'filepath',inputdir);
% %adjust dipfit settings paths for hpg because it was saying files in
% %windows format didn't exist
% %     if (contains(EEG.dipfit.hdmfile,'Z:\') && customHDMflag ==1)
% %         newname = strcat(maindir,extractAfter(EEG.dipfit.hdmfile,'Z:\'));  %remove Z:\ and replace with main hpg directory
% %         EEG.dipfit.hdmfile  = join(split(newname,"\"),"/"); %change to unix
% %         disp(EEG.dipfit.hdmfile)
% %         %do same name change for MRI file
% %         newname = strcat(maindir,extractAfter(EEG.dipfit.mrifile,'Z:\'));  %remove Z:\ and replace with main hpg directory
% %         EEG.dipfit.mrifile  = join(split(newname,"\"),"/"); %change to unix
% %         disp(EEG.dipfit.mrifile)
% %         %do same name change for chanfile
% %         newname = strcat(maindir,extractAfter(EEG.dipfit.chanfile,'Z:\'));  %remove Z:\ and replace with main hpg directory
% %         EEG.dipfit.chanfile  = join(split(newname,"\"),"/"); %change to unix
% %         disp(EEG.dipfit.chanfile)
% %     end

% EEG = pop_multifit(EEG, [1:size(EEG.icaweights,1)] ,'threshold',RV,'rmout','on');
% EEG.setname = strcat(EEG.subject,' dipfit');
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_saveset( EEG, 'filename',strcat(extractBefore(EEG.filename,'.set'),'_dipfit'),'filepath',outputFolder,'version','7.3');

% end
fprintf('Finished dipole fitting');
end