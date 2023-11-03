%prepare_HPG_AMICA_func(EEG,fileNameNoExt,mainAMICAdir,AMICAdir_unix,params,emailStr)  
% creates all files necessary to run AMICA on hipergator
%
% Required inputs:
%   EEG                  - EEG dataset structures you want to save on M drive and run AMICA on via hipergator
%	fileNameNoExt 	     - what you want to save the EEG file as (string without .set extension)
%	mainAMICAdir		 - where you want to store the eeg data set and accompanying parameter files (string)  
%					       You just need to supply the general location (e.g. 'Z:\jacobsen.noelle'). It will
%					       search for or create an "AMICA" folder with a subfolder containing the date
%   outputdir_unix       - name of main folder where you want output files
%                          stored (string; e.g. '/blue/dferris/YOURFOLDER/AMICA/';
%   emailStr 			 - string with your email if you want to be notified when HPG processing is finished
% Optional inputs:
%	params				 - structure containing AMICA parameters with the following fields: 
%								num_time 			 - Wall time hh:mm:ss (string)
%								avgRefPCAReduction   - Number of components to reduce to using PCA. Possible values: 0,1,2,3 
%													    (0 if you reref w/o losing rank, 1 reref EEG, 2 if EEG and EMG separately avg ref,
%													     3 if additional rank was lost through rereferencing; default = 1);
%						  	    pcaOverride			 - %PCA value to use (default = 0)
%							    max_iter			 - maxium number of AMICA iterations (defualt = 1500)
%								channels			 - string of channel types to analyze.'EEG' - scalp channels only, 'EEG+EMG'- scalp and EMG channels,
%													   'EMG'- emg only, 'all'- all channels (default = 'EEG')
%						   leave empty (params = [];) if you want to use default values. 
%
% Example: 
%     params.num_time = '00:45:00'; % Wall time hh:mm:ss
%     params.avgRefPCAReduction = 0; %1 or 2 or maybe even 3 or more if you want to get crazy (2 if EEG and EMG separately avg ref, 1 otherwise); or 0 if for some reason you never average ref or did the verison of avg ref where you don't lose data rank
%     params.pcaOverride =0; %PCA value to use
%     params.channels = 'EEG';
%     mainAMICAdir = 'Z:\jacobsen.noelle';
%     outputdir_unix =  '/blue/dferris/jacobsen.noelle/AMICA/';
%     fileNameNoExt = extractBefore(EEG.filename,'.set');
%     emailStr = 'jacobsen.noelle@ufl.edu';
%     prepare_HPG_AMICA_func(EEG,fileNameNoExt,mainAMICAdir,AMICAdir_unix,params,emailStr)  

% Authors:
% Ryan Downey, Noelle Jacobsen,and other members of Human Neuromenchanics Lab
% University of Florida
% Last updated 2/8/2022

function prepare_HPG_AMICA_func(EEG,fileNameNoExt,mainAMICAdir,outputdir_unix,params,emailStr)  
%check parameters
if isempty(params)
    params.num_time = '00:45:00'; % Wall time hh:mm:ss
    params.pcaOverride =0 ; %PCA value to use
    params.max_iter = 2000; 
    params.channels = 'EEG'; 
else
    if ~isfield(params,'num_time')
        params.num_time = '00:45:00'; % Wall time hh:mm:ss
    end
    if ~isfield(params,'avgRefPCAReduction')
        params.avgRefPCAReduction = 0;
    end
    if ~isfield(params,'pcaOverride')
        params.pcaOverride =0 ; %PCA value to use
    end
    if ~isfield(params,'max_iter')
        params.max_iter = 2000; 
    end
     if ~isfield(params,'channels')
        params.channels = 'EEG'; %'EEG' - scalp channels only, 'EEG+EMG'- scalp and EMG channels,'EMG'- emg only, 'all'- all channels
    end
end

subDirNum = date;
shortfileName = EEG.subject;
AMICAdir_unix =[outputdir_unix,subDirNum,'/',shortfileName,'/'];
%% select channels
getchantypes;
    if isempty(neckEMG_chans) || isempty(EEG_chans)
        msg = 'Error : index is empty';
        error(msg)
    end
    
%select channels
if strcmp(params.channels, 'EEG')  %scalp only
    EEG = pop_select( EEG,'channel',sort([EEG_chans])); 
elseif strcmp(channels, 'EMG')  %scalp only
    EEG = pop_select( EEG,'channel',sort([neckEMG_chans])); 
elseif strcmp(channels, 'EEG+EMG')%scalp and EMG
    EEG = pop_select( EEG,'channel',sort([EEG_chans neckEMG_chans]));
else
   fprintf('Keeping all channel types')
end

num_chans = EEG.nbchan;
num_frames = length(EEG.times);

cd(mainAMICAdir)
AMICAdir_local = [strcat(mainAMICAdir,'\AMICA\',subDirNum,'\',shortfileName,'\')];
%AMICAdir_local = ['\AMICA\',num2str(subDirNum),'\',shortfileName,'\'];
if ~exist(AMICAdir_local , 'dir')
  mkdir(AMICAdir_local )
end

fprintf('\nParticipant %s', EEG.subject);
% Save float file for running ICA
tmpdata = EEG.data;
fprintf('\n');
disp(' Converting data to double...');
tmpdata = double(tmpdata);

%[SUCCESS,MESSAGE,MESSAGEID] = mkdir(AMICAdir_local);
%write file (make sure eeglab writes set and float not just set or you get error)
EEG = pop_saveset( EEG, 'filepath', AMICAdir_local, 'filename', EEG.filename, 'version', '7.3');

%% Save AMICA parameter file 
floatFileName = [strcat(fileNameNoExt,'.fdt')];
floatFileName = char(floatFileName);
num_nodes = 128; % How many nodes to request %%changed from 128 to 264 for big data
num_tasks = num_nodes; % Number of MPI jobs
num_procs = 1; % number of cores per job
num_mem =1000; % memory per cpu used (in MB) (set to 175)
num_tasks_per_node = 1; %ryan edit: noticed the code ran best when you only used a single CPU from each node so i would generally set this to 1
qos_name = 'dferris-b'; % Use 'dferris' or 'dferris-b' (burst allocation)
if params.pcaOverride ~= 0
    pca_keep = params.pcaOverride;
else
   % pca_keep = min([rank(tmpdata), EEG.nbchan-params.avgRefPCAReduction]);
   % % PCA value to use
     pca_keep= getrank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))), params.pcaOverride);
end
display(['Num ch = ', num2str(EEG.nbchan)]);
display(['PCA Reduction = ', num2str(pca_keep)]);
fid = fopen(fullfile(AMICAdir_local, 'input_hipergator.param'),'w');

fprintf(fid,['files ' AMICAdir_unix floatFileName '\n']);
fprintf(fid,['outdir ' AMICAdir_unix '\n']);
fprintf(fid,['num_models 1\n']);
fprintf(fid,['num_mix_comps 3\n']);
fprintf(fid,['pdftype 0\n']);
fprintf(fid,['block_size 128\n']); %noelle
fprintf(fid,['max_iter %d\n'],params.max_iter); %ryan
fprintf(fid,['num_samples 1\n']);
fprintf(fid,['data_dim %d\n'],num_chans);
fprintf(fid,['field_dim %d\n'],num_frames);
fprintf(fid,['field_blocksize 1\n']);
fprintf(fid,['share_comps 0\n']);
fprintf(fid,['share_start 100\n']);
fprintf(fid,['comp_thresh 0.990000\n']);
fprintf(fid,['share_iter 100\n']);
fprintf(fid,['lrate 0.100000\n']);
fprintf(fid,['minlrate 1.000000e-08\n']);
fprintf(fid,['lratefact 0.500000\n']);
fprintf(fid,['rholrate 0.050000\n']);
fprintf(fid,['rho0 1.500000\n']);
fprintf(fid,['minrho 1.000000\n']);
fprintf(fid,['maxrho 2.000000\n']);
fprintf(fid,['rholratefact 0.500000\n']);
fprintf(fid,['kurt_start 3\n']);
fprintf(fid,['num_kurt 5\n']);
fprintf(fid,['kurt_int 1\n']);
fprintf(fid,['do_newton 1\n']);
fprintf(fid,['newt_start 50\n']);
fprintf(fid,['newt_ramp 10\n']);
fprintf(fid,['newtrate 1.000000\n']);
fprintf(fid,['do_reject 0\n']);
fprintf(fid,['numrej 15\n']);
fprintf(fid,['rejsig 3.000000\n']);
fprintf(fid,['rejstart 1\n']);
fprintf(fid,['rejint 1\n']);
%note: there can be memory issues if your dataset is large and num_tasks is
%greater than 1. You can either 1) use PCA reduction on the data or 2) set 
%max_threads to 1 and set num_tasks to be a large number with one task per 
%node. AMICA runs quickly this way
%fprintf(fid,['max_threads %d\n'],num_tasks); 
fprintf(fid,['max_threads %d\n'],1); 
fprintf(fid,['writestep 10\n']);
fprintf(fid,['write_nd 0\n']);
fprintf(fid,['write_LLt 0\n']);
fprintf(fid,['decwindow 1\n']);
fprintf(fid,['max_decs 3\n']);
fprintf(fid,['update_A 1\n']);
fprintf(fid,['update_c 1\n']);
fprintf(fid,['update_gm 1\n']);
fprintf(fid,['update_alpha 1\n']);
fprintf(fid,['update_mu 1\n']);
fprintf(fid,['update_beta 1\n']);
fprintf(fid,['invsigmax 100.000000\n']);
fprintf(fid,['invsigmin 0.000000\n']);
fprintf(fid,['do_rho 1\n']);
fprintf(fid,['load_rej 0\n']);
fprintf(fid,['load_W 0\n']);
fprintf(fid,['load_c 0\n']);
fprintf(fid,['load_gm 0\n']);
fprintf(fid,['load_alpha 0\n']);
fprintf(fid,['load_mu 0\n']);
fprintf(fid,['load_beta 0\n']);
fprintf(fid,['load_rho 0\n']);
fprintf(fid,['load_comp_list 0\n']);
fprintf(fid,['do_mean 1\n']);
fprintf(fid,['do_sphere 1\n']);
fprintf(fid,['doPCA 1\n']);
fprintf(fid,['pcakeep %d\n'],pca_keep);
fprintf(fid,['pcadb 30.000000\n']);
fprintf(fid,['byte_size 4\n']);
fprintf(fid,['doscaling 1\n']);
fprintf(fid,['scalestep 1\n']);
fclose(fid);

%% Save hipergator job script
fid = fopen(fullfile(AMICAdir_local, 'runAMICA_hipergator.sh'),'w');
fprintf(fid,['#!/bin/sh\n']);
fprintf(fid,['#SBATCH --job-name=%s_ICA # Job name\n'],[EEG.subject]);
fprintf(fid,['#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n']);
fprintf(fid,['#SBATCH --mail-user=%s # Where to send mail\n'],[emailStr]);	%email to send updates to
fprintf(fid,['#SBATCH --nodes=%d                    # Use one node\n'],num_nodes);
fprintf(fid,['#SBATCH --ntasks=%d                   # Run a single task\n'],num_tasks);	
fprintf(fid,['#SBATCH --ntasks-per-node=%d          # number of tasks per node\n'],num_tasks_per_node); %ryan edit: added this to have more control over how tasks are split across nodes
fprintf(fid,['#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n'],num_procs);
fprintf(fid,['#SBATCH --mem-per-cpu=%dmb                  # Total memory limit\n'],num_mem);
fprintf(fid,['#SBATCH --time=%s              # Time limit hrs:min:sec\n'],params.num_time);
fprintf(fid,['#SBATCH --output=%sAmicaOUT.out     # Standard output and error log\n'], AMICAdir_unix);
fprintf(fid,['#SBATCH --account=dferris	     # Account name\n']);
fprintf(fid,['#SBATCH --qos=%s		     # Quality of service name\n\n'],qos_name);
%fprintf(fid,['#SBATCH --constraint=haswell|skylake|dhabi|broadwell|sandy-bridge']); 
%fprintf(fid,['#SBATCH --constraint=haswell|skylake|dhabi|broadwell|sandy-bridge']); 
fprintf(fid,['# Run your program with correct path and command line options\n']);
%fprintf(fid,['#SBATCH --partition=hpg2-compute \t# Run your program with correct path and command line options\n']); %noelle added 5/12/21
%fprintf(fid,['module load intel/2020 openmpi/4.0.3\n']); %module load intel openmpi
fprintf(fid,['module load ufrc\n']);
fprintf(fid,['module load intel/2020 openmpi/4.1.5\n']);
% fprintf(fid,['mpirun /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sParticipant%s_input_hipergator.param\n'],hipergator_dir, num2str(Participant));
% fprintf(fid,['srun --mpi=pmix /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir);
% fprintf(fid,['srun --mpi=pmix_v2 /ufrc/dferris/s.peterson/test/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir);
%fprintf(fid,['srun --mpi=pmix_v2 /ufrc/dferris/jacobsen.noelle/AMICA/AMICA_15/amica15ub %sinput_hipergator.param\n'],hipergator_dir); %%%only difference from theresas..  maybe cause of memory errors?
fprintf(fid,['srun --mpi=pmix_v3 /blue/dferris/jacobsen.noelle/AMICA/AMICA_15/amica15ub %sinput_hipergator.param\n'],AMICAdir_unix);
fclose(fid);

disp('Use MobaXTerm to submit your job (use right click to paste)');
disp(newline);
disp(['cd ' AMICAdir_unix]);
disp(['sbatch runAMICA_hipergator.sh']);
disp('squeue -A dferris');

end



% getrank() -- subfunction from pop_runica
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
% See also: runica(), binica(), jader(), fastica()
% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
function tmprank2 = getrank(tmpdata, pca_opt)
    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffman
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        if pca_opt ~= 0
            fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
        end
        tmprank2 = max(tmprank, tmprank2);
    end
end
            
            