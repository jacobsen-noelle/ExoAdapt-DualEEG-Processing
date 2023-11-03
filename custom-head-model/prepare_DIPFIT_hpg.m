%Prepare DIPFIT for Hipergator

%folders
maindir = fullfile('\blue','dferris'); %group directory name
mydir = fullfile(maindir,'jacobsen.noelle','DIPFIT'); %my directory with run_dipfit_hpg.m
DIPFITdir = uigetdir(fullfile(mydir),'Select folder with EEG datasets to undergo dipole fitting');
folder = split(extractAfter(DIPFITdir,'Z:\'),'\');

%parameters
num_time = '00:45:00'; % Wall time hh:mm:ss , 45 mins for 30 sub
RV_thresh = 15; %residual variance threshold (default is 100 to fit all dipoles, typically set to 15 in EEG analysis)
customHDMflag = 0; %1 = using a custom head model,0=standard BEM
% setup path names
DIPFITdir_local = DIPFITdir;
DIPFITdir_unix = ['/blue/dferris/jacobsen.noelle/DIPFIT/',folder{end},'/'];
mydir_unix = char(join(split(mydir,"\"),"/"));


%hipergator_dir = AMICAdir_unix;
num_nodes = 1; % How many nodes to request %%changed from 128 to 264 for big data
num_tasks = num_nodes; % Number of MPI jobs
num_procs = 1; % number of cores per job
num_mem =8000; % memory per cpu used (in MB) (set to 175)
num_tasks_per_node = 1; %ryan edit: noticed the code ran best when you only used a single CPU from each node so i would generally set this to 1
qos_name = 'dferris-b'; % Use 'dferris' or 'dferris-b' (burst allocation)

%% Write .sh file
fid = fopen(fullfile(DIPFITdir_local, 'run_dipfit.sh'),'w');
fprintf(fid,['#!/bin/sh\n']);
fprintf(fid,['#SBATCH --job-name=Dipfit_%s # Job name\n'],[folder{end}]); %Don't have spaces in job name
fprintf(fid,['#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n']);
fprintf(fid,['#SBATCH --mail-user=jacobsen.noelle@ufl.edu  # Where to send mail\n']);	%CHANGE EMAIL
fprintf(fid,['#SBATCH --nodes=%d                    # Use one node\n'],num_nodes);
fprintf(fid,['#SBATCH --ntasks=%d                   # Run a single task\n'],num_tasks);	
fprintf(fid,['#SBATCH --ntasks-per-node=%d          # number of tasks per node\n'],num_tasks_per_node); %ryan edit: added this to have more control over how tasks are split across nodes
fprintf(fid,['#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n'],num_procs);
fprintf(fid,['#SBATCH --mem-per-cpu=%dmb                  # Total memory limit\n'],num_mem);
fprintf(fid,['#SBATCH --time=%s              # Time limit hrs:min:sec\n'],num_time);
fprintf(fid,['#SBATCH --output=%sDIPFIT_OUT.out     # Standard output and error log\n'], DIPFITdir_unix);
fprintf(fid,['#SBATCH --account=dferris	     # Account name\n']);
fprintf(fid,['#SBATCH --qos=%s		     # Quality of service name\n\n'],qos_name);

fprintf(fid,['module load matlab/2018a \n']);
fprintf(fid,['cd %s\n'], char(join(split(mydir,"\"),"/")));
fprintf(fid,['matlab -nodisplay -r ''run_dipfit_hpg ''%s'' %i %i; exit;'' '],folder{end},customHDMflag,RV_thresh);
fclose(fid);  

disp(['cd /blue/dferris/jacobsen.noelle/DIPFIT/',folder{end}]);
disp('sbatch run_dipfit.sh')