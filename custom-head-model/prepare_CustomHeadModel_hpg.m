%Prepare Custom Head Model for Hipergator
clc; %clear all; close all;
%folders
maindir = '\\exasmb.rc.ufl.edu\blue\dferris\jacobsen.noelle\';; %group directory name, blue drive
mydir =  uigetdir(maindir,'Select custom head model parent folder on blue drive containing makeheadmodel_hpg'); %my directory with makeheadmodel_hpg.m
local_data_dir = 'R:\Ferris-Lab\jacobsen.noelle\Exo Adaptation\Data\raw_data\'; %...\SXX\Head Model contains elec_aligned.mat and headmodel_fem.mat
% CustomHDMdir = uigetdir(fullfile(mydir),'Select folder with compute_headmodel_hpg function');


%subjectList =  [18:44 46:48];
%outputFolder = strcat(mydir,'\Data\, datestr(now, 'yyyy-mm-dd')); %blue drive folder to copy data to
outputFolder = [maindir,'\CustomHeadModel\Data\2023-09-22-Exo'];

if ~exist(outputFolder, 'dir') %check to see if output folder exists, if not make new one
    mkdir(outputFolder);
end

%parameters
user_email = 'jacobsen.noelle@ufl.edu';
num_time = '03:00:00'; % Wall time hh:mm:ss , 50 mins for 1 sub
% setup path names
outputFolder_unix = [extractAfter(outputFolder,'\\exasmb.rc.ufl.edu\')]; 
outputFolder_unix  = char(join(split(outputFolder_unix ,"\"),"/"));


%hipergator_dir = AMICAdir_unix;
num_nodes = 1; % How many nodes to request %%changed from 128 to 264 for big data
num_tasks = num_nodes; % Number of MPI jobs
num_procs = 24; % number of cores per job
num_mem =8000; % memory per cpu used (in MB) (set to 175)
num_tasks_per_node = 1; %ryan edit: noticed the code ran best when you only used a single CPU from each node so i would generally set this to 1
qos_name = 'dferris-b'; % Use 'dferris' or 'dferris-b' (burst allocation)

localSubfolders = dir(local_data_dir);
%% copy data over to Blue drive
for s = 1:length(subjectList)
    subject = subjectList(s);
   if subject <10
        subject = ['0',num2str(subject)];
    else
        subject = num2str(subject);
    end
    subject = ['S',subject];
    subject = string(subject);
    subject_dir = strcat(outputFolder,'/',subject);

    if ~exist(subject_dir, 'dir') %check to see if output folder exists, if not make new one
    mkdir(subject_dir);
    end

   
    %copy vol/headmodel_fem and elec_aligned to blue drivefolder
	hm_folder = [local_data_dir,localSubfolders(contains({localSubfolders.name},subject) & [localSubfolders.isdir] ==1).name,'\Head Model'];
    cd(hm_folder)

    try
        load(strcat(subject,'_headmodel_fem.mat'));
        vol = headmodel_fem;
    catch
        load(strcat(subject,'vol.mat'));
    end
    load(strcat(subject,'_elec_aligned.mat'));
    cd(subject_dir)
    save('vol', 'vol');
    save('elec_aligned','elec_aligned');
end


%% Write .sh file
folder = split(extractAfter(outputFolder,'\\exasmb.rc.ufl.edu\'),'\');
fid = fopen(fullfile(outputFolder, 'run_compute_headmodel.sh'),'w');
fprintf(fid,['#!/bin/sh\n']);
fprintf(fid,['#SBATCH --job-name=HeadModel-%s # Job name\n'],[folder{end}]); %Don't have spaces in job name
fprintf(fid,['#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)\n']);
fprintf(fid,['#SBATCH --mail-user=%s  # Where to send mail\n'],[user_email])
fprintf(fid,['#SBATCH --nodes=%d                    # Use one node\n'],num_nodes);
fprintf(fid,['#SBATCH --ntasks=%d                   # Run a single task\n'],num_tasks);	
fprintf(fid,['#SBATCH --ntasks-per-node=%d          # number of tasks per node\n'],num_tasks_per_node); %ryan edit: added this to have more control over how tasks are split across nodes
fprintf(fid,['#SBATCH --cpus-per-task=%d            # Number of CPU cores per task\n'],num_procs);
fprintf(fid,['#SBATCH --mem-per-cpu=%dmb                  # Total memory limit\n'],num_mem);
fprintf(fid,['#SBATCH --time=%s              # Time limit hrs:min:sec\n'],num_time);
fprintf(fid,['#SBATCH --output=%s/diary.out     # Standard output and error log\n'], outputFolder_unix);
fprintf(fid,['#SBATCH --account=dferris	     # Account name\n']);
fprintf(fid,['#SBATCH --qos=%s		     # Quality of service name\n\n'],qos_name);

fprintf(fid,['module load matlab/2022a \n']);
%fprintf(fid,['cd %s\n'], char(join(split(mydir,"\"),"/")));
%fprintf(fid,['matlab -nodisplay -r ''makeheadmodel_hpg ''%s''; exit;'' '],folder{end});
fprintf(fid,'cd /blue/dferris/jacobsen.noelle/CustomHeadModel/1_compute_headmodel\n');
fprintf(fid,'./compute_headmodel_hpg "%s"', outputFolder_unix)
fclose(fid);  


