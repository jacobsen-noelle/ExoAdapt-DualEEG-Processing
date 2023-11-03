%Compute_headmodel_func
% previously "makeheadmodel_hpg
% mydir		parent directory containing subject folders e.g jacobsen.noelle/CustomHeadModel/Data/BATCH-XXX/ which contains folders for each subject

function Compute_headmodel_hpg(mydir)

folderList = dir(fullfile(mydir));

for folderi = 1:size(folderList,1)
	if ~contains(folderList(folderi).name,'.') && folderList(folderi).isdir == 1
		folder = strcat(mydir,'/',folderList(folderi).name);
	%addpath('/blue/dferris/jacobsen.noelle/fieldtrip-20210614') % add fieldtrip path located on HPG,compile functin
	%     does not support addpath
	%ft_defaults;
	%mydir = strcat('/blue/dferris/jacobsen.noelle/CustomHeadModel/',folder); %directory containing sourcemodel.mat, headmodel_fem.mat, elec_aligned.mat
	cd(folder)

	%load variables
	%load ('sourcemodel.mat','-mat', 'sourcemodel');
	load('vol.mat','-mat', 'vol');
	load('elec_aligned.mat', '-mat', 'elec_aligned');
	%disp(sourcemodel)
	%disp(headmodel_fem)
	disp(elec_aligned)
	elec_aligned    = rmfield(elec_aligned,'tra');% Remove this field to force average referencing of leadfield matrix
	
	
	
	%% 9 Compute the leadfield 
	% Note- this step takes a while! When I ran 128 channels, ft_prepare_leadfield took ~8.5 hrs
	%compute the transfer matrix. When using parfor, 2 workers, took ~5 hrs. Need to use parallel processing
	tic;
	parpool(20)
	[headmodel_fem_tr, elec_aligned_fem] = ft_prepare_vol_sens(vol, elec_aligned); % took ~4hrs
	hmend = toc;
	save('headmodel_fem_tr', 'headmodel_fem_tr', '-v7.3')
	save('elec_aligned_fem.mat','elec_aligned_fem'); %is this neccessary?

	 try
			p = gcp('nocreate');
			disp(p.NumWorkers)
	 end

	delete(p) % Turn parallel processing off to avoid out-of-memory issue

 %Create the sourcemodel
    cfg = [];
    cfg.resolution = 7.5;
    cfg.threshold = 0.1;
    cfg.smooth = 5;
    cfg.headmodel = headmodel_fem_tr;
	cfg.unit        = 'mm';
    cfg.inwardshift = 1; %shifts dipoles away from surfaces
    sourcemodel = ft_prepare_sourcemodel(cfg);
    save('sourcemodel.mat','sourcemodel')
	
	tic;
	cfg = [];
	cfg.sourcemodel = sourcemodel; %NJ change cfg.grid to cfg.sourcemodel
	cfg.headmodel= headmodel_fem_tr;
	cfg.elec = elec_aligned;
	cfg.reducerank = 3;%modify the leadfields by reducing the rank (i.e. remove the weakest orientation),default = 3 for EEG
	leadfield = ft_prepare_leadfield(cfg); %took ~3hrs
	save('leadfield','leadfield', '-v7.3')
	lfend = toc;
	
	 try
			p = gcp('nocreate');
			disp(p.NumWorkers)
	 end

	delete(p) % Turn parallel processing off to avoid out-of-memory issue
	end
end