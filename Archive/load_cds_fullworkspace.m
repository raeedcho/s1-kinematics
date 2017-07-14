%% Load data into CDS file

meta.lab=6;
meta.ranBy='ranByRaeed';
meta.monkey='monkeyHan';
meta.task='taskRW';
meta.array='arrayLeftS1Area2';
meta.folder='C:\Users\rhc307\Projects\limblab\data-preproc\MultiWorkspace\SplitWS\Han\20160322\';
meta.fname='Han_20160322_RW_PM_area2_002';
meta.mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Han\left_S1\SN 6251-001459.cmp';

% Make CDS files

cds = commonDataStructure();
cds.file2cds([meta.folder 'preCDS\' meta.fname],meta.ranBy,meta.array,meta.monkey,meta.lab,'ignoreJumps',meta.task,meta.mapfile);

%% Load marker file

marker_data = load([meta.folder 'ColorTracking\Markers\' 'markers_' meta.fname '.mat']);

%% Get TRC

% go to getTRCfromMarkers and run from there for now.
affine_xform = getTRCfromMarkers(cds,marker_data,[meta.folder 'OpenSim\']);

%% Do openSim stuff and save analysis results to analysis folder

% do this in opensim for now

%% Add kinematic information to CDS

% load joint information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'joint_ang')

% load joint velocities
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'joint_vel')

% load muscle information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'muscle_len')

% load muscle velocities
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'muscle_vel')

%% Save CDS

save([meta.folder 'CDS\' meta.fname '_CDS.mat'],'cds','-v7.3')