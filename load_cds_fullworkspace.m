%% Load data into CDS file

meta.lab=6;
meta.ranBy='ranByRaeed';
meta.monkey='monkeyHan';
meta.task='taskRW';
meta.array='arrayLeftS1Area2';
meta.folder='C:\Users\rhc307\Projects\limblab\data-preproc\MultiWorkspace\FullWS\Han\leftS1\20160315\';
meta.fname='Han_20160315_RW_area2_001';
meta.mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Han\left_S1\SN 6251-001459.cmp';

% Make CDS files

cds = commonDataStructure();
cds.file2cds([meta.folder 'preCDS\' meta.fname],meta.ranBy,meta.array,meta.monkey,meta.lab,'ignoreJumps',meta.task,meta.mapfile);

%% Load marker file

marker_data = load([meta.folder 'ColorTracking\Markers\' 'markers_Han_20160315_RW_001' '.mat']);

%% Get TRC

% go to getTRCfromMarkers and run from there for now.

%% Do openSim stuff and save analysis results to analysis folder

% do this in opensim for now

%% Add kinematic information to CDS

% load joint information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'_Kinematics_q.sto')

% load muscle information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'_MuscleAnalysis_Length.sto')

%% Save CDS

save([meta.folder 'CDS\' meta.fname '_CDS.mat'],'cds','-v7.3')