%% Load data into CDS file

meta.lab=6;
meta.ranBy='ranByRaeed';
meta.monkey='monkeyChips';
meta.task='taskRW';
meta.array='arrayLeftS1Area2';
meta.folder='C:\Users\rhc307\Projects\limblab\data-preproc\MultiWorkspace\FullWS\Chips\leftS1\20151203\';
meta.fname='Chips_20151203_RW_002';
meta.mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Chips\left_S1\SN 6251-001455.cmp';

% Make CDS files

cds = commonDataStructure();
cds.file2cds([meta.folder 'preCDS\' meta.fname],meta.ranBy,meta.array,meta.monkey,meta.lab,'ignoreJumps',meta.task,meta.mapfile);

%% Load marker file

marker_data = load([meta.folder 'ColorTracking\Markers\' 'markers_' meta.fname '.mat']);

%% Get TRC

% go to getTRCfromMarkers and run from there for now.
getTRCfromMarker(cds,marker_dta)

%% Do openSim stuff and save analysis results to analysis folder

% do this in opensim for now

%% Add kinematic information to CDS

% load joint information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'joint_ang')

% load muscle information
cds.loadOpenSimData([meta.folder 'OpenSim\Analysis\'],'muscle_len')

%% Save CDS

save([meta.folder 'CDS\' meta.fname '_CDS.mat'],'cds','-v7.3')