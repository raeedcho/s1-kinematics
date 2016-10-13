%% Set up data

lab=6;
ranBy='ranByRaeed';
monkey='monkeyChips';
task='taskRW';
array='arrayLeftS1Area2';
%note the .nev extension is not necessary when providing the file name:
folder='/home/raeed/Projects/limblab/data-raeed/MultiWorkspace/SplitWS/Chips/20151202/';
DL_fname='Chips_20151202_RW_DL_001';
PM_fname='Chips_20151202_RW_PM_003';
mapfile='mapFile/home/raeed/Projects/limblab/data-raeed/MultiWorkspace/SplitWS/Chips/mapfile/left_S1/SN 6251-001455.cmp';

%% Make CDS files

DL_cds = commonDataStructure();
DL_cds.file2cds([folder DL_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

PM_cds = commonDataStructure();
PM_cds.file2cds([folder PM_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

%% pass data to experiments

% DL stuff
DL_ex = experiment();

% set variables to load from cds
DL_ex.meta.hasLfp=false;
DL_ex.meta.hasKinematics=true;
DL_ex.meta.hasForce=false;
DL_ex.meta.hasUnits=true;
DL_ex.meta.hasTrials=true;
DL_ex.meta.hasAnalog=true;

% add session to experiment
DL_ex.addSession(DL_cds);

% PM stuff
PM_ex = experiment();

% set variables to load from cds
PM_ex.meta.hasLfp=false;
PM_ex.meta.hasKinematics=true;
PM_ex.meta.hasForce=false;
PM_ex.meta.hasUnits=true;
PM_ex.meta.hasTrials=true;
PM_ex.meta.hasAnalog=true;

% add session to experiment
PM_ex.addSession(DL_cds);

%% Bin experiment data

DL_ex.firingRateConfig.cropType='tightCrop';
DL_ex.firingRateConfig.offset=-0.015;

DL_ex.binData()

PM_ex.firingRateConfig.cropType='tightCrop';
PM_ex.firingRateConfig.offset=-0.015;

PM_ex.binData()