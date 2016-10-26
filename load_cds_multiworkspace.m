%% Load multiworkspace data into CDS files
addpath(genpath('C:\Users\Raeed\Projects\limblab\ClassyDataAnalysis'))

lab=6;
ranBy='ranByRaeed';
monkey='monkeyHan';
task='taskRW';
array='arrayLeftS1Area2';
% folder='C:\Users\Raeed\Projects\limblab\data-raeed\MultiWorkspace\SplitWS\Han\20160322\area2\preCDS\';
folder = '/home/raeed/Projects/limblab/data-raeed/MultiWorkspace/SplitWS/Han/20160322/area2/preCDS/';
DL_fname='Han_20160322_RW_DL_area2_001';
PM_fname='Han_20160322_RW_PM_area2_002';
mapfile='mapFileC:\Users\Raeed\Projects\limblab\data-raeed\MultiWorkspace\SplitWS\Han\mapfile\left_S1\SN 6251-001459.cmp';

% Make CDS files

DL_cds = commonDataStructure();
DL_cds.file2cds([folder DL_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

PM_cds = commonDataStructure();
PM_cds.file2cds([folder PM_fname],ranBy,array,monkey,lab,'ignoreJumps',task,mapfile);

%% Save CDS files
% savepath = 'C:\Users\Raeed\Projects\limblab\data-raeed\MultiWorkspace\SplitWS\Han\20160322\area2\CDS\';
savepath = '/home/raeed/Projects/limblab/data-raeed/MultiWorkspace/SplitWS/Han/20160322/area2/CDS/';
save([savepath DL_fname '_cds.mat'],'DL_cds')
save([savepath PM_fname '_cds.mat'],'PM_cds')