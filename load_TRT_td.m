%%

meta.lab=6;
meta.ranBy='ranByRaeed';
meta.monkey='monkeyChips';
meta.task='taskTRT';
meta.array='arrayLeftS1Area2';
meta.folder='C:\Users\rhc307\Projects\limblab\data-preproc\MultiWorkspace\TRT\Chips\20170906\';
meta.fname='Chips_20170906_TRT_area2_001';
meta.mapfile='mapFileC:\Users\rhc307\Projects\limblab\data-preproc\Meta\Mapfiles\Chips\left_S1\SN 6251-001455.cmp';

% Make CDS files
cds = commonDataStructure();
cds.file2cds([meta.folder 'preCDS\' meta.fname],meta.ranBy,meta.array,meta.monkey,meta.lab,'ignoreJumps',meta.task,meta.mapfile);

%%

params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = [0,255];
params.event_list = {'spaceNum';'targetStartTime'};
params.trial_results = {'R','A','F','I'};
meta = struct('task','TRT');
params.meta = meta;

trial_data = parseFileByTrial(cds,params);