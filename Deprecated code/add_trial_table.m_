function bdf_new = add_trial_table(bdf_old)
% Adds random walk trial table
% 
% Each row of the table coresponds to a single trial.  Columns are as
% follows:
%    1: Start time
%    2: max number of targets
%    3: number of targets attempted
%    4: x offset
%    5: y offset
%    6: target size - tolerance
%    [7->6+(2*num_tgts)]: [Go_types, Go_times] (Go_types: 0=Center_target_on, 1=Go_cue, 2=Catch_trial)
%    (6+2*num_tgts)+1   : Trial End time
%    (6+2*num_tgts)+2   : Trial result    -- R, A, F, I or N (N coresponds to no-result)

tt = rw_trial_table(bdf_old);
bdf_new = bdf_old;
bdf_new.trial_table = tt;

end