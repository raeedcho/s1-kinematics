function [behaviors_new] = extract_behaviors(behaviors,times)
% Get new behaviors struct with arm data subsampled from original behaviors
% between times(i,1) and times(i,2). Meant to be used after
% extract_workspace_times to extract a specific workspace of reaches.

