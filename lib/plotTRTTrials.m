%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotTRTTrials(trial_data,params)
%   For a TRT trial_data struct, plot the target centers and the handle
%   position throughout each given trial on one figure. Meant to be used
%   for only a few trials at a time, or trials in different parts of the
%   workspace
% 
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .targetcolor    :   plotting color for targets, used in scatter()
%                           default: [1 0 0]
%       .targetsize     :   size of targets, used in scatter()
%                           default: 50
%       .linespec       :   plotting spec for handle position, used in plot()
%                           default: 'k-'
%       .linewidth      :   linewidth for handle position, used in plot()
%                           default: 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTRTTrials(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
targetcolor = [1 0 0];
targetsize = 100;
linespec = 'k-';
linewidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>1
    assignParams(who,params); % overwrite parameters
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get hold state
held = ishold;

% Loop through trials
for trial_idx = 1:length(trial_data)
    % first plot targets with scatter
    % (these targets are wrong because of an error in behavior calculation)
    targets = trial_data(trial_idx).target_center-repmat([0 35],4,1);
    scatter(targets(:,1),targets(:,2),targetsize,targetcolor,'filled')

    hold on

    % then plot handle position with plot
    pos = trial_data(trial_idx).pos;
    plot(pos(:,1),pos(:,2),linespec,'linewidth',linewidth)
end
axis equal

% restore hold state
if ~held
    hold off
end
