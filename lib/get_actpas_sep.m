function [separability,mdl] = get_actpas_sep(td,params)
% GET_ACTPAS_SEP uses LDA to classify active and passive trials from signals
% Inputs:
%   td - trial_data structure
%   params - params struct
%       use_trials - use these trials to calculate model
%       signals - signals to use for LDA
%       do_plot - whether to plot separability plot
%       mdl - if not empty, will use this model instead of fitting
%           a new one
% Outputs:
%   separability - classification accuracy of LDA on signals
%   mdl - output model from LDA fitting
%
% Note: does not cross-validate by itself. If you want crossval, fit
% a model with training data, then pass in test data along with the
% model into params.mdl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETER VALUES
use_trials      =  1:length(td);
signals         =  getTDfields(td,'spikes');
do_plot         =  false;
fig_handle      = [];
mdl             =  [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra parameters you can change that aren't described in header
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process and prepare inputs
signals = check_signals(td(1),signals);
if iscell(use_trials) % likely to be meta info
    use_trials = getTDidx(td,use_trials{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals = check_signals(td,signals);

% ideally, this would work like trialAverage, where the function would take
% in a condition
% But...I'm not sure how to deal with more than two values for the condition
[~,td_act] = getTDidx(td,'ctrHoldBump',false);
[~,td_pas] = getTDidx(td,'ctrHoldBump',true);

signal_act = get_vars(td_act,signals);
signal_pas = get_vars(td_pas,signals);

% Find total separability
signal = cat(1,signal_act,signal_pas);
actpas = [ones(length(signal_act),1);zeros(length(signal_pas),1)];

% get model
if isempty(mdl)
    mdl = fitcdiscr(signal,actpas);
end
class = predict(mdl,signal);
separability = sum(class == actpas)/length(actpas);

if do_plot
    if ~isempty(fig_handle)
        figure(fig_handle)
        hold on
    else
        % figure
    end
    % plot active as filled, passive as open
    % bump_colors = linspecer(4);
    % act_dir_idx = floor(cat(1,td_act.target_direction)/(pi/2))+1;
    % pas_dir_idx = floor(cat(1,td_pas.bumpDir)/90)+1;
    act_color = [114 191 111]/256;
    pas_color = [88 137 176]/256;
    
    w = mdl.Sigma\diff(mdl.Mu)';
    signal_sep = signal*w;

    % get basis vector orthogonal to w for plotting
    null_sep = null(w');
    signal_null_sep = signal*null_sep;
    [~,signal_null_sep_scores] = pca(signal_null_sep);

    % plot for act/pas separability
    hold all
    scatter3(signal_sep(actpas==1),signal_null_sep_scores(actpas==1,1),signal_null_sep_scores(actpas==1,2),100,act_color,'filled')
    scatter3(signal_sep(actpas==0),signal_null_sep_scores(actpas==0,1),signal_null_sep_scores(actpas==0,2),100,pas_color,'filled')
    ylim = get(gca,'ylim');
    zlim = get(gca,'zlim');
    plot3([0 0],ylim,[0 0],'--k','linewidth',2)
    plot3([0 0],[0 0],zlim,'--k','linewidth',2)
    set(gca,'box','off','tickdir','out')
    view([0 0])
    axis equal
    axis off
end
