% function fitPredictSplitWS(td)
% For multiworkspace files, with dl and pm workspaces:
%   * Fits three different coordinate frame models to data from both workspaces
%   * Predicts firing rates from models on test data in both workspaces
%   * Compare predicted firing rates to actual firing rates in test data

%% Prep inputs by binning and running PCA on muscles
td = binTD(trial_data,5);

% do PCA on muscles
PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
                    'do_plot',true);
[td,pca_info] = getPCA(td,PCAparams);

% temporary hack to allow us to do PCA on velocity too
for i=1:length(td)
    td(i).opensim_len_pca = td(i).opensim_pca;
end

% get velocity PCA
PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
                    'do_plot',true);
[td,pca_info_vel] = getPCA(td,PCAparams_vel);

% temporary hack to allow us to save into something useful
for i=1:length(td)
    td(i).opensim_muscVel_pca = td(i).opensim_pca;
end

% get rid of superfluous PCA
td = rmfield(td,'opensim_pca');

%% split into epochs
[~,td_dl] = getTDidx(td,'epoch','DL','result','R');
[~,td_pm] = getTDidx(td,'epoch','PM','result','R');
minsize = min(length(td_dl),length(td_pm));
td_dl = td_dl(1:minsize);
td_pm = td_pm(1:minsize);

%% Get training and testing sets
% num_folds = 10;
[train_idx,test_idx] = crossvalind('HoldOut',minsize,0.2);

%% Plot handle positions
% figure
% pos_dl = cat(1,td_dl.pos);
% plot(pos_dl(:,1),pos_dl(:,2),'r')
% mean(pos_dl)
% hold on
% pos_pm = cat(1,td_pm.pos);
% plot(pos_pm(:,1),pos_pm(:,2),'b')
% axis equal

%% Fit models
% set up model fits
td_train = [td_dl(train_idx) td_pm(train_idx)];

glm_ext_params = struct('model_type','glm',...
                        'model_name','ext_model',...
                        'in_signals',{{'pos',[1 2]; 'vel',[1 2]}},...
                        'out_signals',{{'S1_spikes'}});
glm_musc_params = struct('model_type','glm',...
                        'model_name','musc_model',...
                        'in_signals',{{'opensim_len_pca',1:5;'opensim_muscVel_pca',1:5}},...
                        'out_signals',{'S1_spikes'});
[~,glm_info_ext] = getModel(td_train,glm_ext_params);
[~,glm_info_musc] = getModel(td_train,glm_musc_params);

%% Predict firing rates and evaluate models
td_dl_test = td_dl(test_idx);
td_pm_test = td_pm(test_idx);
td_dl_test = getModel(td_dl_test,glm_info_ext);
td_pm_test = getModel(td_pm_test,glm_info_ext);
td_dl_test = getModel(td_dl_test,glm_info_musc);
td_pm_test = getModel(td_pm_test,glm_info_musc);

% evaluate models
eval_params = glm_info_ext;
eval_params.eval_metric = 'pr2';
td_ext_dl_eval = squeeze(evalModel(td_dl_test,eval_params));
td_ext_pm_eval = squeeze(evalModel(td_pm_test,eval_params));

eval_params = glm_info_musc;
eval_params.eval_metric = 'pr2';
td_musc_dl_eval = squeeze(evalModel(td_dl_test,eval_params));
td_musc_pm_eval = squeeze(evalModel(td_pm_test,eval_params));

% package test sets together
td_test = cat(2,td_dl_test,td_pm_test)

%% Plot example
% neuron 1 has decent pseudo R2
neuron_idx = 1;

temp_vel = cat(1,td_dl_test.vel);
temp_spikes = cat(1,td_dl_test.S1_spikes);
temp_pred_ext = cat(1,td_dl_test.glm_ext_model);
temp_pred_musc = cat(1,td_dl_test.glm_musc_model);

figure
ax1 = subplot(2,1,1);
plot(temp_vel(:,1),'b','linewidth',2)
hold on
plot(temp_vel(:,2),'g','linewidth',2)
set(gca,'box','off','tickdir','out')

ax2 = subplot(2,1,2);
plot(temp_spikes(:,neuron_idx),'k','linewidth',2)
hold on
plot(temp_pred_ext(:,neuron_idx),'r','linewidth',2)
plot(temp_pred_musc(:,neuron_idx),'c','linewidth',2)
set(gca,'box','off','tickdir','out')

linkaxes([ax1 ax2],'x')

%% Plot pR2s against each other
av_pR2_ext_dl = mean(td_ext_dl_eval,2);
av_pR2_musc_dl = mean(td_musc_dl_eval,2);
av_pR2_ext_pm = mean(td_ext_pm_eval,2);
av_pR2_musc_pm = mean(td_musc_pm_eval,2);

good_neurons = td_ext_dl_eval(:,1) > 0 & td_ext_pm_eval(:,1) > 0 & td_musc_pm_eval(:,1) > 0 & td_musc_dl_eval(:,1) > 0;

figure
plot(av_pR2_ext_dl(good_neurons),av_pR2_musc_dl(good_neurons),'ro','linewidth',2)
hold on
plot([-1 1],[-1 1],'k--','linewidth',2)
plot([0 0],[-1 1],'k-','linewidth',2)
plot([-1 1],[0 0],'k-','linewidth',2)
set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])

figure
plot(av_pR2_ext_pm(good_neurons),av_pR2_musc_pm(good_neurons),'bo','linewidth',2)
hold on
plot([-1 1],[-1 1],'k--','linewidth',2)
plot([0 0],[-1 1],'k-','linewidth',2)
plot([-1 1],[0 0],'k-','linewidth',2)
set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])

