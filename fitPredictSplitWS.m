function fitPredictSplitWS(td)
% For multiworkspace files, with dl and pm workspaces:
%   * Fits three different coordinate frame models to data from both workspaces
%   * Predicts firing rates from models on test data in both workspaces
%   * Compare predicted firing rates to actual firing rates in test data

%% First, split td into pm and dl
[~,td_dl] = getTDidx(td,'epoch','dl','result','R');
[~,td_pm] = getTDidx(td,'epoch','pm','result','R');
minsize = min(length(td_dl),length(td_pm));
td_dl = td_dl(1:minsize);
td_pm = td_pm(1:minsize);

%% Get training and testing sets
num_folds = 10;
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
                        'in_signals',{'pos',[1 2];'vel',[1 2]},...
                        'out_signals',{'S1_spikes'});
[td_ext,glm_info_ext] = getModel(td_train,glm_ext_params)
