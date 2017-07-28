%% load file
load('/home/raeed/Projects/limblab/data-td/MultiWorkspace/FullWS/Han/20160325/Han_20160325_FullWS_area2_TD-notracking.mat')

%% Split into train and test
td = binTD(trial_data,5);
[train_idx,test_idx] = crossvalind('HoldOut',length(td),0.2);

td_train = td(train_idx);

glm_ext_params = struct('model_type','glm',...
                        'model_name','ext_model',...
                        'in_signals',{{'pos',[1 2];'vel',[1 2]}},...
                        'out_signals',{{'S1_spikes'}});
[~,glm_info_ext] = getModel(td_train,glm_ext_params);

td_test = td(test_idx);
td_test = getModel(td_test,glm_info_ext);

%% evaluate model
eval_params = glm_info_ext;
eval_params.eval_metric = 'pr2';
td_ext_eval = squeeze(evalModel(td_test,eval_params))

%% plot evals
meaneval = mean(td_ext_eval,2);
err = diff(td_ext_eval,1,2);
neur_num = (1:length(err))';
errorbar(neur_num,meaneval,err,'o','linewidth',2)
hold on
plot(neur_num,zeros(size(neur_num)),'-k','linewidth',2)
title('Pseudo R^2 for neurons using position+velocity GLM')
xlabel('Neuron number')
ylabel('Pseudo R^2')
set(gca,'box','off','tickdir','out','xlim',[0 neur_num(end)+1])
