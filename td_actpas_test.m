%%
td_act_avg = trimTD(td_act,{'idx_movement_on',-15},{'idx_movement_on',30});
td_pas_avg = trimTD(td_pas,{'idx_movement_on',-15},{'idx_movement_on',30});

td_act_avg = trialAverage(td_act_avg,'tgtDir');
td_pas_avg = trialAverage(td_pas_avg,'bumpDir');

%%
num_vars = 1;
dir_colors = linspecer(4);
figure
for i=1:4
    for varnum = 1:num_vars
        subplot(num_vars,1,varnum)
        plot(td_act_avg(i).speed(:,varnum),'-','linewidth',2,'color',dir_colors(i,:))
        hold on
        plot(td_pas_avg(i).speed(:,varnum),'--','linewidth',2,'color',dir_colors(i,:))
    end
end

%%
td_temp = td_act;
sig_to_plot = 'plane_force_norm';
figure
for trialnum = 1:30
    clf
    subplot(2,1,1)
    plot(td_temp(trialnum).(sig_to_plot),'k','linewidth',2)
    hold on
    plot(repmat(td_temp(trialnum).idx_movement_on,2,1),ylim,'--b')
    plot(repmat(td_temp(trialnum).idx_bumpTime,2,1),ylim,'--r')
    plot(repmat(td_temp(trialnum).idx_goCueTime,2,1),ylim,'--g')
    plot(repmat(td_temp(trialnum).idx_endTime,2,1),ylim,'--k')
    subplot(2,1,2)
    plot(td_temp(trialnum).(strcat('d',sig_to_plot)),'k','linewidth',2)
    hold on
    plot(repmat(td_temp(trialnum).idx_movement_on,2,1),ylim,'--b')
    plot(repmat(td_temp(trialnum).idx_bumpTime,2,1),ylim,'--r')
    plot(repmat(td_temp(trialnum).idx_goCueTime,2,1),ylim,'--g')
    plot(repmat(td_temp(trialnum).idx_endTime,2,1),ylim,'--k')
    plot(xlim,zeros(2,1),'-k')
    waitforbuttonpress
end

%%
[~,td_act_bin] = getTDidx(td_bin,'ctrHoldBump',false);
[~,td_pas_bin] = getTDidx(td_bin,'ctrHoldBump',true);

signame = 'marker_vel';
figure
scatter3(getSig(td_act_bin,{signame,1}),getSig(td_act_bin,{signame,2}),getSig(td_act_bin,{signame,3}),[],'k','filled')
hold on
scatter3(getSig(td_pas_bin,{signame,1}),getSig(td_pas_bin,{signame,2}),getSig(td_pas_bin,{signame,3}),[],'k')
axis equal

%% try classifier on kinematics
num_sigs = 3;
[train_idx,test_idx] = crossvalind('HoldOut',length(td_bin),0.2);
mdl = fitcdiscr(getSig(td_bin(train_idx),{'marker_vel',1:num_sigs}),cat(1,td_bin(train_idx).ctrHoldBump));
sep = sum(predict(mdl,getSig(td_bin(test_idx),{'marker_vel',1:num_sigs})) == cat(1,td_bin(test_idx).ctrHoldBump))/sum(test_idx)

%%
[~,td_act_bin] = getTDidx(td_bin,'ctrHoldBump',false);
[~,td_pas_bin] = getTDidx(td_bin,'ctrHoldBump',true);

figure
scatter(getSig(td_act_bin,{'vel',2}),getSig(td_act_bin,{'marker_vel',2})*100,[],'k','filled')
hold on
scatter(getSig(td_pas_bin,{'vel',2}),getSig(td_pas_bin,{'marker_vel',2})*100,[],'k')
plot(xlim,xlim,'--k')
axis equal

%%
figure
ax = zeros(2,1);
td_temp = td_pas;
for signum = 1:2
    ax(signum) = subplot(2,1,signum);
    plot(getSig(td_temp,{'vel',signum}),'-k')
    hold on
    plot(getSig(td_temp,{'marker_vel',signum})*100,'-r')
end
linkaxes(ax,'x')