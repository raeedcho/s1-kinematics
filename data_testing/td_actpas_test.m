%%
td_act_avg = trimTD(td_act,{'idx_movement_on',-15},{'idx_movement_on',30});
td_pas_avg = trimTD(td_pas,{'idx_movement_on',-15},{'idx_movement_on',30});

td_act_avg = trialAverage(td_act_avg,'tgtDir');
td_pas_avg = trialAverage(td_pas_avg,'bumpDir');

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
td_temp =  {td_act,td_pas};
sig_to_plot = 'speed';
figure
for trialnum = 1:30
    clf
    for i = 1:2
        subplot(2,2,i)
        plot(td_temp{i}(trialnum).(sig_to_plot),'k','linewidth',2)
        hold on
        plot(repmat(td_temp{i}(trialnum).idx_movement_on,2,1),ylim,'--b')
        plot(repmat(td_temp{i}(trialnum).idx_bumpTime,2,1),ylim,'--r')
        plot(repmat(td_temp{i}(trialnum).idx_goCueTime,2,1),ylim,'--g')
        plot(repmat(td_temp{i}(trialnum).idx_endTime,2,1),ylim,'--k')
        subplot(2,2,i+2)
        plot(td_temp{i}(trialnum).(strcat('d',sig_to_plot)),'k','linewidth',2)
        hold on
        plot(repmat(td_temp{i}(trialnum).idx_movement_on,2,1),ylim,'--b')
        plot(repmat(td_temp{i}(trialnum).idx_bumpTime,2,1),ylim,'--r')
        plot(repmat(td_temp{i}(trialnum).idx_goCueTime,2,1),ylim,'--g')
        plot(repmat(td_temp{i}(trialnum).idx_endTime,2,1),ylim,'--k')
        plot(xlim,zeros(2,1),'-k')
    end
    waitforbuttonpress
end

%%
[~,td_act_bin] = getTDidx(td_bin,'ctrHoldBump',false);
[~,td_pas_bin] = getTDidx(td_bin,'ctrHoldBump',true);

[dirs,~,dir_idx_act] = unique(cat(1,td_act_bin.tgtDir));
[~,~,dir_idx_pas] = unique(cat(1,td_pas_bin.bumpDir));
dir_colors = linspecer(max(dir_idx_act));
forcename = {'force',1:3};
kinname = {'marker_vel',7:9};

figure
subplot(1,2,1)
sig_act = getSig(td_act_bin,kinname);
sig_pas = getSig(td_pas_bin,kinname);
scatter3(sig_act(:,1),sig_act(:,2),sig_act(:,3),[],dir_colors(dir_idx_act,:),'filled')
hold on
scatter3(sig_pas(:,1),sig_pas(:,2),sig_pas(:,3),[],dir_colors(dir_idx_pas,:))
axis equal
subplot(1,2,2)
sig_act = getSig(td_act_bin,forcename);
sig_pas = getSig(td_pas_bin,forcename);
scatter3(sig_act(:,1),sig_act(:,2),sig_act(:,3),[],dir_colors(dir_idx_act,:),'filled')
hold on
scatter3(sig_pas(:,1),sig_pas(:,2),sig_pas(:,3),[],dir_colors(dir_idx_pas,:))
axis equal
suptitle(strrep(filenames{filenum},'_','-'))

%% try classifier on kinematics
num_sigs = 3;
[train_idx,test_idx] = crossvalind('HoldOut',length(td_bin),0.2);
mdl = fitcdiscr(getSig(td_bin(train_idx),{'marker_vel',1:num_sigs}),cat(1,td_bin(train_idx).ctrHoldBump));
sep = sum(predict(mdl,getSig(td_bin(test_idx),{'marker_vel',1:num_sigs})) == cat(1,td_bin(test_idx).ctrHoldBump))/sum(test_idx)

%%
[~,td_act_bin] = getTDidx(td_bin,'ctrHoldBump',false);
[~,td_pas_bin] = getTDidx(td_bin,'ctrHoldBump',true);

figure
scatter(getSig(td_act_bin,{'vel',1}),getSig(td_act_bin,{'marker_vel',7})*100,[],'k','filled')
hold on
scatter(getSig(td_pas_bin,{'vel',1}),getSig(td_pas_bin,{'marker_vel',7})*100,[],'k')
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

%%
    td = getNorm(td,struct('signals',{{'marker_vel',25:27}},'norm_name','elbow_speed'));
    td = getNorm(td,struct('signals',{{'marker_vel',7:9}},'norm_name','hand_speed'));
    figure('defaultaxesfontsize',18)
    for trial = 1:length(td)
        timevec = ((1:length(td(trial).elbow_speed))-td(trial).idx_movement_on)*td(trial).bin_size;
        if td(trial).ctrHoldBump
            plot(timevec,td(trial).elbow_speed,'r')
        else
            plot(timevec,td(trial).elbow_speed,'k')
        end
        hold on
    end
    plot(zeros(2,1),ylim,'--k','linewidth',2)
    hold on
    plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
    xlabel('Time from movement onset (s)')
    ylabel('Elbow speed (m/s)')
    set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
    set(gcf,'renderer','Painters')
    suptitle(strrep(strrep(filenames{filenum},'_COactpas_TD.mat',''),'_','-'))

