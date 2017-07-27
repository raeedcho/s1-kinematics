function plotPredCurves(td_test)
% Plots actual and predicted tuning curves for ext and musc models
% (assuming td_test has predicted firing rates in it already)

%% Separate test sets
[~,td_dl_test] = getTDidx(td_test,'epoch','DL');
[~,td_pm_test] = getTDidx(td_test,'epoch','PM');

%% Get tuning curves
vel_dl = cat(1,td_dl_test.vel);
spikes_dl = cat(1,td_dl_test.S1_spikes);
pred_ext_dl = cat(1,td_dl_test.glm_ext_model);
pred_musc_dl = cat(1,td_dl_test.glm_musc_model);

vel_pm = cat(1,td_pm_test.vel);
spikes_pm = cat(1,td_pm_test.S1_spikes);
pred_ext_pm = cat(1,td_pm_test.glm_ext_model);
pred_musc_pm = cat(1,td_pm_test.glm_musc_model);

[dl_curves,bins] = getTuningCurves(vel_dl,spikes_dl);
[dl_ext_curves] = getTuningCurves(vel_dl,pred_ext_dl);
[dl_musc_curves] = getTuningCurves(vel_dl,pred_musc_dl);

[pm_curves,bins] = getTuningCurves(vel_pm,spikes_pm);
[pm_ext_curves] = getTuningCurves(vel_pm,pred_ext_pm);
[pm_musc_curves] = getTuningCurves(vel_pm,pred_musc_pm);

%% Get PDs
params = struct('move_corr','vel','num_boot',10);

params.out_signals = 'S1_spikes';
dl_pds = getTDPDs(td_dl_test,params);
pm_pds = getTDPDs(td_pm_test,params);

params.out_signals = 'glm_ext_model';
dl_ext_pds = getTDPDs(td_dl_test,params);
pm_ext_pds = getTDPDs(td_pm_test,params);

params.out_signals = 'glm_musc_model';
dl_musc_pds = getTDPDs(td_dl_test,params);
pm_musc_pds = getTDPDs(td_pm_test,params);

%% Plot tuning curves
figure
for neuron_idx = 1:height(dl_curves)
    clf;
    maxFR = max([dl_curves(neuron_idx,:).CIhigh, pm_curves(neuron_idx,:).CIhigh, ...
                    dl_ext_curves(neuron_idx,:).CIhigh, pm_ext_curves(neuron_idx,:).CIhigh, ...
                    dl_musc_curves(neuron_idx,:).CIhigh, pm_musc_curves(neuron_idx,:).CIhigh]);
    subplot(2,2,1)
    plotTuning(bins,dl_pds(neuron_idx,:),dl_curves(neuron_idx,:),maxFR,[1 0 0]);
    plotTuning(bins,pm_pds(neuron_idx,:),pm_curves(neuron_idx,:),maxFR,[0.6 0.5 0.7]);
    title 'actual curves'
    
    subplot(2,2,2)
    plotTuning(bins,dl_ext_pds(neuron_idx,:),dl_ext_curves(neuron_idx,:),maxFR,[1 0 0]);
    plotTuning(bins,pm_ext_pds(neuron_idx,:),pm_ext_curves(neuron_idx,:),maxFR,[0.6 0.5 0.7]);
    title 'ext curves'
    
    subplot(2,2,3)
    plotTuning(bins,dl_musc_pds(neuron_idx,:),dl_musc_curves(neuron_idx,:),maxFR,[1 0 0]);
    plotTuning(bins,pm_musc_pds(neuron_idx,:),pm_musc_curves(neuron_idx,:),maxFR,[0.6 0.5 0.7]);
    title 'musc curves'

    subplot(2,2,4)
    title(['Neuron ' num2str(neuron_idx)])

    waitforbuttonpress;
end

%% Plot PD shift summary
CIwidth_dl = diff(dl_pds.velDirCI-repmat(dl_pds.velDir,1,2),1,2);
CIwidth_pm = diff(pm_pds.velDirCI-repmat(pm_pds.velDir,1,2),1,2);
tuned_idx = CIwidth_dl<pi/4 & CIwidth_pm<pi/4;
pdDiff = remove_wrapping(dl_pds.velDir-pm_pds.velDir);
pdDiff_ext = remove_wrapping(dl_ext_pds.velDir-pm_ext_pds.velDir);
pdDiff_musc = remove_wrapping(dl_musc_pds.velDir-pm_musc_pds.velDir);

pdDiff = pdDiff(tuned_idx);
pdDiff_ext = pdDiff_ext(tuned_idx);
pdDiff_musc = pdDiff_musc(tuned_idx);

% find linear models of pdDiff
ext_diff_lm = fitlm(pdDiff,pdDiff_ext)
musc_diff_lm = fitlm(pdDiff,pdDiff_musc)
ext_pred_diff = ext_diff_lm.predict((-pi:pi/2:pi)');
musc_pred_diff = musc_diff_lm.predict((-pi:pi/2:pi)');

figure
ext_color = [1 0.5 0];
musc_color = [0.12 0.69 0.67];
plot([-pi pi],[0 0],'-k','linewidth',2)
hold on
plot([0 0],[-pi pi],'-k','linewidth',2)
plot([-pi pi],[-pi pi],'--k','linewidth',2)
scatter(pdDiff,pdDiff_ext,50,ext_color,'x','linewidth',2)
scatter(pdDiff,pdDiff_musc,50,musc_color,'o','linewidth',2)
plot(-pi:pi/2:pi,ext_pred_diff,'--','Color',ext_color,'linewidth',2)
plot(-pi:pi/2:pi,musc_pred_diff,'--','Color',musc_color,'linewidth',2)
set(gca,'xlim',[-pi pi],'ylim',[-pi pi],'box','off','tickdir','out')
axis square

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

