%% get bdfs
NEVNSx_DL.NS2 = [];
NEVNSx_PM.NS2 = [];
NEVNSx_PM.NS3 = [];
NEVNSx_PM.NS4 = [];
NEVNSx_PM.NS5 = [];
NEVNSx_DL.NS3 = [];
NEVNSx_DL.NS4 = [];
NEVNSx_DL.NS5 = [];
NEVNSx_DL.NEV = NEV_DL;
NEVNSx_PM.NEV = NEV_PM;

bdf_DL = get_nev_mat_data(NEVNSx_DL,6);
bdf_PM = get_nev_mat_data(NEVNSx_PM,6);

%% get units
bdf_DL = testAllTuning(bdf_DL);
bdf_PM = testAllTuning(bdf_PM);
[tuning_vec_DL,tuned_single_DL] = get_tuning_vec(bdf_DL);
[tuning_vec_PM,tuned_single_PM] = get_tuning_vec(bdf_PM);
tuned_single_both = intersect(tuned_single_DL, tuned_single_PM, 'rows');

%% get PDs
[pos_DL,vel_DL,acc_DL] = sample_bdf_kin(bdf_DL,50);
[pos_PM,vel_PM,acc_PM] = sample_bdf_kin(bdf_PM,50);
glm_input_DL = [pos_DL vel_DL sqrt(vel_DL(:,1).^2+vel_DL(:,2).^2)];
glm_input_PM = [pos_PM vel_PM sqrt(vel_PM(:,1).^2+vel_PM(:,2).^2)];

keep_idx_DL = 1:1.73e4;
keep_idx_PM = 1:1.5e4;

for ctr = 1:length(tuned_single_both)
    s_DL = bin_spikes(bdf_DL,50,tuned_single_both(ctr,1),tuned_single_both(ctr,2));
    s_PM = bin_spikes(bdf_PM,50,tuned_single_both(ctr,1),tuned_single_both(ctr,2));
    
    for repct = 1:100
        idx_rand_DL = uint32(1+(length(keep_idx_DL)-1)*rand(length(keep_idx_DL),1));
        b_bootstrap_DL = glmfit(glm_input_DL(idx_rand_DL,:),s_DL(idx_rand_DL),'poisson');
        PD_boot_DL = atan2d(b_bootstrap_DL(5),b_bootstrap_DL(4));
        
        idx_rand_PM = uint32(1+(length(keep_idx_PM)-1)*rand(length(keep_idx_PM),1));
        b_bootstrap_PM = glmfit(glm_input_PM(idx_rand_PM,:),s_PM(idx_rand_PM),'poisson');
        PD_boot_PM = atan2d(b_bootstrap_PM(5),b_bootstrap_PM(4));
        
        PD_diff_temp = PD_boot_PM-PD_boot_DL;
        PD_diff_temp(PD_diff_temp<-180) = PD_diff_temp(PD_diff_temp<-180)+360;
        PD_diff_temp(PD_diff_temp>180) = PD_diff_temp(PD_diff_temp>180)-360;
        
        PD_diff_boot(ctr,repct) = PD_diff_temp;
    end
    
    b_DL = glmfit(glm_input_DL(keep_idx_DL,:),s_DL(keep_idx_DL),'poisson');
    b_PM = glmfit(glm_input_PM(keep_idx_PM,:),s_PM(keep_idx_PM),'poisson');
    
    PD_DL(ctr) = atan2d(b_DL(5),b_DL(4));
    PD_PM(ctr) = atan2d(b_PM(5),b_PM(4));
    
    moddepth_DL(ctr) = sqrt(sum(b_DL(4:5).^2));
    moddepth_PM(ctr) = sqrt(sum(b_PM(4:5).^2));
end

PD_diff = PD_PM-PD_DL;
PD_diff(PD_diff<-180) = PD_diff(PD_diff<-180)+360;
PD_diff(PD_diff>180) = PD_diff(PD_diff>180)-360;

%% plot
max_moddepth = max([moddepth_DL moddepth_PM]);
for ctr = 1:14
%     figure(1)
%     subplot(4,4,ctr)
%     polar(0,-1,'.')
%     hold on
%     polar([PD_DL(ctr) PD_DL(ctr)]*pi/180,[0 moddepth_DL(ctr)/max_moddepth],'b')
%     polar([PD_DL(ctr)]*pi/180,[moddepth_DL(ctr)/max_moddepth],'bo')
%     polar([PD_PM(ctr) PD_PM(ctr)]*pi/180,[0 moddepth_PM(ctr)/max_moddepth],'r')
%     polar([PD_PM(ctr)]*pi/180,[moddepth_PM(ctr)/max_moddepth],'ro')
    
    figure(2)
    subplot(3,5,ctr)
    ksdensity(PD_diff_boot(ctr,:));
    axis([-90 90 0 0.2])
%     
%     figure(3)
%     subplot(4,4,ctr)
%     polar(0,-1,'.')
%     hold on
%     PD_diff_mean = mean(PD_diff_boot(ctr,:));
%     polar([PD_diff_mean PD_diff_mean]*pi/180,[0 (moddepth_PM(ctr)+moddepth_DL(ctr))/(2*max_moddepth)],'b')
end