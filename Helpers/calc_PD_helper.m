function [out_tuning] = calc_PD_helper(bootfunc,endpoint_kin,FR)
% Helper function for plot_PD_predictions and iris plots

    boot_tuning = bootstrp(100,@(X,y) {bootfunc(X,y)}, endpoint_kin, FR);
    
    %extract coefficiencts from boot_tuning
    boot_coef = cell2mat(cellfun(@(x) x.Coefficients.Estimate',boot_tuning,'uniformoutput',false));
    
    %coefficient covariance
%     coef_cov = cov(boot_coef);
    
    %get coefficient means
    coef_means = mean(boot_coef);
    
    %get 95% CIs for coefficients
    coef_CIs = prctile(boot_coef,[2.5 97.5]); 
    
    % bootstrap directions
    boot_dirs = atan2(boot_coef(:,5),boot_coef(:,4));
    % recenter boot_dirs
    mean_dir = atan2(coef_means(5),coef_means(4));
    centered_boot_dirs = boot_dirs-mean_dir;
    while(sum(centered_boot_dirs<-pi))
        centered_boot_dirs(centered_boot_dirs<-pi) = centered_boot_dirs(centered_boot_dirs<-pi)+2*pi;
    end
    while(sum(centered_boot_dirs>pi))
        centered_boot_dirs(centered_boot_dirs>pi) = centered_boot_dirs(centered_boot_dirs>pi)-2*pi;
    end

    % Calculate dir CI
    dir_CI = prctile(centered_boot_dirs,[2.5 97.5]);

    % uncenter CI
    dir_CI = dir_CI+mean_dir;
    
    out_tuning.dir = mean_dir;
    out_tuning.dir_CI = dir_CI;
    out_tuning.moddepth = sqrt(sum(coef_means(4:5).^2));
end