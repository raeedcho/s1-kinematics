function [curve] = get_single_tuning_curve(move_corr,FR)
% GET_SINGLE_TUNING_CURVE Gets single tuning curve for firing rate against
% move_corr. Outputs curve in 8 bins, along with high and low CI.

% check FR
if(size(FR,2)>1)
    error('This function is only for single tuning curves')
end

dir = atan2(move_corr(:,2),move_corr(:,1));
spd = sqrt(sum(move_corr.^2,2));

% bin directions
dir_bins = round(dir/(pi/4))*(pi/4);
dir_bins(dir_bins==-pi) = pi;
bins = -3*pi/4:pi/4:pi;
bins = bins';

% find FR in each bin, along with CI
for i = 1:length(bins)
    FR_in_bin = FR(dir_bins==bins(i));
    spd_in_bin = spd(dir_bins==bins(i));
    
    % Mean binned FR has normal-looking distribution (checked with
    % bootstrapping)
    binned_FR(i,:) = mean(FR_in_bin); % mean firing rate
    binned_stderr = std(FR_in_bin)/sqrt(length(FR_in_bin)); % standard error
    tscore = tinv(0.975,length(FR_in_bin)-1); % t-score for 95% CI
    binned_CI_high(i,:) = binned_FR(i,:)+tscore*binned_stderr; %high CI
    binned_CI_low(i,:) = binned_FR(i,:)-tscore*binned_stderr; %low CI
end

curve.bins = {bins};
curve.FR = {binned_FR};
curve.CI_high = {binned_CI_high};
curve.CI_low = {binned_CI_low};