function [curves,bins] = getTuningCurves(moveCorr,FR)
% GETTUNINGCURVES Gets tuning curves for firing rate against
% moveCorr. Outputs curve in 8 bins, along with high and low CI.

dir = atan2(moveCorr(:,2),moveCorr(:,1));
spd = sqrt(sum(moveCorr.^2,2));

% bin directions
dir_bins = round(dir/(pi/4))*(pi/4);
dir_bins(dir_bins==-pi) = pi;
bins = -3*pi/4:pi/4:pi;
% bins = bins';

% find FR in each bin, along with CI
for i = 1:length(bins)
    % get FR when moveCorr is in the direction of bin
    % Also transpose FR so that rows are neurons and columns are observations
    FR_in_bin = FR(dir_bins==bins(i),:)';

    % Mean binned FR has normal-looking distribution (checked with
    % bootstrapping on a couple S1 neurons)
    binnedFR(:,i) = mean(FR_in_bin,2); % mean firing rate
    binned_stderr = std(FR_in_bin,0,2)/sqrt(size(FR_in_bin,2)); % standard error
    tscore = tinv(0.975,size(FR_in_bin,2)-1); % t-score for 95% CI
    binned_CIhigh(:,i) = binnedFR(:,i)+tscore*binned_stderr; %high CI
    binned_CIlow(:,i) = binnedFR(:,i)-tscore*binned_stderr; %low CI
end

% set up output struct
curves = table(binnedFR,binned_CIlow,binned_CIhigh,'VariableNames',{'binnedFR','CIlow','CIhigh'});
% curve.bins = {bins};
% curve.FR = {binnedFR};
% curve.CIhigh = {binned_CIhigh};
% curve.CIlow = {binned_CIlow};
