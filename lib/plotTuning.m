function [handle] = plotTuning(bins,modeledTuning,curve,maxRadius,color)
% PLOT_TUNING makes a single figure showing the tuning curve and PD with
% confidence intervals. Leave either entry blank to skip plotting it. Color
% is a 3 element vector for the color of the plotted tuning curve and PD.

% plot initial point
h=polar(0,maxRadius);
set(h,'color','w')
hold all

% tuning curve
if(~isempty(curve))
    h=polar(repmat(bins,1,2),repmat(curve.binnedFR,1,2));
    set(h,'linewidth',2,'color',color)
    th_fill = [fliplr(bins) bins(end) bins(end) bins];
    r_fill = [fliplr(curve.CIhigh) curve.CIhigh(end) curve.CIlow(end) curve.CIlow];
    [x_fill,y_fill] = pol2cart(th_fill,r_fill);
    patch(x_fill,y_fill,color,'facealpha',0.3,'edgealpha',0);
end

% PD
if(~isempty(modeledTuning))
    dir = atan2(modeledTuning(2),modeledTuning(1));
    h=polar(repmat(dir,2,1),maxRadius*[0;1]);
    set(h,'linewidth',2,'color',color)
%     th_fill = [modeledTuning.dir_CI(2) modeledTuning.dir modeledTuning.dir_CI(1) 0];
%     r_fill = [maxRadius maxRadius maxRadius 0];
%     [x_fill,y_fill] = pol2cart(th_fill,r_fill);
%     patch(x_fill,y_fill,color,'facealpha',0.3);
end