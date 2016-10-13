function [handle] = plot_tuning(modeled_tuning,curve,max_FR,color)
% PLOT_TUNING makes a single figure showing the tuning curve and PD with
% confidence intervals. Leave either entry blank to skip plotting it. Color
% is a 3 element vector for the color of the plotted tuning curve and PD.

% plot initial point
h=polar(0,max_FR);
set(h,'color','w')
hold all

% tuning curve
if(~isempty(curve))
    h=polar(repmat(curve.bins{1},2,1),repmat(curve.FR{1},2,1));
    set(h,'linewidth',2,'color',color)
    th_fill = [flipud(curve.bins{1}); curve.bins{1}(end); curve.bins{1}(end); curve.bins{1}];
    r_fill = [flipud(curve.CI_high{1}); curve.CI_high{1}(end); curve.CI_low{1}(end); curve.CI_low{1}];
    [x_fill,y_fill] = pol2cart(th_fill,r_fill);
    patch(x_fill,y_fill,color,'facealpha',0.3,'edgealpha',0);
end

% PD
if(~isempty(modeled_tuning))
    h=polar(repmat(modeled_tuning.dir,2,1),max_FR*[0;1]);
    set(h,'linewidth',2,'color',color)
    th_fill = [modeled_tuning.dir_CI(2) modeled_tuning.dir modeled_tuning.dir_CI(1) 0];
    r_fill = [max_FR max_FR max_FR 0];
    [x_fill,y_fill] = pol2cart(th_fill,r_fill);
    patch(x_fill,y_fill,color,'facealpha',0.3);
end