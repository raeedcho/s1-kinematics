function iris_plot(tuning_PM,tuning_DL)
% make iris plot from tuning structure given by calc_PD_helper

% extract relevant information
angs_PM = [tuning_PM.dir];
dir_CI_PM = [tuning_PM.dir_CI];
angs_DL = [tuning_DL.dir];
dir_CI_DL = [tuning_DL.dir_CI];

% calculate CI widths
DL_CI_width = diff(dir_CI_DL,1,2); % get CI widths
PM_CI_width = diff(dir_CI_PM,1,2);
DL_CI_width(DL_CI_width<0) = DL_CI_width(DL_CI_width<0)+2*pi;
PM_CI_width(PM_CI_width<0) = PM_CI_width(PM_CI_width<0)+2*pi;
tuned_neurons = DL_CI_width<pi/4 & PM_CI_width<pi/4;
angs_PM_tuned = angs_PM(tuned_neurons);
angs_DL_tuned = angs_DL(tuned_neurons);

%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

% plot changes with alpha dependent on CI width
for unit_ctr = 1:length(angs_PM_tuned)
    h=polar(linspace(angs_PM_tuned(unit_ctr),angs_DL_tuned(unit_ctr),2),linspace(0.5,1,2));
    set(h,'linewidth',2,'color',[0.1 0.6 1])
end

%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

set(findall(gcf, 'String','  0.2','-or','String','  0.4','-or','String','  0.6','-or','String','  0.8',...
        '-or','String','  1') ,'String', ' '); % remove a bunch of labels from the polar plot; radial and tangential