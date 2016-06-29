%% Plot example PD change diagram assuming extrinsic
h=figure('name','extrinsic_PD_diff');
figure_handles=[figure_handles h];
%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

% plot changes with alpha dependent on CI width
angs_ext = linspace(0,2*pi,16);
for unit_ctr = 1:15
    h=polar(linspace(angs_ext(unit_ctr),angs_ext(unit_ctr),2),linspace(0.5,1,2));
    set(h,'linewidth',2,'color',[0.1 0.6 1])
end
%plot circles again
% h=polar(linspace(-pi,pi,1000),ones(1,1000));
% set(h,'linewidth',2,'color',[1 0 0])
% hold all
% h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
% set(h,'linewidth',2,'color',[0.6 0.5 0.7])

set(findall(gcf, 'String','  0.2','-or','String','  0.4','-or','String','  0.6','-or','String','  0.8',...
        '-or','String','  1') ,'String', ' '); % remove a bunch of labels from the polar plot; radial and tangential

title('Plot of PD changes (extrinsic)')

%% Plot example PD change diagram assuming sort of intrinsic
h=figure('name','intrinsic_PD_diff');
figure_handles=[figure_handles h];
%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

% plot changes with alpha dependent on CI width
angs_ext = linspace(0,2*pi,16);
for unit_ctr = 1:15
    h=polar(linspace(angs_ext(unit_ctr),angs_ext(unit_ctr)-pi/6,2),linspace(0.5,1,2));
    set(h,'linewidth',2,'color',[0.1 0.6 1])
end
%plot circles again
% h=polar(linspace(-pi,pi,1000),ones(1,1000));
% set(h,'linewidth',2,'color',[1 0 0])
% hold all
% h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
% set(h,'linewidth',2,'color',[0.6 0.5 0.7])

set(findall(gcf, 'String','  0.2','-or','String','  0.4','-or','String','  0.6','-or','String','  0.8',...
        '-or','String','  1') ,'String', ' '); % remove a bunch of labels from the polar plot; radial and tangential

title('Plot of PD changes (intrinsic)')