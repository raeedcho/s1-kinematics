function [tunedNeurons] = irisPlot( PM_pdData, DL_pdData, which_neurons )
%IRISPLOT Create iris plot for PM and DL
%   Inputs - PM and DL pdData tables with PDs calculated (extract from
%   binnedData object), logical array of neurons to plot (if tuned)
%   Outputs - logical array of which neurons were plotted, based on
%   original data

%extract relevant information
angsPM = PM_pdData.velDir;
dirCIPM = PM_pdData.velDirCI;
angsDL = DL_pdData.velDir;
dirCIDL = DL_pdData.velDirCI;

% calculate CI widths
DLCIwidth = diff(dirCIDL,1,2); % get CI widths
PMCIwidth = diff(dirCIPM,1,2);
DLCIwidth(DLCIwidth<0) = DLCIwidth(DLCIwidth<0)+2*pi;
PMCIwidth(PMCIwidth<0) = PMCIwidth(PMCIwidth<0)+2*pi;
tunedNeurons = DLCIwidth<pi/4 & PMCIwidth<pi/4;

if(~isempty(which_neurons))
    tunedNeurons = tunedNeurons & which_neurons;
end

angsPMtuned = angsPM(tunedNeurons);
angsDLtuned = angsDL(tunedNeurons);

%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

% plot changes with alpha dependent on CI width
for unit_ctr = 1:length(angsPMtuned)
    h=polar(linspace(angsPMtuned(unit_ctr),angsDLtuned(unit_ctr),2),linspace(0.5,1,2));
    set(h,'linewidth',2,'color',[0.1 0.6 1])
end

%plot circles
h=polar(linspace(-pi,pi,1000),ones(1,1000));
set(h,'linewidth',2,'color',[1 0 0])
hold all
h=polar(linspace(-pi,pi,1000),0.5*ones(1,1000));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])

set(findall(gcf, 'String','  0.2','-or','String','  0.4','-or','String','  0.6','-or','String','  0.8',...
        '-or','String','  1') ,'String', ' '); % remove a bunch of labels from the polar plot; radial and tangential'

end

