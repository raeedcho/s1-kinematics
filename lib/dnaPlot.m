function dnaPlot( PM_pdData, DL_pdData)
% Flattened iris plot
%   Inputs - PM and DL pdData tables with PDs calculated (extract from
%   binnedData object), logical array of neurons to plot (if tuned)
%   Outputs - logical array of which neurons were plotted, based on
%   original data

%extract relevant information
angsPM = PM_pdData.velPD;
dirCIPM = PM_pdData.velPDCI;
angsDL = DL_pdData.velPD;
dirCIDL = DL_pdData.velPDCI;

% check tuned neurons
% isTuned_params = struct('move_corr','vel','CIthresh',pi/3);
% tunedNeurons = checkIsTuned(PM_pdData,isTuned_params)...
%             & checkIsTuned(DL_pdData,isTuned_params);
% 
% if(~isempty(which_neurons))
%     tunedNeurons = tunedNeurons & which_neurons;
% end
% angsPM = angsPM(which_neurons);
% angsDL = angsDL(which_neurons);

%plot backbones
plot([-pi pi],[1 1],'linewidth',4,'color',[1 0 0])
hold all
plot([-pi pi],[0 0],'linewidth',4,'color',[0.6 0.5 0.7])

% plot changes with alpha dependent on CI width
for unit_ctr = 1:length(angsPM)
    % figure out whether to wrap forwards or backwards
    if angsDL(unit_ctr)-angsPM(unit_ctr)>pi
        angsDL(unit_ctr) = angsDL(unit_ctr)-2*pi;
    elseif angsDL(unit_ctr)-angsPM(unit_ctr)<-pi
        angsDL(unit_ctr) = angsDL(unit_ctr)+2*pi;
    end
    plot([angsPM(unit_ctr),angsDL(unit_ctr)],[0 1],'linewidth',2,'color',[0.1 0.6 1]);
    plot([angsPM(unit_ctr)+2*pi,angsDL(unit_ctr)+2*pi],[0 1],'linewidth',2,'color',[0.1 0.6 1]);
    plot([angsPM(unit_ctr)-2*pi,angsDL(unit_ctr)-2*pi],[0 1],'linewidth',2,'color',[0.1 0.6 1]);
end

%plot backbones
plot([-pi pi],[1 1],'linewidth',4,'color',[1 0 0])
hold all
plot([-pi pi],[0 0],'linewidth',4,'color',[0.6 0.5 0.7])
axis equal

set(gca,'box','off','tickdir','out','ytick',[0 1],...
        'yticklabel',{'PM','DL'},'xtick',-pi:pi/2:pi,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},...
        'xlim',[-pi pi])

% swap x and y
view([90 -90])

end

