function  [outliers] = scatterlims(xlims,ylims,xdata,ydata,size,color,varargin)
%SCATTERLIMS plots a scatterplot where any data outside xlims and ylims is indicated by an arrow inside the limits
%   Detailed explanation goes here

    idx_outliers_x = xdata > xlims(2) | xdata < xlims(1);
    idx_outliers_y = ydata > ylims(2) | ydata < ylims(1);
    
    scatter(...
        xdata(~idx_outliers_x & ~idx_outliers_y),...
        ydata(~idx_outliers_x & ~idx_outliers_y),...
        size,color,varargin{:})
    hold on
    
    xrange = diff(xlims);
    yrange = diff(ylims);
    
    % plot x outliers
    % for all xdata < xlims(1) & ~idx_outliers_y
    idx_to_plot = xdata < xlims(1) & ~idx_outliers_y;
    plot(...
        repmat([xlims(1)+0.05*xrange xlims(1)],sum(idx_to_plot),1)',...
        repmat(ydata(idx_to_plot),1,2)',...
        'linewidth',2,'color',color)
    % for all xdata > xlims(2) & ~idx_outliers_y
    idx_to_plot = xdata > xlims(2) & ~idx_outliers_y;
    plot(...
        repmat([xlims(2)-0.05*xrange xlims(2)],sum(idx_to_plot),1)',...
        repmat(ydata(idx_to_plot),1,2)',...
        'linewidth',2,'color',color)
    % for all ydata > ylims(2) & ~idx_outliers_x
    idx_to_plot = ydata > ylims(2) & ~idx_outliers_x;
    plot(...
        repmat(xdata(idx_to_plot),1,2)',...
        repmat([ylims(2)-0.05*yrange ylims(2)],sum(idx_to_plot),1)',...
        'linewidth',2,'color',color)
    % for all ydata < ylims(1) & ~idx_outliers_x
    idx_to_plot = ydata < ylims(1) & ~idx_outliers_x;
    plot(...
        repmat(xdata(idx_to_plot),1,2)',...
        repmat([ylims(1)+0.05*yrange ylims(1)],sum(idx_to_plot),1)',...
        'linewidth',2,'color',color)
    
    if any(idx_outliers_x & idx_outliers_y)
        scatter(xlims(1),ylims(1),size,color,'^',varargin{:})
        
        outliers = struct();
        outliers.x = xdata(idx_outliers_x & idx_outliers_y);
        outliers.y = ydata(idx_outliers_x & idx_outliers_y);
    end
end

