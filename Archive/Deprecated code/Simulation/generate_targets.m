function targets = generate_targets(num_targets,boundaries,varargin)
% GENERATE_TARGETS generates a list of target locations in a workspace.
%   targets = GENERATE_TARGETS(num_targets,boundaries,spacing_limits)
%
%   Outputs-
%       targets - 2D array for position of target (x position in first
%       column, y position in second column).
%   Inputs-
%       num_targets - number of targets to generate
%       boundaries - boundaries of workspace defined as [xmin xmax ymin
%       ymax]
%       spacing_limits - minimum and maximum spacing between targets in
%       format [min_spacing max_spacing] (optional)

xmin = boundaries(1);
xmax = boundaries(2);
ymin = boundaries(3);
ymax = boundaries(4);

if length(varargin)<1 % no spacing requirements specified
    target_x = rand(num_targets,1)*(xmax-xmin)+xmin;
    target_y = rand(num_targets,1)*(ymax-ymin)+ymin;
    targets = [target_x target_y];
else
    min_spacing = varargin{1}(1);
    max_spacing = varargin{1}(2);
    
    if(min_spacing>=max_spacing)
        error('min_spacing must be less than max_spacing');
    end
    
    targets = zeros(num_targets,2);
    
    %first target
    target_x = rand(1,1)*(xmax-xmin)+xmin;
    target_y = rand(1,1)*(ymax-ymin)+ymin;
    targets(1,:) = [target_x target_y];
    
    %rest of the targets
    for i = 2:num_targets
        %try a target
        target_x = rand(1,1)*(xmax-xmin)+xmin;
        target_y = rand(1,1)*(ymax-ymin)+ymin;
        temp_target = [target_x target_y];
        target_dist = sqrt(sum((temp_target-targets(i-1,:)).^2));
        
        %check the target
        while(target_dist<min_spacing || target_dist>max_spacing)
            target_x = rand(1,1)*(xmax-xmin)+xmin;
            target_y = rand(1,1)*(ymax-ymin)+ymin;
            temp_target = [target_x target_y];
            target_dist = sqrt(sum((temp_target-targets(i-1,:)).^2));
        end
        
        targets(i,:) = temp_target;
    end
end