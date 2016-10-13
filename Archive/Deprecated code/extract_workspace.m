function [times] = extract_workspace(bdf,bottom_left,top_right)
% Extract workspace from BDF with corners given by bottom_left and
% top_right (coordinates are in [x y] format). Outputs new bdf and set of
% times [start_times stop_times] for reaches in workspace

    if(length(bottom_left(:))~=2 || length(top_right(:))~=2)
        error('Coordinates are not in correct format')
    end
    
    bdf_new = bdf;
    t = bdf.pos(:,1);
    
    ind = bdf.pos(:,2)>bottom_left(1) & bdf.pos(:,2)<top_right(1) & bdf.pos(:,3)>bottom_left(2) & bdf.pos(:,3)<top_right(2);
    
    % find starts and stops of all reaches in each workspace
    iStart = find(diff(ind)>0);
    iStop  = find(diff(ind)<0);
    % A little Kluge to eleminate any partial trajectories at the beginning or
    % end of the file
    if iStart(1) > iStop(1)
        iStop = iStop(2:end);
    end

    if length(iStart) > length(iStop)
        iStart = iStart(1:length(iStop));
    end
    
    % discard some paths less than minimum length
    lMin = 3;
    keepers = true(size(iStart));
    for i = 1:length(keepers)
        snip = bdf_new.pos(iStart(i):iStop(i), 2:3);

        % Reject paths that are too short
        steps = diff(snip);
        len = sum(sqrt(steps(:,1).^2+steps(:,2).^2));
        keepers(i) = keepers(i) & len > lMin;
    end
    
    % Dump all the rejected trajectories
    iStart = iStart(keepers);
    iStop = iStop(keepers);
    
%     % compile new fields for bdf_new
%     fr = zeros(length(iStart),length(unit_list(bdf)));
%     new_pos = [];
%     new_vel = [];
%     new_acc = [];
%     new_force = [];
%     new_good_flag = [];
%     for i=1:length(iStart)
%         new_pos = [new_pos;bdf.pos(iStart(i):iStop(i),:)];
%         new_vel = [new_vel;bdf.vel(iStart(i):iStop(i),:)];
%         new_acc = [new_acc;bdf.acc(iStart(i):iStop(i),:)];
%         if(isfield(bdf,'force'))
%             new_force = [new_force;bdf.force(iStart(i):iStop(i),:)];
%         end
%         if(isfield(bdf,'good_kin_data'))
%             new_good_flag = [new_good_flag;bdf.good_kin_data(iStart(i):iStop(i),:)];
%         end
%     end
%     bdf_new.pos = new_pos;
%     bdf_new.vel = new_vel;
%     bdf_new.acc = new_acc;
%     if(isfield(bdf,'force'))
%         bdf_new.force = new_force;
%     end
%     if(isfield(bdf,'good_kin_data'))
%         bdf_new.good_kin_data = new_good_flag;
%     end
%     
%     % get new spikes
%     for uid = 1:length(bdf.units)
%         spikes = [];
%         s = bdf.units(uid).ts;
%         for i=1:length(iStart)
%             spikes = [spikes;s(s < t(iStop(i)) & s > t(iStart(i)))];
%         end
%         bdf_new.units(uid).ts = spikes;
%     end % foreach unit
    
    % extract actual times of workspaces [start_times stop_times]
    times = [t(iStart,1) t(iStop,1)];
end