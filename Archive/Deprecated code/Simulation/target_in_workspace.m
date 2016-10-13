function is_inside = target_in_workspace(target,boundaries)
% TARGET_IN_WORKSPACE checks if target is inside workspace boundaries
%   Outputs-
%       is_inside - logical variable or vector (true if inside, false if not)
%   Inputs-
%       target - [x y] position of target (each row is a different target)
%       boundaries - [xmin xmax ymin ymax] boundaries of workspace

xmin = boundaries(1);
xmax = boundaries(2);
ymin = boundaries(3);
ymax = boundaries(4);

targetx = target(:,1);
targety = target(:,2);

is_inside = targetx<=xmax & targety<=ymax & targetx>=xmin & targety>=ymin;