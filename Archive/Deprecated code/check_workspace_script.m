%% test split workspaces
figure
rose(atan2(vel(:,3),vel(:,2)))
title 'Movement direction distribution of whole workspace'

DL_ind = pos(:,2)<mean(pos(:,2)) & pos(:,3)>mean(pos(:,3));
PM_ind = pos(:,2)>mean(pos(:,2)) & pos(:,3)<mean(pos(:,3));
figure
subplot 221
rose(atan2(vel(DL_ind,3),vel(DL_ind,2)))
subplot 224
rose(atan2(vel(PM_ind,3),vel(PM_ind,2)))

%%
figure
plot(pos(:,2),pos(:,3))