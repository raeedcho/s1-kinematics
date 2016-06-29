%%
for i=1:length(bdf_PM)
    subplot(length(bdf_PM),1,i)
    
    pos_PM = sample_bdf_kin(bdf_PM{i},50);

    plot(pos_PM(:,1),pos_PM(:,2))
    axis equal
    title(num2str(i))
end

%%
figure
for i=1:length(bdf_PM)
    subplot(length(bdf_PM),1,i)
    
    pos_PM = sample_bdf_kin(bdf_PM{i},50);

    plot(pos_PM(:,1)>0)
    title(num2str(i))
end

%%
figure
pos_PM = sample_bdf_kin(bdf_PM{7},50);

plot(pos_PM(:,1),pos_PM(:,2))
axis equal