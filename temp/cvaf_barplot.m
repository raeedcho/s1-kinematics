% Make a cVAF bar plot like old times...
sessionnum = 1;
cvaf = [];
for monkeynum = 1:length(monkey_names)
    avg_shift_vaf = neuronAverage(shift_vaf{monkeynum,sessionnum},...
        struct('keycols',{{'monkey','date','task','crossvalID'}},'do_ci',false));
    [~,cols] = ismember(strcat(models_to_plot,'_vaf'),avg_shift_vaf.Properties.VariableNames);
    cvaf = [cvaf avg_shift_vaf{:,cols}];
end
mean_cvaf = mean(cvaf);
var_cvaf = var(cvaf);
correction = 1/100 + 1/4;
std_err_cvaf = sqrt(correction*var_cvaf);
num_cols = length(monkey_names)*length(models_to_plot); % this should be the number columns in cvaf now
model_colors_rep = repmat(getModelColors(models_to_plot),length(monkey_names),1);
% model_x = [1.5 2 2.5 4.5 5 5.5 7.5 8 8.5]/10;
monk_x = (2:3:((length(monkey_names)-1)*3+2))/10;
template_x = linspace(-0.5,0.5,length(models_to_plot))/10;
model_spacing = mode(diff(template_x));
model_x = [];
for i = 1:length(monk_x)
    model_x = [model_x template_x+monk_x(i)];
end
figure('defaultaxesfontsize',18)
for colnum = 1:num_cols
    % scatter(repmat(colnum/10,size(cvaf,1),1),cvaf(:,colnum),50,model_colors_rep(colnum,:),'filled')
    % hold on
    % plot(colnum/10,mean(cvaf(:,colnum)),'.','color',model_colors_rep(colnum,:),'linewidth',3,'markersize',40)
    bar(model_x(colnum),mean(cvaf(:,colnum)),model_spacing,'facecolor',model_colors_rep(colnum,:),'edgecolor','none')
    hold on
    % errorbar(colnum/10,mean_cvaf(:,colnum),std_err_cvaf(:,colnum),'color',model_colors_rep(colnum,:),'linewidth',3)
end
errorbar(model_x,mean_cvaf,std_err_cvaf,'k.','linewidth',3)
set(gca,'tickdir','out','box','off','xtick',monk_x,'xticklabel',{'Monkey H','Monkey C','Monkey L'},'ytick',[0 1])
% axis equal
ylim([0 1])
xlim([0 1])
ylabel('Circular VAF')


