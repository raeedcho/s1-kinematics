%%
figure
for i = 1:height(output_data_pred.muscle_tuning_DL)
    clf
    maxFR = max([output_data_pred.muscle_curves_PM(i,:).CI_high{1};output_data_pred.muscle_curves_DL(i,:).CI_high{1}]);
    plot_tuning(output_data_pred.muscle_tuning_PM(i,:),output_data_pred.muscle_curves_PM(i,:),maxFR,[0.6 0.5 0.7])
    plot_tuning(output_data_pred.muscle_tuning_DL(i,:),output_data_pred.muscle_curves_DL(i,:),maxFR,[1 0 0])
    waitforbuttonpress
end