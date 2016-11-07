function [predOutputTable ,pseudoR2] = neuronPredict(neuronModel,testDataTable)
% NEURONPREDICT Function to predict neural output given neural models
% (usually GLM). Predict using the test data table in the same format as
% the table used to train it.

neurIDX = contains(testDataTable.Properties.VariableNames,'LeftS1Area2CH');
predOutputTable = testDataTable;
predNeur = zeros(height(testDataTable),sum(neurIDX));
pseudoR2 = array2table(zeros(1,sum(neurIDX)),'VariableNames',testDataTable.Properties.VariableNames(neurIDX));
trueNeur = testDataTable{:,neurIDX};

for neuronctr = 1:length(neuronModel)
    model = neuronModel{neuronctr};
    predNeur(:,neuronctr) = model.predict(testDataTable);
    
    pseudoR2{1,neuronctr} = compute_pseudo_R2(trueNeur(:,neuronctr),predNeur(:,neuronctr),mean(trueNeur(:,neuronctr)));
end
    
predOutputTable{:,neurIDX} = predNeur;

end