function [winners,model_pairs] = compareEncoderPR2(encoderResults,models)
    % comareEncoderPR2 - get outcomes of all pairwise
    % model comparisons

    alpha = 0.05;

    model_pairs = nchoosek(models,2);
    bonferroni_alpha = alpha/size(model_pairs,1);

    neuron_ids = unique(encoderResults.crossEval.signalID,'rows');

    winners = cell(size(model_pairs,1),size(neuron_ids,1));
    for neuronnum = 1:size(neuron_ids,1)
        [~,neuron_eval] = getNTidx(encoderResults.crossEval,'signalID',neuron_ids(neuronnum,:));
        for pairnum = 1:size(model_pairs,1)
            winners{pairnum,neuronnum} = test_model_pair(...
                neuron_eval,...
                model_pairs(pairnum,:),...
                struct('alpha',bonferroni_alpha,...
                    'num_folds',encoderResults.params.num_folds,...
                    'num_repeats',encoderResults.params.num_repeats));
        end
    end
end

function winner = test_model_pair(neuron_eval,model_pair,params)
    % test a specific model pair for a neuron with a t-test
    % at alpha-level. winner is name of model that wins or
    % empty for a draw

    alpha = 0.05;
    num_folds = 5;
    num_repeats = 20;
    assignParams(who,params);

    crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
    diffstat = neuron_eval.(sprintf('glm_%s_model_eval',model_pair{2}))-neuron_eval.(sprintf('glm_%s_model_eval',model_pair{1}));

    alphaup = 1-alpha/2;
    alphalow = alpha/2;

    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    upp = tinv(alphaup,num_folds*num_repeats-1);
    low = tinv(alphalow,num_folds*num_repeats-1);

    dpR2CI_lo = mudiff + low * sqrt(crossval_correction*vardiff);
    dpR2CI_hi = mudiff + upp * sqrt(crossval_correction*vardiff);
    
    % classify neuron
    if dpR2CI_lo>0
        winner = model_pair{2};
    elseif dpR2CI_hi<0
        winner = model_pair{1};
    else
        winner = '';
    end
end
