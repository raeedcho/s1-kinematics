function [CI_lo,CI_hi] = crossval_errorbars(vals,params)
    % this function calculates the 95% confidence bars for a crossvalidation table

    alpha = 0.05;
    num_repeats = 0;
    num_folds = 0;
    assignParams(who,params)

    crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
    meanvals = mean(vals);
    varvals = var(vals);
    upp = tinv(1-alpha/2,num_folds*num_repeats-1);
    low = tinv(alpha/2,num_folds*num_repeats-1);
    CI_lo = meanvals + low * sqrt(crossval_correction*varvals);
    CI_hi = meanvals + upp * sqrt(crossval_correction*varvals);
