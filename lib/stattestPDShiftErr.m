function [h,p] = stattestPDShiftErr(err,modelcompare,tails,num_repeats,num_folds)
% run a bonferroni corrected statistical test on models assuming repeaded k-fold
% crossvalidation

assert(size(modelcompare,2) == 2,'Cannot compare more than two models at a time');
assert(size(modelcompare,1) == size(tails,1),'Tails needs to be same length as modelcompare')

% bonferroni correction
bonferroni_factor = size(modelcompare,1);
alpha = 0.05/bonferroni_factor;

% correction for repeated k-fold cross-val being correlated samples
correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);

h = false(size(modelcompare,1),1);
p = zeros(size(modelcompare,1),1);
for i = 1:size(modelcompare,1)
    diffstat = err.(modelcompare{i,1})-err.(modelcompare{i,2});
    std_err = sqrt(var(diffstat)*correction);
    t_score = mean(diffstat)/std_err;

    switch(tails{i})
    case 'left'
        p(i) = tcdf(t_score,num_folds*num_repeats-1);
    case 'right'
        p(i) = tcdf(-t_score,num_folds*num_repeats-1);
    case 'both'
        p(i) = 2*tcdf(-abs(t_score),num_folds*num_repeats-1);
    otherwise
        error('Unexpected value for tail: %s',tails{i})
    end

    h(i) = p(i)<alpha;
end
