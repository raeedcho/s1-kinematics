function [rots,scales,shifted_curves] = getEmpiricalShift(curves)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   comares tuning between different conditions with empirical
%   tuning curves and tries to rotate and scale first curve to
%   match second.
%   Inputs - 
%       curves - cell array of tuning curve tables, one table
%                   per condition
%   Outputs - 
%       rots - column vector of rotations from condition 1 to
%                condition 2
%       scales - column vector of scale factors from condition 1 to
%                condition 2
%       shifted_curves - new tuning curve table containing shifted
%                           and scaled condition1 tuning curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check to make sure inputs are valid
if size(curves) ~= 2
    error('curves must be a 2 element cell array of tuning curve tables')
end

% setup
num_neurons = height(curves{1});
shifted_curves = curves{1};
rots = zeros(num_neurons,1);
scales = ones(num_neurons,1);

% for each neural tuning curve
for neuronidx = 1:num_neurons
    % compute complex curve representation as complex column vectors
    cond1 = curves{1}.binnedResponse(neuronidx,:)'.*exp(1i*curves{1}.bins(neuronidx,:)');
    cond2 = curves{2}.binnedResponse(neuronidx,:)'.*exp(1i*curves{2}.bins(neuronidx,:)');

    % assuming cond2 = cond1*complex_scale + epsilon:
    % -> complex_scale = (cond1'*cond1)\cond1'*cond2
    % by OLS solution with complex variables
    complex_scale = (cond1'*cond1)\cond1'*cond2;
    rots(neuronidx) = angle(complex_scale);
    scales(neuronidx) = abs(complex_scale);

    fit_cond = cond1*complex_scale;
    shifted_curves(neuronidx).bins = angle(fit_cond);
    shifted_curves(neuronidx).binnedResponse = abs(fit_cond);
end
