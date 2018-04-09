%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function curveTable = getCurvesFromCrossval(crossvalTable)
%
%   Gets curve table from crossvalidation table in analyzeTRT
%
% INPUTS:
%   crossvalTable : the neuron table that contains weights calculated from model
%                   and curves from model
%
% OUTPUTS:
%   curveTable : calculated curve table
%
% Written by Raeed Chowdhury. Updated Nov 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function curveTable = getCurvesFromCrossval(crossvalTable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get key column indices
key_cols = contains(crossvalTable.Properties.VariableDescriptions,'meta');

% get average table first
avgTable = neuronAverage(crossvalTable,key_cols);

% get curve column indices
curve_cols = endsWith(avgTable.Properties.VariableNames,'Curve');
ciLo_cols = endsWith(avgTable.Properties.VariableNames,'CurveCILo');
ciHi_cols = endsWith(avgTable.Properties.VariableNames,'CurveCIHi');
% get key column indices again
key_cols = contains(avgTable.Properties.VariableDescriptions,'meta');
% get bins column
bin_cols = endsWith(avgTable.Properties.VariableNames,'bins');


% Extract only curve and CI columns
curveTable = avgTable(:,key_cols | bin_cols | curve_cols | ciLo_cols | ciHi_cols);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
