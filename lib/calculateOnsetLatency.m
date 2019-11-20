function [latency_table] = calculateOnsetLatency(trial_data,params)
%CALCULATEONSETLATENCY finds the changepoint of signals, averaged over trials, related to some trial event
%   Based on an algorithm in West RW, Ogden RT. Continuous-time estimation of a change-point in a Poisson process. J Stat Comp 56: 293–302, 1997.
%   This function will average over the given conditions to find the earliest
%   change point in the signals across the conditions. It returns a NeuronTable
%   Inputs:
%       trial_data - the data to check
%       params - parameters structure with fields:
%           signals - (char array or cell array) Signals to check (ex. 'S1_spikes' or 'S1_FR')
%           out_signal_names - names of signals to be used in the output table
%           conditions - (char array or cell array) field names in trial_data to average over
%           ref_idx - (char array) name of idx for event to reference latency to
%           meta - (struct) any conditions that should be tagged in the output table
%               (ex. struct('isPassive',1) for passive trials)
%   Outputs:
%       latency_table - (table) neuron-table containing the earliest change point for neurons
%           across the conditions averaged over, along with the 95% confidence intervals,
%           calculated by bootstrapping. Table will have columns designating the trial_data
%           meta information.

% parameter setup
signals = '';
out_signal_names = '';
conditions = 'all';
ref_idx = 'idx_startTime';
meta = [];
assignParams(who,params)

% trial average over conditions
td_avg = trialAverage(trial_data,struct('conditions',conditions));

latency_table = makeNeuronTableStarter(trial_data,struct('out_signal_names',out_signal_names,'meta',meta));
end

