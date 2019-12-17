function [latency_table] = calculateOnsetLatency(trial_data,params)
    %CALCULATEONSETLATENCY finds the changepoint of signals, averaged over trials, related to some trial event
    %   Based on an algorithm in West RW, Ogden RT. Continuous-time estimation of a change-point in a Poisson process. J Stat Comp 56: 293–302, 1997.
    %   This function will average over the given conditions to find the earliest
    %   change point in the signals across the conditions. It returns a NeuronTable
    %   Inputs:
    %       trial_data - the data to check
    %       params - parameters structure with fields:
    %           signals - (char array or cell array) Signals to check (default: 'S1_spikes')
    %           signal_names - names of signals to be used in the output table
    %           conditions - (char array or cell array) field names in trial_data to average over
    %           ref_idx - (char array) name of idx for event to reference latency to
    %           meta - (struct) any conditions that should be tagged in the output table
    %               (ex. struct('isPassive',1) for passive trials)
    %           num_boots - number of bootstrap samples to take for CI (if <=1, doesn't bootstrap)
    %               (default: 1)
    %           verbose - whether to print out messages while bootstrapping (default: false)
    %   Outputs:
    %       latency_table - (table) neuron-table containing the earliest change point for neurons
    %           across the conditions averaged over, along with the 95% confidence intervals,
    %           calculated by bootstrapping. Table will have columns designating the trial_data
    %           meta information.
    
    % parameter setup
    meta = [];
    num_boots = 1;
    verbose = false;
    assignParams(who,params)
    
    % bootstrap?
    if num_boots > 1
        if verbose
            fprintf('Starting bootstrap...\n')
        end
        boot_latency_cell = cell(num_boots,1);
        boot_tic = tic;
        for bootnum = 1:num_boots
            % choose a random set of trials to work through (with replacement)
            boot_idx = randsample(length(trial_data),length(trial_data),true);

            % run through change point analysis
            if isempty(meta)
                meta = struct();
            end
            meta.bootnum = bootnum;

            boot_latency_cell{bootnum} = get_latency_table(trial_data(boot_idx),params);
            if verbose && mod(bootnum,floor(num_boots/100)) == 0
                fprintf('Evaluated bootstrap %d of %d at time %f\n',bootnum,num_boots,toc(boot_tic))
            end
        end
        boot_latency_table = vertcat(boot_latency_cell{:});

        % calculate mean over all bootnum with CI
        if verbose
            fprintf('Compiling bootstrap table...\n')
        end
        keycols = strcmpi(boot_latency_table.Properties.VariableDescriptions,'meta') & ...
            ~strcmpi(boot_latency_table.Properties.VariableNames,'bootnum');
        latency_table = neuronAverage(boot_latency_table,struct(...
            'keycols',keycols,...
            'do_ci',true));
        if verbose
            fprintf('Done\n')
        end
    else
        % get latency for full data set, without CIs
        latency_table = get_latency_table(trial_data,params);
    end
end

function latency_table = get_latency_table(trial_data,params)
    % mini function to get the latency table on a given set of trials
    signals = 'S1_spikes';
    signal_names = trial_data(1).S1_unit_guide;
    conditions = 'all';
    ref_idx = 'idx_movement_on';
    meta = [];
    assignParams(who,params)

    % parameter check
    assert(~isempty(signal_names),'Must provide signal_names! Usually this is just the unit_guide from trial_data')
    assert(size(signals,1)==1,'Onset latency code currently does not work for more than one signal field')
    if ~iscell(conditions), conditions = {conditions}; end
    signals = check_signals(trial_data,signals);

    % trial average over conditions
    [td_avg,cond_idx] = trialAverage(trial_data,struct('conditions',conditions));

    latency_cell = cell(length(td_avg),1);
    for trialnum = 1:length(td_avg)
        % set up neuron-table here...
        if isempty(meta)
            meta = struct();
        end
        for condnum = 1:length(conditions)
            meta.(conditions{condnum}) = td_avg(trialnum).(conditions{condnum});
        end
        latency_cell{trialnum} = makeNeuronTableStarter(trial_data,struct(...
            'out_signal_names',signal_names,...
            'meta',meta));

        change_point_rel = zeros(length(signals{1,2}),1);
        for signum = 1:length(signals{1,2})
            % get average data, but make sure we have integers by multiplying by number of trials
            % (round to remove floating point div errors)
            data = td_avg(trialnum).(signals{1,1})(:,signals{1,2}(signum));
            data = round(data*length(cond_idx{trialnum}));

            % find the changepoint
            change_point = find_change_point(data);

            % reference change point to ref_idx
            change_point_rel(signum) = change_point - td_avg(trialnum).(ref_idx);
        end
        % load into table (and change to seconds instead of bins)
        change_point_table = table(change_point_rel*td_avg(trialnum).bin_size,...
            'VariableNames',{'onsetLatency'});
        change_point_table.Properties.VariableDescriptions = {'linear'};
        change_point_table.Properties.VariableUnits = {'s'};
        latency_cell{trialnum} = horzcat(latency_cell{trialnum},change_point_table);
    end

    % put it all together
    latency_table = vertcat(latency_cell{:});
end

function change_point = find_change_point(data)
    % finds the change point for given data, assuming Poisson process

    % check if data is vector
    assert(size(data,2)==1 && size(data,1)>1, 'data must be a column vector')

    % grid search over all time points
    change_point = 0;
    max_logL = -inf;
    for timepoint = 1:length(data)-2
        theta0 = mean(data(1:timepoint));
        theta1 = mean(data(timepoint+2:end));
        tau = timepoint + (data(timepoint+1)-theta1)/(theta0-theta1);
        p_tau = rem(tau,1);

        % calculate log likelihood for these ML params for this interval
        logL = -tau*theta0 + sum(data(1:timepoint))*log(theta0) + ...
            -(length(data)-(tau+1))*theta1 + sum(data(timepoint+2:end))*log(theta1) + ...
            -p_tau*theta0 - (1-p_tau)*theta1 + ...
            data(timepoint+1)*log(p_tau*theta0+(1-p_tau)*theta1) - log(sum(factorial(data)));

        if logL > max_logL
            max_logL = logL;
            change_point = tau;
        end
    end
end

