%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pdTable = getTDPDs(trial_data,params)
%
%   Gets PD table for given out_signal. You need to define the out_signal
% and move_corr parameters at input.
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .out_signals  : which signals to calculate PDs for
%       .trial_idx    : (NOT IMPLEMENTED) trials to evaluate. Ways to use:
%                     1) 1:end treats each trial separately
%                     2) 1:N:end predicts in bins of size N trials
%                     3) [1,end] returns a single value for predicting all trials
%                         DEFAULT: [1,length(trial_data]
%       .move_corr    : (string) name of behavior correlate for PD
%                           'vel' : velocity of handle
%                           'acc' : acceleration of handle
%                           'force'  : force on handle
%       .block_trials : (NOT IMPLEMENTED) if true, takes input of trial indices and pools
%                       them together for a single eval. If false, treats the trial indices
%                       like a list of blocked testing segments
%       .num_boots    : # bootstrap iterations to use (if <2, doesn't bootstrap)
%
% OUTPUTS:
%   pdTable : calculated velocity PD table with CIs
%                Note: will return relative metric if model_name is 1x2 cell of names
%
% Written by Raeed Chowdhury. Updated Jul 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pdTable = getTDPDs(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT PARAMETERS
out_signals      =  [];
trial_idx        =  [1,length(trial_data)];
move_corr      =  '';
block_trials     =  false;
num_boots        =  1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some undocumented parameters
td_fn_prefix     =  '';    % prefix for fieldname
if nargin > 1, assignParams(who,params); end % overwrite parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
possible_corrs = {'vel','acc','force'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process inputs
if isempty(out_signals), error('Need to provide output signal'); end
if isempty(move_corr), error('Must provide movement correlate.'); end
if ~any(ismember(move_corr,possible_corrs)), error('Correlate not recognized.'); end
out_signals = check_signals(trial_data(1),out_signals);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PD
response_var = get_vars(trial_data,out_signals);
bootfunc = @(data) fitglm(data(:,2:end),data(:,1),'Distribution','Poisson');
tic;
for uid = 1:size(response_var,2)
    disp(['  Bootstrapping GLM PD computation(ET=',num2str(toc),'s).'])
    %bootstrap for firing rates to get output parameters
    if block_trials
        % not implemented currently, look at evalModel for how block trials should be implemented
        error('getTDPDs:noBlockTrials','Block trials option is not implemented yet')
    else
        data_arr = [response_var(:,uid) cat(1,trial_data.(move_corr))];
        boot_tuning = bootstrp(num_boots,@(data) {bootfunc(data)}, data_arr);
        boot_coef = cell2mat(cellfun(@(x) x.Coefficients.Estimate',boot_tuning,'uniformoutput',false));

        if size(boot_coef,2) ~= 3
            error('getTDPDs:moveCorrProblem','GLM doesn''t have correct number of inputs')
        end

        dirs = atan2(boot_coef(:,3),boot_coef(:,2));
        %handle wrap around problems:
        centeredDirs=dirs-mean(dirs);
        while(sum(centeredDirs<-pi))
            centeredDirs(centeredDirs<-pi) = centeredDirs(centeredDirs<-pi)+2*pi;
        end
        while(sum(centeredDirs>pi))
            centeredDirs(centeredDirs>pi) = centeredDirs(centeredDirs>pi)-2*pi;
        end
        dirArr(uid,:)=mean(dirs);
        dirCIArr(uid,:)=prctile(centeredDirs,[2.5 97.5])+mean(dirs);
    end
end

% package output
pdTable = table(dirArr,dirCIArr,'VariableNames',{[move_corr 'Dir'],[move_corr 'DirCI']});



end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
