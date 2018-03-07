%% Estimate Jacobian for each epoch
    % prep trial data by getting only rewards and trimming to only movements
    [~,td] = getTDidx(trial_data,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % bin data at 50ms
    td = binTD(td,5);
    
    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);
    
    % get model of muscles to hand (fit Jacobian)
    muscVel_pm = get_vars(td_pm,{'opensim',find(contains(td(1).opensim_names,'_muscVel'))});
    muscVel_dl = get_vars(td_dl,{'opensim',find(contains(td(1).opensim_names,'_muscVel'))});
    handVel_pm = get_vars(td_pm,{'opensim',find(contains(td(1).opensim_names,'_handVel'))});
    handVel_dl = get_vars(td_dl,{'opensim',find(contains(td(1).opensim_names,'_handVel'))});
    handleVel_pm = get_vars(td_pm,{'vel',1:2})/100; % handleVel is in cm/s in TD
    handleVel_dl = get_vars(td_dl,{'vel',1:2})/100;
    
    % muscVel = handVel*J
    J_pm = handVel_pm\muscVel_pm;
    J_dl = handVel_dl\muscVel_dl;
    
    % rearrange to only get x and y components
    J_pm = J_pm([3;1],:);
    J_dl = J_dl([3;1],:);
    
    % get change in Jacobian
    delJ = sum(sqrt(sum((J_dl-J_pm).^2,1)));
    