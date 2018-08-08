%% Setup
    % Load data
    datadir = '/home/raeed/Projects/limblab/data-td/MultiWorkspace';
    load(sprintf('%s/Han_20171101_TD.mat',datadir))

    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));

    % get only rewards
    [~,td] = getTDidx(td,'result','R');

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

%% Set up helper variables
    % to predict elbow from both neurons and hand
    elbow_idx = 28:30;
    hand_idx = 1:3;
    other_hand_idx = 4:6;

%% Find best lag for elbow to hand by xcorr of the speed profiles
    hand_vel = get_vars(td,{'marker_vel',hand_idx});
    hand_speed = sqrt(sum(hand_vel.^2,2));
    other_hand_vel = get_vars(td,{'marker_vel',other_hand_idx});
    other_hand_speed = sqrt(sum(other_hand_vel.^2,2));
    elbow_vel = get_vars(td,{'marker_vel',elbow_idx});
    elbow_speed = sqrt(sum(elbow_vel.^2,2));

    maxlag = 10;
    timevec = linspace(-td(1).bin_size*maxlag,td(1).bin_size*maxlag,2*maxlag+1);
    hand_x_elbow = xcorr(hand_speed,elbow_speed,maxlag);
    hand_x_hand = xcorr(hand_speed,other_hand_speed,maxlag);

    [~,maxcorr_idx] = max(hand_x_elbow);
    [~,maxcorr_hand_idx] = max(hand_x_hand);

    figure
    subplot(2,1,1)
    plot(timevec,hand_x_elbow,'o-k','linewidth',2)
    hold on
    ylim = get(gca,'ylim');
    plot(repmat(timevec(maxcorr_idx),1,2),ylim,'r','linewidth',2)
    set(gca,'box','off','tickdir','out')
    title 'Cross-correlation between hand and elbow speed'

    subplot(2,1,2)
    plot(timevec,hand_x_hand,'o-k','linewidth',2)
    hold on
    ylim = get(gca,'ylim');
    plot(repmat(timevec(maxcorr_hand_idx),1,2),ylim,'r','linewidth',2)
    set(gca,'box','off','tickdir','out')
    title 'Cross-correlation between hand and hand'
    xlabel 'Lag time (s)'


