function plotDecodeAccuracy(decode_accuracy_table)
    % This function plots target classification accuracy from kinematics and neurons

    % set up with parameters that we know (hard-coded in the decoder calculation)
    timevec = -150:50:150;
    timevec = timevec+25;

    % figure out monkeys
    monkey_names = unique(decode_accuracy_table.monkey);

    % plot out active decodes
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:length(monkey_names)
        % figure out what sessions we have for this monkey
        [~,monkey_decodes] = getNTidx(decode_accuracy_table,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_decodes.date);

        for sessionnum = 1:length(session_dates)
            % if there aren't the same number of sessions for each monkey, this might look weird
            subplot(length(monkey_names),length(session_dates),length(monkey_names)*(monkeynum-1)+sessionnum)
            [~,session_decodes] = getNTidx(monkey_decodes,'date',session_dates{sessionnum});

            % plot kinematic decodes
            plot([0 0],[0 1],'-k','linewidth',2)
            hold on
            plot([timevec(1) timevec(end)],[0.25 0.25],'--k','linewidth',2)
            plot(repmat(timevec,height(session_decodes),1)',session_decodes.kin_accuracy_act','-k')
            plot(repmat(timevec,height(session_decodes),1)',session_decodes.neur_accuracy_act','-b')
            set(gca,'box','off','tickdir','out','ylim',[0 1])
            xlabel('Time from movement onset (ms)')
            ylabel('Classification accuracy')
            title(sprintf('%s %s',monkey_names{monkeynum},session_dates{sessionnum}))
        end
    end
    suptitle('Active trials')
    
    % plot out passive decodes
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:length(monkey_names)
        % figure out what sessions we have for this monkey
        [~,monkey_decodes] = getNTidx(decode_accuracy_table,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_decodes.date);

        for sessionnum = 1:length(session_dates)
            % if there aren't the same number of sessions for each monkey, this might look weird
            subplot(length(monkey_names),length(session_dates),length(monkey_names)*(monkeynum-1)+sessionnum)
            [~,session_decodes] = getNTidx(monkey_decodes,'date',session_dates{sessionnum});

            % plot kinematic decodes
            plot([0 0],[0 1],'-k','linewidth',2)
            hold on
            plot([timevec(1) timevec(end)],[0.25 0.25],'--k','linewidth',2)
            plot(repmat(timevec,height(session_decodes),1)',session_decodes.kin_accuracy_pas','-k')
            plot(repmat(timevec,height(session_decodes),1)',session_decodes.neur_accuracy_pas','-b')
            set(gca,'box','off','tickdir','out','ylim',[0 1])
            xlabel('Time from movement onset (ms)')
            ylabel('Classification accuracy')
            title(sprintf('%s %s',monkey_names{monkeynum},session_dates{sessionnum}))
        end
    end
    suptitle('Passive trials')

    % plot out active v passive neural decode accuracy
    figure('defaultaxesfontsize',18)
    alpha = 0.05;
    for monkeynum = 1:length(monkey_names)
        % figure out what sessions we have for this monkey
        [~,monkey_decodes] = getNTidx(decode_accuracy_table,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_decodes.date);

        for sessionnum = 1:length(session_dates)
            % if there aren't the same number of sessions for each monkey, this might look weird
            subplot(length(monkey_names),length(session_dates),length(monkey_names)*(monkeynum-1)+sessionnum)
            [~,session_decodes] = getNTidx(monkey_decodes,'date',session_dates{sessionnum});

            % figure out confidence intervals
            num_repeats = double(max(session_decodes.crossvalID(:,1)));
            num_folds = double(max(session_decodes.crossvalID(:,2)));
            crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
            mean_acc_act = mean(session_decodes.neur_accuracy_act);
            mean_acc_pas = mean(session_decodes.neur_accuracy_pas);
            var_acc_act = var(session_decodes.neur_accuracy_act);
            var_acc_pas = var(session_decodes.neur_accuracy_pas);
            upp = tinv(1-alpha/2,num_folds*num_repeats-1);
            low = tinv(alpha/2,num_folds*num_repeats-1);
            CI_lo_act = mean_acc_act + low * sqrt(crossval_correction*var_acc_act);
            CI_lo_pas = mean_acc_pas + low * sqrt(crossval_correction*var_acc_pas);
            CI_hi_act = mean_acc_act + upp * sqrt(crossval_correction*var_acc_act);
            CI_hi_pas = mean_acc_pas + upp * sqrt(crossval_correction*var_acc_pas);

            % plot kinematic decodes
            plot([0 0],[0 1],'-k','linewidth',2)
            hold on
            plot([timevec(1) timevec(end)],[0.25 0.25],'--k','linewidth',2)
            % plot(repmat(timevec,height(session_decodes),1)',session_decodes.neur_accuracy_act','-k')
            % plot(repmat(timevec,height(session_decodes),1)',session_decodes.neur_accuracy_pas','-r')
            patch([timevec fliplr(timevec)],[CI_lo_act fliplr(CI_hi_act)],'k','facealpha',0.2)
            patch([timevec fliplr(timevec)],[CI_lo_pas fliplr(CI_hi_pas)],'r','facealpha',0.2)
            plot(timevec,mean_acc_act,'-k')
            plot(timevec,mean_acc_pas,'-r')
            set(gca,'box','off','tickdir','out','ylim',[0 1])
            xlabel('Time from movement onset (ms)')
            ylabel('Classification accuracy')
            title(sprintf('%s %s',monkey_names{monkeynum},session_dates{sessionnum}))
        end
    end
    suptitle('Neural direction classification accuracy')

