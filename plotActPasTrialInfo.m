%% Set up meta info
    if ispc
        homefolder = 'C:\Users\rhc307';
    else
        homefolder = '/home/raeed';
    end
    
    % for data loading
    datadir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','td-library');
    file_info = dir(fullfile(datadir,'*COactpas*'));
    filenames = horzcat({file_info.name})';
    arrayname = 'S1';
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Separability';
    run_date = char(datetime('today','format','yyyyMMdd'));

    % misc
    monkey_names = {'Chips';'Han'};
    session_dates = {{'2017/9/12','2017/9/13'};{'2017/1/5','2017/2/3'}};
    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Loop through files
    pdTable = cell(4,1);
    for filenum = 1:4
        %% load and preprocess data
            td = load(fullfile(datadir,[filenames{filenum}]));
    
            % rename trial_data for ease
            td = td.trial_data;
    
            % first process marker data
            % find times when markers are NaN and replace with zeros temporarily
            for trialnum = 1:length(td)
                markernans = isnan(td(trialnum).markers);
                td(trialnum).markers(markernans) = 0;
                td(trialnum) = smoothSignals(td(trialnum),struct('signals','markers'));
                td(trialnum).markers(markernans) = NaN;
                clear markernans
            end
    
            % get marker velocity
            td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
            
            % get speed and ds
            td = getNorm(td,struct('signals','vel','norm_name','speed'));
            td = getDifferential(td,struct('signals','speed','alias','dspeed'));
            
            % remove unsorted neurons
            unit_ids = td(1).S1_unit_guide;
            unsorted_units = (unit_ids(:,2)==0);
            new_unit_guide = unit_ids(~unsorted_units,:);
            for trialnum = 1:length(td)
                td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
                
                spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
                spikes(:,unsorted_units) = [];
                td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
            end
    
            % add firing rates in addition to spike counts
            td = addFiringRates(td,struct('array',arrayname));
    
            % prep trial data by getting only rewards and trimming to only movements
            % split into trials
            td = splitTD(...
                td,...
                struct(...
                    'split_idx_name','idx_startTime',...
                    'linked_fields',{{...
                        'trialID',...
                        'result',...
                        'bumpDir',...
                        'tgtDir',...
                        'ctrHoldBump',...
                        'ctrHold',...
                        }},...
                    'start_name','idx_startTime',...
                    'end_name','idx_endTime'));
            [~,td] = getTDidx(td,'result','R');
            td = reorderTDfields(td);
            
            % clean nans out...?
            nanners = isnan(cat(1,td.tgtDir));
            td = td(~nanners);
            fprintf('Removed %d trials because of missing target direction\n',sum(nanners))
            biggers = cat(1,td.ctrHoldBump) & abs(cat(1,td.bumpDir))>360;
            td = td(~biggers);
            fprintf('Removed %d trials because bump direction makes no sense\n',sum(biggers))
    
            % remove trials where markers aren't present
            bad_trial = false(length(td),1);
            for trialnum = 1:length(td)
                if any(any(isnan(td(trialnum).markers)))
                    bad_trial(trialnum) = true;
                end
            end
            td(bad_trial) = [];
            fprintf('Removed %d trials because of missing markers\n',sum(bad_trial))
            
            % remove trials where muscles aren't present
            bad_trial = false(length(td),1);
            for trialnum = 1:length(td)
                if any(any(isnan(td(trialnum).muscle_len) | isnan(td(trialnum).muscle_vel)))
                    bad_trial(trialnum) = true;
                end
            end
            td(bad_trial) = [];
            fprintf('Removed %d trials because of missing muscles\n',sum(bad_trial))
            
            % for Chips_20170912, file has bumps on 100% of trials
            if strcmpi(td(1).monkey,'Chips') && contains(td(1).date_time,'2017/9/12')
                td_copy = td;
                [td_copy.ctrHoldBump] = deal(false);
                td = cat(2,td,td_copy);
                clear td_copy
            end
            
            % split into active and passive
            [~,td_act] = getTDidx(td,'ctrHoldBump',false);
            [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    
            % find the relevant movmement onsets
            td_act = getMoveOnsetAndPeak(td_act,struct(...
                'start_idx','idx_goCueTime',...
                'start_idx_offset',20,...
                'peak_idx_offset',20,...
                'end_idx','idx_endTime',...
                'method','peak',...
                'peak_divisor',10,...
                'min_ds',1));
            td_pas = getMoveOnsetAndPeak(td_pas,struct(...
                'start_idx','idx_bumpTime',...
                'start_idx_offset',-5,... % give it some wiggle room
                'peak_idx_offset',-5,... % give it some wiggle room
                'end_idx','idx_goCueTime',...
                'method','peak',...
                'peak_divisor',10,...
                'min_ds',1));
            % throw out all trials where bumpTime and movement_on are more than 3 bins apart
            bad_trials = isnan(cat(1,td_pas.idx_movement_on)) | abs(cat(1,td_pas.idx_movement_on)-cat(1,td_pas.idx_bumpTime))>3;
            td_pas = td_pas(~bad_trials);
    
            % even out sizes and put back together
            minsize = min(length(td_act),length(td_pas));
            td_act = td_act(1:minsize);
            td_pas = td_pas(1:minsize);
            td = cat(2,td_act,td_pas);
    
            % trim to just movements
            td = trimTD(td,{'idx_movement_on',-50},{'idx_movement_on',60});
    
        %% Plot out hand speed
            figure('defaultaxesfontsize',18)
            for trial = 1:length(td)
                timevec = ((1:length(td(trial).speed))-td(trial).idx_movement_on)*td(trial).bin_size;
                if td(trial).ctrHoldBump
                    plot(timevec,td(trial).speed,'r')
                else
                    plot(timevec,td(trial).speed,'k')
                end
                hold on
            end
            plot(zeros(2,1),ylim,'--k','linewidth',2)
            hold on
            plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
            xlabel('Time from movement onset (s)')
            ylabel('Hand speed (cm/s)')
            set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
            set(gcf,'renderer','Painters')
            suptitle(strrep(strrep(filenames{filenum},'_COactpas_TD.mat',''),'_','-'))
            % saveas(gcf,fullfile(figdir,sprintf('%s_handspeed_run%s.pdf',strrep(filenames{filenum},'_TD.mat',''),run_date)))
    
        %% plot out average hand speed
            % td_avg = trialAverage(td,struct('conditions','ctrHoldBump','add_std',true));
            % figure
            % for trial = 1:length(td_avg)
            %     if td_avg(trial).ctrHoldBump
            %         plot(td_avg(trial).speed,'r')
            %     else
            %         plot(td_avg(trial).speed,'k')
            %     end
            %     hold on
            % end
            % plot(repmat(td(1).idx_movement_on,2,1),ylim,'--k','linewidth',2)
            % hold on
            % plot(repmat(td(1).idx_movement_on+12,2,1),ylim,'--k','linewidth',2)
            % set(gca,'box','off','tickdir','out')
    
        %% Plot out example rasters for each direction
            dirs = unique(cat(1,td.tgtDir));
            figure('defaultaxesfontsize',18)
            for dirnum = 1:length(dirs)
                % pick a random active and random passive trial with this direction
                act_idx = getTDidx(td,'tgtDir',dirs(dirnum),'ctrHoldBump',false,'rand',1);
                pas_idx = getTDidx(td,'bumpDir',dirs(dirnum),'ctrHoldBump',true,'rand',1);
                td_temp = td([act_idx pas_idx]);
    
                for trialnum = 1:length(td_temp)
                    spikes = getSig(td_temp(trialnum),'S1_spikes')';
                    timevec = ((1:size(spikes,2))-td_temp(trialnum).idx_movement_on)*td_temp(trialnum).bin_size;
                    % active on left, passive on right
                    subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+trialnum)
                    % neurons
                    for neuronnum = 1:size(spikes,1)
                        spike_times = timevec(spikes(neuronnum,:)>0);
                        scatter(spike_times,repmat(neuronnum,size(spike_times)),5,'k','filled')
                        hold on
                    end
                    plot(zeros(1,2),[0 size(spikes,1)+1],'--k')
                    plot(ones(1,2)*0.12,[0 size(spikes,1)+1],'--k')
                    xlabel('Time from movement onset (s)')
                    set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5],'ytick',[])
                end
                subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+1)
                ylabel(sprintf('Direction %f',dirs(dirnum)))
            end
            suptitle(strrep(strrep(filenames{filenum},'_COactpas_TD.mat',''),'_','-'))
            % saveas(gcf,fullfile(figdir,sprintf('%s_exampleraster_run%s.pdf',strrep(filenames{filenum},'_TD.mat',''),run_date)))
    
        %% Plot out average rasters
            % dirs = unique(cat(1,td.tgtDir));
            % act_idx = getTDidx(td,'ctrHoldBump',false);
            % pas_idx = getTDidx(td,'ctrHoldBump',true);
            % td_act_avg = trialAverage(td(act_idx),'tgtDir');
            % td_pas_avg = trialAverage(td(pas_idx),'bumpDir');
            % clim = [min(min(getSig(td_avg,'S1_spikes'))) max(max(getSig(td_avg,'S1_spikes')))];
            % figure('defaultaxesfontsize',18)
            % for dirnum = 1:length(dirs)
            %     % pick a random active and random passive trial with this direction
            %     act_idx = getTDidx(td_act_avg,'tgtDir',dirs(dirnum));
            %     pas_idx = getTDidx(td_pas_avg,'bumpDir',dirs(dirnum));
            %     td_temp = td([act_idx pas_idx]);
    
            %     for trialnum = 1:length(td_temp)
            %         spikes = getSig(td_temp(trialnum),'S1_spikes')';
            %         timevec = ((1:size(spikes,2))-td_temp(trialnum).idx_movement_on)*td_temp(trialnum).bin_size;
            %         % active on left, passive on right
            %         subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+trialnum)
            %         % neurons
            %         imagesc(spikes,clim)
            %         hold on
            %         plot(repmat(td_temp(trialnum).idx_movement_on,1,2),[0 size(spikes,1)+1],'--w')
            %         plot(repmat(td_temp(trialnum).idx_movement_on+12,1,2),[0 size(spikes,1)+1],'--w')
            %         xlabel('Time from movement onset (s)')
            %         set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
            %     end
            %     subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+1)
            %     ylabel(sprintf('Direction %f',dirs(dirnum)))
            % end
            % suptitle(strrep(strrep(filenames{filenum},'_COactpas_TD.mat',''),'_','-'))
    
        %% Get PDs for a scatter plot
            td_tuning = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',11});
            td_tuning = binTD(td_tuning,'average');
            % add a direction vector term for PD calculations
            [td_tuning.movevec] = deal([0 0]);
            for trialnum = 1:length(td_tuning)
                if td_tuning(trialnum).ctrHoldBump
                    td_tuning(trialnum).movevec = [cosd(td_tuning(trialnum).bumpDir) sind(td_tuning(trialnum).bumpDir)];
                else
                    td_tuning(trialnum).movevec = [cosd(td_tuning(trialnum).tgtDir) sind(td_tuning(trialnum).tgtDir)];
                end
            end
            act_idx = getTDidx(td_tuning,'ctrHoldBump',false);
            pas_idx = getTDidx(td_tuning,'ctrHoldBump',true);
            pd_params = struct(...
                'do_plot',false,...
                'num_test_dirs',length(unique(cat(1,td_tuning.tgtDir))),...
                'in_signals','movevec',...
                'out_signals','S1_FR',...
                'out_signal_names',td_tuning(1).S1_unit_guide);
            pd_params.prefix = 'act';
            pdTable_act = getTDClassicalPDs(td_tuning(act_idx),pd_params);
            pd_params.prefix = 'pas';
            pdTable_pas = getTDClassicalPDs(td_tuning(pas_idx),pd_params);
            pdTable{filenum} = join(pdTable_act,pdTable_pas);
            % pdTable = join(pdTable_act,pdTable_pas);
            % isTuned = pdTable.act_movevecTuned & pdTable.pas_movevecTuned;
            % figure
            % scatter(180/pi*pdTable.act_movevecPD(isTuned),180/pi*pdTable.pas_movevecPD(isTuned),[],'k','filled')
            % hold on
            % plot([-180 180],[-180 180],'--k')
            % axis equal
    end
    pdTable = vertcat(pdTable{:});

%% make PD plot
    % set xvals for histogram dot plot
    num_categories = 4;
    category_x = (2:3:((num_categories-1)*3+2))/10;

    % start figure
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:length(monkey_names)
        % set template values for sessions based on monkey
        template_x = linspace(-0.3,0.3,length(session_dates{monkeynum}))/10;
        hist_xvals = category_x'+template_x;

        subplot(2,length(monkey_names),monkeynum)
        plot([-180 180],[0 0],'-k','linewidth',2)
        hold on
        plot([0 0],[-180 180],'-k','linewidth',2)
        plot([-180 180],[-180 180],'--k','linewidth',2)
        for sessionnum = 1:length(session_dates{monkeynum})
            [~,session_table] = getNTidx(pdTable,'monkey',monkey_names{monkeynum},'date',session_dates{monkeynum}{sessionnum});
            isTuned = session_table.act_movevecTuned & session_table.pas_movevecTuned;
            
            actPD = 180/pi*session_table.act_movevecPD;
            pasPD = 180/pi*session_table.pas_movevecPD;

            actCI = 180/pi*session_table.act_movevecPDCI;
            pasCI = 180/pi*session_table.pas_movevecPDCI;
            actCI(actCI(:,2)<actCI(:,1),2) = actCI(actCI(:,2)<actCI(:,1),2)+360;
            pasCI(pasCI(:,2)<pasCI(:,1),2) = pasCI(pasCI(:,2)<pasCI(:,1),2)+360;

            % isTuned = diff(actCI,1,2)<=90 & diff(pasCI,1,2)<=90;
            % make scatter plot
            subplot(2,length(monkey_names),monkeynum)
            scatter(actPD(isTuned),pasPD(isTuned),[],session_colors(sessionnum,:),'filled')
            plot(repmat(actPD(isTuned),1,2)',pasCI(isTuned,:)','-','linewidth',1,'color',session_colors(sessionnum,:))
            plot(actCI(isTuned,:)',repmat(pasPD(isTuned),1,2)','-','linewidth',1,'color',session_colors(sessionnum,:))
            % plot CIs again for wrap around
            plot(repmat(actPD(isTuned),1,2)',pasCI(isTuned,:)'-360,'-','linewidth',1,'color',session_colors(sessionnum,:))
            plot(actCI(isTuned,:)'-360,repmat(pasPD(isTuned),1,2)','-','linewidth',1,'color',session_colors(sessionnum,:))

            % make dot histogram of neural tuning properties
            subplot(2,length(monkey_names),length(monkey_names)+monkeynum)
            num_act = sum(session_table.act_movevecTuned & ~session_table.pas_movevecTuned);
            num_pas = sum(~session_table.act_movevecTuned & session_table.pas_movevecTuned);
            num_non = sum(~session_table.act_movevecTuned & ~session_table.pas_movevecTuned);
            num_actpas = sum(isTuned);
            session_x = hist_xvals(:,sessionnum);
            scatter(repmat(session_x(1),1,num_act),1:num_act,[],session_colors(sessionnum,:),'filled')
            hold on
            scatter(repmat(session_x(2),1,num_pas),1:num_pas,[],session_colors(sessionnum,:),'filled')
            scatter(repmat(session_x(3),1,num_actpas),1:num_actpas,[],session_colors(sessionnum,:),'filled')
            scatter(repmat(session_x(4),1,num_non),1:num_non,[],session_colors(sessionnum,:),'filled')
        end
        % set axes
        subplot(2,length(monkey_names),monkeynum)
        axis equal
        set(gca,'box','off','tickdir','out','xlim',[-180 180],'ylim',[-180 180],'xtick',[-180 180],'ytick',[-180 180])
        title(monkey_names{monkeynum})
        xlabel('Active PD')
        ylabel('Passive PD')

        % plot out the number of tuned neurons
        subplot(2,length(monkey_names),length(monkey_names)+monkeynum)
        set(gca,'box','off','tickdir','out',...
            'xlim',[category_x(1)+template_x(1)-0.1 category_x(end)+template_x(end)+0.1],...
            'xtick',category_x,...
            'ylim',[0 25],...
            'ytick',0:5:25,...
            'xticklabel',{'Active-only','Passive-only','Active-Passive','None'})
        ylabel('Number of neurons')
    end

