%% Load replay data and aggregate them in the Multi_Day structure
for hide = 1
    clr;    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_Replays_*'));   RP = [];
    sessions2include = {...
        };
    for nc = 1:length(FileList)
        if contains(FileList(nc).name,sessions2include)
            cd(FileList(nc).folder);
            load(FileList(nc).name);
            fieldNames = fieldnames(detected_replays);
            if nc == 1
                RP = detected_replays;
            else
                for i = 1:length(fieldNames)
                    RP.(fieldNames{i}) = [RP.(fieldNames{i}) detected_replays.(fieldNames{i})];
                end
            end
            disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
            %         disp([numel(RP.comp) numel(RP.replay_type)]);
            if numel(detected_replays.comp) ~= numel(detected_replays.replay_type)
                disp(['incompatible size: ' FileList(nc).name])
            end
        else
            disp(['omitted: ' FileList(nc).name])
        end
    end
    cd(Folder); clearvars -except RP sessions2include
end
%% load session data 
for hide = 1
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'session_info_*'));   sessions_info={};
    for nc = 1:length(FileList)
        if contains(FileList(nc).name,sessions2include)
            cd(FileList(nc).folder);
            load(FileList(nc).name);
            fieldNames = fieldnames(session_info);
%             if nc == 1
%                 sessions_info = session_info;
%             else
            for i = 1:length(fieldNames)
%                 if strcmp(fieldNames{i},'r')
%                     continue
%                 end
                sessions_info.(fieldNames{i}){nc} = {session_info.(fieldNames{i})};
            end
%             end
        else
            disp(['omitted: ' FileList(nc).name])
        end
    end
    cd(Folder); clearvars -except RP sessions_info
end

%% preparing variables
for hide = 1
%     RP_velocity = RP.velocity;

    RP_flight_duration = cell2mat(RP.flight_duration);
    RP_flight_length = cell2mat(RP.flight_length);
    RP_clus_id = cell2mat(RP.cluster_id);
    RP_corr = RP.corr;
    RP_comp = RP.comp;
    weighted_corr_replay = RP.corr;
    replay_scores_replay = RP.scores;
    slope_replay = RP.slopes;
    avg_jump_distance_replay = RP.avg_j;
    
    
    RP_bat_id = [];
    for i = 1: numel(RP.bat_id)
        RP_bat_id = vertcat (RP_bat_id,  str2num(RP.bat_id{i}));
    end
    

end

%% segmenting replay
for hide = 1
RP_post_trimmed = {};
for i = 1:numel(RP.post)
    current_post = RP.post{i};
    if isempty(current_post)
        continue
    end
    [new_st, new_en] =  horiz_segment_v3(current_post,10); % cut if nan area is 10% of time bins long
    RP_post_trimmed {i} = current_post([new_st:new_en],:);
end
% calculating replay metrics after segmentation
% replay_scores_replay_old = replay_scores_replay;

for hide = 1 % not needed because replay metrics were saved
    weighted_corr_replay_s = [];
    avg_jump_distance_replay_s  = [];
    max_jump_distance_replay_s=[];
    replay_scores_replay_s =  [];
    slope_replay_s = [];
    num_peaks_replay_s = [];
    fitted_y_replay_s = {};
    posterior_spread_replay_s = [];
    com_jump_replay_s = [];
    spatial_coverage_s = [];
    segment_len_frac_s = [];
    horizontal_percent_s = [];
    
    warning('off');
    for i = 1:numel(RP.post)
        current_posteriors = RP_post_trimmed{i};
        [weighted_corr_replay_s(i),max_jump_distance_replay_s(i), avg_jump_distance_replay_s(i), num_peaks_replay_s(i), replay_scores_replay_s(i), slope_replay_s(i),fitted_y_replay_s{i}, posterior_spread_replay_s(i), com_jump_replay_s(i),segment_len_frac_s(i), horizontal_percent_s(i),spatial_coverage_s(i)] = evaluate_candidate_event_v7(current_posteriors);

    end
    warning('on');
    
end
% calculating the duration of events using the number of bins
RP_dur = [];
% for i = 1:numel(RP.post)
for i = 1:numel(RP_post_trimmed)
%     RP_dur = [RP_dur size(RP.post{i},1)*5];
    RP_dur = [RP_dur size(RP_post_trimmed{i},1)*5];
end
figure; hist(RP_dur,100); xlabel ("Replay duration (ms)"); ylabel ("count")
end
    RP_velocity = cell2mat(RP.flight_length) ./(RP_dur/1000);
    RP_velocity(find(RP_velocity == Inf)) = NaN;
%% calculating replay metrics & shuffling [SLOW]
% replay_scores_replay_old = replay_scores_replay;
tic
calculate_shuffle = 0;
replay_padding = 40;
center_only = 1;
% posts = RP_post_trimmed;
posts = RP.post;
if ~calculate_shuffle 
    load("shuffle_sig")
end
for hide = 1 % not needed because replay metrics were saved
    weighted_corr_replay = [];
    avg_jump_distance_replay  = [];
    max_jump_distance_replay=[];
    replay_scores_replay =  [];
    slope_replay = [];
    num_peaks_replay = [];
    fitted_y_replay = {};
    posterior_spread_replay = [];
    com_jump_replay = [];
    spatial_coverage = [];
    segment_len_frac = [];
    horizontal_percent = [];
    % variables for storing p values
    % weighted_corr_p_replay = [];
    % avg_jump_p_replay = [];
    % replay_score_p_replay = [];
    % slope_p_replay = [];
    % weighted_corr_shuffled = [];
    % avg_jump_shuffled =[];
    % slope_shuffled = [];
    % replay_score_shuffled =[];
    
    warning('off');
    for i = 1:numel(posts)
        current_posteriors = posts{i};
%         current_posteriors = RP.post;
        if center_only 
            current_posteriors = current_posteriors(1+replay_padding:size(current_posteriors,1)-replay_padding,:);
        end
        [weighted_corr_replay(i),max_jump_distance_replay(i), avg_jump_distance_replay(i), num_peaks_replay(i), replay_scores_replay(i), slope_replay(i),fitted_y_replay{i}, posterior_spread_replay(i), com_jump_replay(i),segment_len_frac(i), horizontal_percent(i),spatial_coverage(i)] = evaluate_candidate_event_v7(current_posteriors);
        if calculate_shuffle
            [weighted_corr_p_replay(i,:), avg_jump_p_replay(i,:), replay_score_p_replay(i,:), slope_p_replay(i,:), weighted_corr_shuffled(i,:), avg_jump_shuffled(i,:), slope_shuffled(i,:), replay_score_shuffled(i,:)] = shuffle_validation_v2(current_posteriors);
        end
    end
    warning('on');
    if calculate_shuffle
        p_combined = [weighted_corr_p_replay replay_score_p_replay];
        shuffle_sig = max(p_combined')<0.05;
        numel(find(shuffle_sig))/numel(shuffle_sig)
        save('shuffle_sig','shuffle_sig')
    end
    toc
end
%% calculating replay metrics 
for hide = 1
%=== filtering out events that are horizontal lines
replay_horiz = [];
empty_count = 0;
for i = 1:numel(RP.post)
    horiz_perc = 0.1;% if position moved less than this percent in the selected timebins, then erase
    min_dur = 5; % segment only horizontal if it's horizontal by at least this long
    current_post = RP.post{i};
    if isempty(current_post)
        empty_count = empty_count + 1;
        continue
    end
    [~,predicted_position] = max(current_post');
    rescale_pos =  rescale([1:size(current_post,2)]);
    predicted_pos = rescale(predicted_position);
    predicted_pos_diff = diff(predicted_pos);
    predicted_pos_diff = [predicted_pos_diff predicted_pos_diff(end)];
    horiz_segment = zeros (1,size(current_post,1));
    horiz_segment (abs(predicted_pos_diff)<= horiz_perc) = 1;
    
    % rule out spikes that are close to one another briefly
    event_st = strfind(horiz_segment, [0 1]);
    event_en = strfind(horiz_segment, [1 0]);
    if horiz_segment(1)==1
        event_st = [1 event_st];
    end
    if horiz_segment(end)==1
        event_en = [event_en numel(horiz_segment)];
    end
    event_dur = event_en - event_st;
    event_erase_idx = find(event_dur<min_dur);
    for ii = event_erase_idx
        horiz_segment([event_st(ii): event_en(ii)])=0;
    end
    event_erase_idx = find(event_dur>=min_dur);
    for ii = event_erase_idx
        if (abs(mean(predicted_pos_diff([event_st(ii)+1: event_en(ii)])))>horiz_perc*2)
            horiz_segment([event_st(ii): event_en(ii)])=0;
        end
    end
    % calculate percent of time bins that are horizontal
    horiz_sum = sum(horiz_segment);
    perc_horiz = horiz_sum/size(current_post,1);
    if perc_horiz < 0.4
        replay_horiz(i) = 0;
    else
        replay_horiz(i) = 1;
    end
    %    figure; imagesc(current_post'); axis xy; hold on; plot([1:size(current_post,1)], horiz_segment*size(current_post,2));
    
end
%=== not recalculating metrics 
weighted_corr_replay = RP.corr;
avg_jump_distance_replay  = RP.avg_j;
max_jump_distance_replay=RP.max_j;
replay_scores_replay =  RP.scores;
% slope_replay = RP.slopes;
% posterior_spread_replay = [];
% com_jump_replay = [];
end
%% selecting for top replays 
top_replay = RP_clus_id~=1 & abs(weighted_corr_replay) > 0.4 & replay_scores_replay>0.4 & segment_len_frac>0.7 & shuffle_sig;
disp(['Number of top replays:   '  num2str(numel(find(top_replay)))])
disp(['Percent forward replay:   ' num2str(sum(weighted_corr_replay(top_replay)>0)/numel(find(top_replay))) ])
good_replay = RP_clus_id~=1 & abs(weighted_corr_replay) > 0.4 & replay_scores_replay>0.4 & segment_len_frac>0.5 & shuffle_sig;
disp(['Number of good replays:   '  num2str(numel(find(good_replay)))])
disp(['Percent forward replay:   ' num2str(sum(weighted_corr_replay(good_replay)>0)/numel(find(good_replay))) ])
%% general visualization
for hide = 1
    figure('units','normalized','outerposition',[0 .5 0.9 0.5]) ;
    tiledlayout(2,5);
    nexttile; hist(unique(cell2mat(RP.flight_length)),100); xlabel('Average Flight length of all sessions (m)')
    nexttile; hist(unique(cell2mat(RP.flight_duration)),100); xlabel('Average Flight duration of all sessions (s)')
    nexttile; hist(RP.corr,100); xlabel('Weighted Correlation');
    nexttile; hist(replay_scores_replay,100); xlabel("replay scores");
    nexttile; hist(RP.velocity(find(weighted_corr_replay>0.8)),100); xlabel('Replay velocity');
%     nexttile; hist(RP.comp,100); xlabel('Replay compression ');
    
    nexttile; hist(unique(cell2mat(RP.flight_length(top_replay))),100); xlabel('Average Flight length of all sessions (m)')
    nexttile; hist(unique(cell2mat(RP.flight_duration(top_replay))),100); xlabel('Average Flight duration of all sessions (s)')
    nexttile; hist(RP.corr(top_replay),100); xlabel('Weighted Correlation of all significant replays');
    nexttile; hist(replay_scores_replay(top_replay),100); xlabel("replay scores");
    nexttile; hist(RP.velocity(top_replay),100); xlabel('Replay velocity of all significant replays');
%     nexttile; hist(RP.comp(top_replay),100); xlabel('Replay compression of all significant replays');
end
%% visualizing flight duration vs replay speed
for hide = 1
    %=== for all replays
    % figure; scatter (abs(RP_velocity), RP_flight_duration,'filled');xlabel('replay speed of all replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
    %=== for top replays (high weighted corr)
    corr_cutoff = 0.9;
    % ids = find(abs(RP.corr)>corr_cutoff);
    ids = find(top_replay);
    % figure;scatter (abs(RP_velocity(ids)), RP_flight_duration(ids),'filled');xlabel('replay speed of high weighted corr replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
    %=== mean replay speeds for all replays
    flight_dur_unique = unique(RP_flight_duration);
    RP_velocity_mean = [];
    for i = 1:length(flight_dur_unique)
        current_flights = find(RP_flight_duration == flight_dur_unique(i));
        RP_velocity_mean = [RP_velocity_mean mean(abs(RP_velocity(current_flights)))];
    end
    % figure; scatter (RP_velocity_mean, flight_dur_unique,'filled');xlabel('mean replay speed of all replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
    %=== mean replay speed for top replays
    flight_dur_unique = unique(RP_flight_duration);
    RP_velocity_mean_top = [];
    for i = 1:length(flight_dur_unique)
        current_flights = find(RP_flight_duration == flight_dur_unique(i) & top_replay);
        RP_velocity_mean_top = [RP_velocity_mean_top mean(abs(RP_velocity(current_flights)))];
    end
    % figure;hold on; scatter (RP_velocity_mean_top, flight_dur_unique,'filled');xlabel('mean replay speed of top replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
    %=== visualizing everything
    figure;
    subplot(2,2,1);
    scatter (abs(RP_velocity), RP_flight_duration,'filled');xlabel('replay speed of all replays (m/s)'); ylabel('flight duration (s)');
    alpha(0.2);
    subplot(2,2,2);
    hold on; scatter (RP_velocity_mean, flight_dur_unique,'filled');xlabel('mean replay speed of all replays (m/s)'); ylabel('flight duration (s)');
    alpha(0.2);
    subplot(2,2,3);
    hold on; scatter (abs(RP_velocity(ids)), RP_flight_duration(ids),'filled');xlabel('replay replay speed of top replays (m/s)'); ylabel('flight duration (s)');
    alpha(0.2);
    subplot(2,2,4);
    hold on; scatter (RP_velocity_mean_top, flight_dur_unique,'filled');xlabel('mean replay speed of top replays (m/s)'); ylabel('flight duration (s)');
    alpha(0.2);
end
%% flight length vs replay speed (color by bat)
for hide = 1
    %=== calculations
    corr_cutoff = 0.9;
    top_ids = find(top_replay);
    bats = unique(RP.bat_id);
    
    RP_v_mean = [];
    RP_v_mean_top = [];
    flight_len_unique = unique(RP_flight_length);
    bat_ids_flight_dur = [];
    for i = 1:length(flight_len_unique)
        current_flights = find(RP_flight_length == flight_len_unique(i));
        current_bat = unique(RP.bat_id(current_flights));
        bat_ids_flight_dur = vertcat(bat_ids_flight_dur, str2num(current_bat{1}));
        RP_v_mean = [RP_v_mean nanmean(abs(RP_velocity(current_flights)))];
        current_flights = find(RP_flight_length == flight_len_unique(i) & top_replay);
        RP_v_mean_top = [RP_v_mean_top nanmean(abs(RP_velocity(current_flights)))];
    end
    
    
    %=== visualizing everything
      figure('units','normalized','outerposition',[0.3 .5 0.25 0.5]) ;
%     figure('units','normalized','outerposition',[0.3 .5 0.5 0.5]) ;
%     subplot(1,2,1);
%     for i = 1:numel(bats)
%         bat = bats(i);
%         bat_index = find(str2num(bat{1}) == RP_bat_id)';
%         hold on; scatter (abs(RP_velocity(bat_index)), RP_flight_length(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('Replay velocity (m/s)'); ylabel('flight length (m)');alpha(0.2);
%     end
%     legend('Location', 'southeast');
%     
%     subplot(1,2,2);
    for i = 1:numel(bats)
        bat = bats(i);
        bat_index = find(str2num(bat{1}) == bat_ids_flight_dur)';
%         hold on; scatter (RP_v_mean(bat_index), flight_len_unique(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('Mean replay velocity of cluster (m/s)'); ylabel('flight length (m)');
        hold on; scatter (RP_v_mean_top(bat_index), flight_len_unique(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('Mean replay velocity of cluster (m/s)'); ylabel('flight length (m)');
%         alpha(0.3);
    end
%     axis square;
%     lsline;
    p = polyfit(RP_v_mean, flight_len_unique,1);
    x = [0: max(RP_v_mean)+1];
    ylim([0 20])
    plot(x,polyval(p,x));
%     lsline();
%     b_reg = mean(flight_len_unique./RP_v_mean);
%     y_reg = b_reg * RP_v_mean;
%     hold on; plot(RP_v_mean, y_reg);
    % legend('Location', 'southeast');
%     
%     subplot(2,2,3);
%     
%     for i = 1:numel(bats)
%         bat = bats(i);
%         bat_index = intersect(find(str2num(bat{1}) == RP_bat_id)', top_ids);
%         hold on; scatter (abs(RP_velocity(bat_index)), RP_flight_length(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('temporal compression of top replays'); ylabel('flight duration (s)');alpha(0.2);
%     end
%     % legend('Location', 'southeast');
%     % scatter (abs(RP_comp(ids)), RP_flight_duration(ids),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('temporal compression of top replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
%     
%     subplot(2,2,4);
%     for i = 1:numel(bats)
%         bat = bats(i);
%         bat_index = find(str2num(bat{1}) == bat_ids_flight_dur)';
%         hold on; scatter (RP_v_mean_top(bat_index), flight_len_unique(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('mean temporal compression of top replays'); ylabel('flight duration (s)');alpha(0.3);
%     end
    % legend('Location', 'southeast');
    % scatter (RP_comp_mean_top, flight_dur_unique,'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('mean temporal compression of top replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
end
%% visualizing flight duration vs temporal compression (color by bat)
for hide = 1
    %=== calculations
    corr_cutoff = 0.9;
    top_ids = find(top_replay);
    bats = unique(RP.bat_id);
    
    
    
    RP_comp_mean = [];
    RP_comp_mean_top = [];
    flight_dur_unique = unique(RP_flight_duration);
    bat_ids_flight_dur = [];
    for i = 1:length(flight_dur_unique)
        current_flights = find(RP_flight_duration == flight_dur_unique(i));
        current_bat = unique(RP.bat_id(current_flights));
        bat_ids_flight_dur = vertcat(bat_ids_flight_dur, str2num(current_bat{1}));
        RP_comp_mean = [RP_comp_mean mean(abs(RP_comp(current_flights)))];
        current_flights = find(RP_flight_duration == flight_dur_unique(i) & top_replay);
        RP_comp_mean_top = [RP_comp_mean_top mean(abs(RP_comp(current_flights)))];
    end
    
    
    %=== visualizing everything
    figure;
    subplot(2,2,1);
    for i = 1:numel(bats)
        bat = bats(i);
        bat_index = find(str2num(bat{1}) == RP_bat_id)';
        hold on; scatter (abs(RP_comp(bat_index)), RP_flight_duration(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('temporal compression of all replays'); ylabel('flight duration (s)');alpha(0.2);
    end
    legend('Location', 'southeast');
    
    subplot(2,2,2);
    for i = 1:numel(bats)
        bat = bats(i);
        bat_index = find(str2num(bat{1}) == bat_ids_flight_dur)';
        hold on; scatter (RP_comp_mean(bat_index), flight_dur_unique(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('mean temporal compression of all replays'); ylabel('flight duration (s)');alpha(0.3);
    end
    % legend('Location', 'southeast');
    
    subplot(2,2,3);
    
    for i = 1:numel(bats)
        bat = bats(i);
        bat_index = intersect(find(str2num(bat{1}) == RP_bat_id)', top_ids);
        hold on; scatter (abs(RP_comp(bat_index)), RP_flight_duration(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('temporal compression of top replays'); ylabel('flight duration (s)');alpha(0.2);
    end
    % legend('Location', 'southeast');
    % scatter (abs(RP_comp(ids)), RP_flight_duration(ids),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('temporal compression of top replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
    
    subplot(2,2,4);
    for i = 1:numel(bats)
        bat = bats(i);
        bat_index = find(str2num(bat{1}) == bat_ids_flight_dur)';
        hold on; scatter (RP_comp_mean_top(bat_index), flight_dur_unique(bat_index),'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('mean temporal compression of top replays'); ylabel('flight duration (s)');alpha(0.3);
    end
    % legend('Location', 'southeast');
    % scatter (RP_comp_mean_top, flight_dur_unique,'filled','DisplayName',['bat ', num2str(bat{1})]);xlabel('mean temporal compression of top replays (m/s)'); ylabel('flight duration (s)');alpha(0.2);
end
%% plotting top examples of replays 
for hide = []
    %===sorting top_idx by selected variable
    top_idx = find(top_replay)
    
    figure('units','normalized','outerposition',[0 .05 0.6 0.8]);
    plot_index = 1;
    for ii = 1:5*6
        i=top_idx(ii);
        if plot_index > 25
            break
        end
        subplot(5,6,plot_index);
        colormap(hot);
        current_posteriors = RP_post_trimmed {i};
        uniform_prob = 1/size(current_posteriors,2);
        [position_prob, decoded_position] = max(current_posteriors');
        current_posteriors(find(position_prob' < 3*uniform_prob),:) = NaN;

        imagesc(current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
   
        if plot_index ==1
            ylabel("Position");
            xlabel("Time (s)");
        end
        set(gca, 'YTick', []);
        set(gca, 'XTick', []);
        
        xlabel([num2str(RP_dur(i)) ' ms']);
        ylabel([num2str(round(RP.flight_length{i},1)) ' m'])
        axis xy;
        plot_index = plot_index +1;
    end
end
%% Visualizing replay location vs replay flight type
for hide = 1
%=== finding unique clusters
sessions = unique(RP.session_id);
flight_ids = [];
current_flight_id = 1;
current_cluster = 1;
current_session = RP.session_id{1};
session_ids_flight = [{current_session}];
cluster_idx = [1];
flight_ids_uni = [current_cluster];
for i = 1: numel(RP.session_id)
    if isequal (current_session, RP.session_id{i}) == 0
        current_flight_id = current_flight_id +1;
        current_session =  RP.session_id{i};
        current_cluster = 1;
        %         disp(RP.session_id {i});
        %         disp(RP.cluster_id{i});
        %         disp(current_flight_id);
        session_ids_flight = [session_ids_flight {current_session}];
        flight_ids_uni = [flight_ids_uni current_flight_id];
        cluster_idx = [cluster_idx i];
    elseif RP.cluster_id{i} ~= current_cluster
        current_flight_id = current_flight_id +1;
        current_cluster = RP.cluster_id{i};
        %         disp(RP.session_id {i});
        %         disp(RP.cluster_id{i});
        %         disp(current_flight_id);
        session_ids_flight = [session_ids_flight {current_session}];
        flight_ids_uni = [flight_ids_uni current_flight_id];
        cluster_idx = [cluster_idx i];
    end
    flight_ids = [flight_ids current_flight_id];
end
cluster_type =  RP.flight_type(cluster_idx); 
cluster_include = strcmp(cluster_type, "nonfeeder-nonfeeder")==0;

%=== visualizing flight type
% figure; histogram(categorical(RP.flight_type(find(top_replay))), "DisplayOrder", "descend"); ylabel ("Replay count"); title ("distribution of replay by flight type")
flight_types = [];
flight_types_uni = unique(cluster_type(cluster_include));
for i = unique(flight_ids)
    flight_types = [flight_types RP.flight_type(find(flight_ids == i,1))];
end
figure; histogram(categorical(flight_types), "DisplayOrder", "descend"); ylabel ("Flight count"); title ("Flight distribution")
% plotting flight vs replay types 

%=== normalizing percent of TOP replays by occupancy
% all_replay_types = unique(RP.replay_type);
all_replay_types = {'remote', 'takeoff','landing'};
replay_occ_norm_all_flights = [];
for i = unique(flight_ids)
    replay_idx = find(flight_ids == i & top_replay); % selecting top replays of that flight
    %     replay_idx = find(flight_ids == i); % selecting top replays of that flight
    replay_types =  RP.replay_type(replay_idx);
    replay_types_uni = unique(replay_types);
    replay_occ = RP.replay_pos_occ(replay_idx);
    replay_count = [];
    replay_occ_uni = [];
    %     % finding the count of each type of replay and occupancy of that type
    %     for ii = replay_types_uni
    %         replay_type_idx = find(strcmpi(replay_types,ii));
    %         replay_count = [replay_count numel(replay_type_idx)];
    %         replay_occ_uni = [replay_occ_uni replay_occ(replay_type_idx(1))]; % takes the first element, bc all of them have the same occupancy value
    %     end
    % finding the count of each type of replay and occupancy of that type
    for ii = all_replay_types
        if ~any(strcmp(ii, replay_types_uni))
            replay_count = [replay_count 0];
            replay_occ_uni = [replay_occ_uni 0]; % takes the first element, bc all of them have the same occupancy value
        else
            replay_type_idx = find(strcmpi(replay_types,ii));
            replay_count = [replay_count numel(replay_type_idx)];
            replay_occ_uni = [replay_occ_uni replay_occ(replay_type_idx(1))]; % takes the first element, bc all of them have the same occupancy value
        end
    end
    replay_clus_norm = (replay_count./replay_occ_uni);
    replay_clus_norm (isnan(replay_clus_norm)) =0;
    replay_clus_norm = replay_clus_norm./sum(replay_clus_norm);
    replay_clus_norm (isnan(replay_clus_norm)) =0;
    replay_occ_norm_all_flights(i,:) =replay_clus_norm;
%     figure; bar(reordercats(categorical(all_replay_types),all_replay_types), replay_clus_norm);ylabel("Fraction replay normalized by occupancy");
end

figure; imagesc(replay_occ_norm_all_flights); xlabel ("replay types"); ylabel("flight #");
set(gca, 'XTick', [1:numel(all_replay_types)], 'XTickLabel', all_replay_types);

%=== visualizing replay type of all replays
% figure; histogram(categorical(RP.replay_type(find(top_replay))), "DisplayOrder", "descend"); ylabel ("Replay count"); title ("Distribution of top replays by replay type")
replay_loc_distribution = sum(replay_occ_norm_all_flights,1)./sum(sum(replay_occ_norm_all_flights,1)).*100;
[replay_loc_distribution, order] = sort(replay_loc_distribution, 'descend')
all_replay_types_sorted = all_replay_types(order)
figure; bar(replay_loc_distribution);ylabel ("Replay %"); title ("Distribution of top replays by replay location")
set(gca, 'XTick', [1:numel(all_replay_types_sorted)], 'XTickLabel', all_replay_types_sorted);
% reorder rows such that they sort y by flight type
% [sortedCategories, sortOrder] = sort(flight_types);
% replay_occ_all_flight_sorted = replay_occ_norm_all_flights(sortOrder);
% figure; imagesc(replay_occ_norm_all_flights); xlabel ("replay types"); ylabel("flight type");
figure;
flight_types_uni = unique(flight_types);
new_i = find(sum(replay_occ_norm_all_flights,2)==1);
replay_occ_norm_all_flights = replay_occ_norm_all_flights(new_i,:);
flight_types = flight_types(new_i);
for i = 1:numel(flight_types_uni)
    subplot(4,1,i);
    current_flights =  find(strcmp(flight_types,flight_types_uni(i)));
    imagesc(replay_occ_norm_all_flights(current_flights,:)); ylabel(flight_types_uni(i));
    set(gca, 'XTick', [1:numel(all_replay_types)], 'XTickLabel', all_replay_types);
    if i == 1
        title("top replays")
    end
end

for hide = []
%=== group the landing groups into a single group and takeoff groups into a single group
% new_replay_types = {'remote', 'feeder ', 'takeoff', 'feeder takeoff landing ', 'landing'};
% all_replay_types = {'remote', 'feeder ', 'takeoff ','feeder takeoff ', 'feeder takeoff landing ', 'landing ', 'feeder landing ', 'takeoff landing '};
new_replay_types = {'remote', 'landing', 'takeoff'};
all_replay_types = new_replay_types;
replay_dist_short = zeros(numel(flight_types), numel(new_replay_types)); % short version of replay_occ_norm_all_flights
for i = 1:numel(all_replay_types)
    if any(strcmp(all_replay_types{i}, new_replay_types))
        new_idx = find(strcmp(all_replay_types{i}, new_replay_types));
        replay_dist_short(:,new_idx) = replay_occ_norm_all_flights(:,i);
    elseif contains(all_replay_types{i}, "landing")
        new_idx = find(strcmp("landing", new_replay_types));
        replay_dist_short(:,new_idx) =  replay_dist_short(:,new_idx) + replay_occ_norm_all_flights(:,i);
    elseif contains(all_replay_types{i}, "takeoff")
        new_idx = find(strcmp("takeoff", new_replay_types));
        replay_dist_short(:,new_idx) =  replay_dist_short(:,new_idx) + replay_occ_norm_all_flights(:,i);
    end
end
figure;
for i = 1:numel(flight_types_uni)
    subplot(4,1,i);
    current_flights =  find(strcmp(flight_types,flight_types_uni(i)));
    imagesc(replay_dist_short(current_flights,:)); ylabel(flight_types_uni(i));
    set(gca, 'XTick', [1:numel(new_replay_types)], 'XTickLabel', new_replay_types);
    if i == 1
        title("top replays")
    end
end

% average across columns
type_distribution_mean = [];
for i = 1:numel(flight_types_uni)
    if isempty(replay_dist_short(current_flights,:))
        continue;
    end
    current_flights =  find(strcmp(flight_types,flight_types_uni(i)));
    type_distribution_mean (i,:) = mean(replay_dist_short(current_flights,:));
end
figure;
imagesc(type_distribution_mean);    
set(gca, 'XTick', [1:numel(new_replay_types)], 'XTickLabel', new_replay_types);
set(gca, 'YTick', [1:numel(flight_types_uni)], 'YTickLabel', flight_types_uni);
colorbar;

% simplify further, combine "other feeder" and remote 
% replay_dist_types_cut1 = {'remote', 'takeoff', 'feeder takeoff landing ', 'landing'};
replay_dist_types_cut1 = {'remote', 'landing', 'takeoff'};
type_distribution_mean_cut1 = [type_distribution_mean(:,1)+type_distribution_mean(:,2) type_distribution_mean(:,3:5)];
% figure;
% imagesc(type_distribution_mean_cut1);    
% set(gca, 'XTick', [1:numel(replay_dist_types_cut1)], 'XTickLabel', replay_dist_types_cut1);
% set(gca, 'YTick', [1:numel(flight_types_uni)], 'YTickLabel', flight_types_uni);
% colorbar;

%=== simplify to remote, takeoff, landing
replay_dist_types_cut2 = {'remote', 'takeoff', 'landing'};
type_distribution_mean_cut2 = [type_distribution_mean_cut1(:,1:2) type_distribution_mean_cut1(:,4)];
type_distribution_mean_cut2(1,2) = (1-type_distribution_mean_cut2(1,1))/2;
type_distribution_mean_cut2(1,3) = (1-type_distribution_mean_cut2(1,1))/2;


figure;
imagesc(type_distribution_mean_cut2);    
set(gca, 'XTick', [1:numel(replay_dist_types_cut2)], 'XTickLabel', replay_dist_types_cut2);
set(gca, 'YTick', [1:numel(flight_types_uni)], 'YTickLabel', flight_types_uni);
colorbar;
end
end
%% distribution of flight
for hide = 1
%=== visualizing flight type
% figure; histogram(categorical(RP.flight_type(find(top_replay))), "DisplayOrder", "descend"); ylabel ("Replay count"); title ("distribution of replay by flight type")
flight_types = [];
cluster_type =  RP.flight_type(cluster_idx); 
cluster_include = strcmp(cluster_type, "nonfeeder-nonfeeder")==0;
flight_types_uni = unique(cluster_type(cluster_include));
for i = unique(flight_ids)
    flight_types = [flight_types RP.flight_type(find(flight_ids == i,1))];
end
figure; histogram(categorical(flight_types), "DisplayOrder", "descend"); ylabel ("Flight count"); title ("Flight distribution")
% count the distribution of flights
summary(categorical(flight_types))
end
%% visualizing flight type vs fragmentation
%=== by flight types 
for hide = []
flight_types_uni = unique(RP.flight_type);
figure('units','normalized','outerposition',[0.3 0.5 0.1 0.5]);
tiledlayout(numel(flight_types_uni),1);
for i = 1:numel(flight_types_uni)
    nexttile;
    current_flights =  find(strcmp(RP.flight_type,flight_types_uni(i)));
    histogram(segment_len_frac(current_flights),10,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    ylabel(flight_types_uni(i))
end
xlabel("segment length")


%=== by replay locations 
replay_types_uni = unique(RP.replay_type);
figure('units','normalized','outerposition',[0.5 0 0.1 1]);
tiledlayout(numel(replay_types_uni),1);
for i = 1:numel(replay_types_uni)
    nexttile;
    current_flights =  find(strcmp(RP.replay_type,replay_types_uni(i)));
    histogram(segment_len_frac(current_flights),10,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    ylabel(replay_types_uni(i))
end
xlabel("segment length")

%=== repaly trajectory length vs replay duration [like davidson]
figure;scatter (replay_dur_bin(top_idx), segment_len_m(top_idx),"filled");alpha(0.1);xlabel("Replay duration (s)");ylabel("Repaly trajectory length (m)");
xlim([0.05,0.25])
end
%% plot peri vs nonperi by session 
for hide = 1
% sessions = unique(RP.session_id);
current_session_id = 1;
current_cluster = 1;
current_session = RP.session_id{1};
session_ids_flight = [{current_session}]; % id of the current session 
session_idx = [1]; % first index of a session 
session_ids = [current_session_id]; %1: # of session
for i = 1: numel(RP.session_id) % loop through session 
    if isequal (current_session, RP.session_id{i}) == 0
        current_session_id = current_session_id +1;
        current_session =  RP.session_id{i};
        session_ids_flight = [session_ids_flight {current_session}];
        session_idx = [session_idx i];
        session_ids = [session_ids current_session_id];
        
    end
end
peri_rates_all = cell2mat(RP.peri_rate(session_idx));
nonperi_rates_all = cell2mat(RP.nonperi_rate(session_idx));

% figure; 
% PlotDistr_AF_v1([peri_rates_all; nonperi_rates_all],turbo(size(flighttype_comp,2)),"# Replays per second",["replay near flight" "replay far from flight"]);
[~,p] = ttest(peri_rates_all,nonperi_rates_all)
p
signrank(peri_rates_all,nonperi_rates_all)
% ylim([0,1])
end
% calculating replay rate vs replay location (for fig 3b)
for hide = 1
% quick visualization 
for hide = []
segment_min = 0;
max_len = 0;
location_types = unique(RP_RTL);
for i = 1:numel(location_types)
    max_len = max(max_len, numel(find (top_replay & strcmp(RP_RTL,location_types{i}) & segment_len_frac > segment_min)));
end
location_replay_rate = NaN(max_len,numel(location_types));
for i = 1:numel(location_types)
    curr_replay_rate = RP.replay_rate_by_type(find(top_replay & strcmp(RP_RTL,location_types{i})& segment_len_frac > segment_min));
    location_replay_rate((1:numel(curr_replay_rate)),i)=curr_replay_rate;
end
% figure;
% PlotDistr_AF_v1(location_replay_rate',turbo(size(location_replay_rate,2)),"Replay rate (#replays/second)",location_types);
end

% loop through each session and flight, find all the unique replay types and their corresponding rate. declare based on RLT
sessions = unique(RP.session_id);
flight_ids = [];                                      % same size as number of replay, shows the flight id
current_flight_id = 1;
current_cluster = 1;
current_session = RP.session_id{1};
session_ids_flight = [{current_session}];             % all the exp session ids 
cluster_idx = [1];                                    % the index of the first value of a cluster 
flight_ids_uni = [current_cluster];
for i = 1: numel(RP.session_id)
    if isequal (current_session, RP.session_id{i}) == 0
        current_flight_id = current_flight_id +1;
        current_session =  RP.session_id{i};
        current_cluster = 1;
        session_ids_flight = [session_ids_flight {current_session}];
        flight_ids_uni = [flight_ids_uni current_flight_id];
        cluster_idx = [cluster_idx i];
    elseif RP.cluster_id{i} ~= current_cluster % loop through cluster 
        current_flight_id = current_flight_id +1;
        current_cluster = RP.cluster_id{i};
        session_ids_flight = [session_ids_flight {current_session}];
        flight_ids_uni = [flight_ids_uni current_flight_id];
        cluster_idx = [cluster_idx i];
    end
    flight_ids = [flight_ids current_flight_id];
end
cluster_type =  RP.flight_type(cluster_idx); 
cluster_include = strcmp(cluster_type, "nonfeeder-nonfeeder")==0;

check = [];
% location_seq = ["R","T","L","TL"];
location_seq = ["R","T","L"];
replay_rate_location = zeros(numel(flight_ids_uni), numel(location_seq));
for i = 1:numel(cluster_idx)
    clus_st = cluster_idx(i);
    if i == numel(cluster_idx)
        clus_en = numel(RP.session_id);
    else 
        clus_en = cluster_idx (i+1)-1;
    end
    current_t_vec = [clus_st:clus_en];
    [category, cat_idx, test]=unique(RP.replay_type(current_t_vec));
%     [category, index, ~]=unique(RP_RTL([clus_st:clus_en]));

    for ii = 1:numel(category)
%         if contains(category(ii), "takeoff landing")
%             replay_rate_location(i,find(location_seq=="TL"))=replay_rate_location(i,find(location_seq=="TL"))+RP.replay_rate_by_type(current_t_vec(cat_idx(ii)));
        if contains(category(ii), "takeoff")
            replay_rate_location(i,find(location_seq=="T"))=replay_rate_location(i,find(location_seq=="T"))+RP.replay_rate_by_type(current_t_vec(cat_idx(ii)));
        elseif contains(category(ii), "landing")
            replay_rate_location(i,find(location_seq=="L"))=replay_rate_location(i,find(location_seq=="L"))+RP.replay_rate_by_type(current_t_vec(cat_idx(ii)));
        elseif contains(category(ii), "feeder")
        replay_rate_location(i,find(location_seq=="R"))=replay_rate_location(i,find(location_seq=="R"))+RP.replay_rate_by_type(current_t_vec(cat_idx(ii)));
        elseif contains(category(ii), "remote")
        replay_rate_location(i,find(location_seq=="R"))=replay_rate_location(i,find(location_seq=="R"))+RP.replay_rate_by_type(current_t_vec(cat_idx(ii)));
        end
    end
    check = [check; numel(category)];
end
end
%% replay rate vs location (from session-level calculations)
for hide = []
cluster_include = ~strcmp(flight_types,'nonfeeder-nonfeeder');
%=== box plot 
figure;
replay_rate_location_i = replay_rate_location(cluster_include,:);
replay_rate_non_zero = replay_rate_location_i;
replay_rate_non_zero(replay_rate_non_zero==0) = nan
PlotDistr_AF_v1(replay_rate_non_zero',turbo(size(replay_rate_non_zero,2)),"Replay rate (#replays/second)",location_seq);

%=== calculating significance 
% perform friedman test 
[p,tbl,stats]=friedman(replay_rate_location_i,1)
mtp_comparisons = multcompare(stats,[],'on');
mtp_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);



%=== line plot 
figure('units','normalized','outerposition',[0.3 .2 0.15 0.4])
v1 = replay_rate_location_i(:,1);
v2 = replay_rate_location_i(:,2);
v3 = replay_rate_location_i(:,3);

for i = 1: numel(v1)
    line([1,2],[v1(i),v2(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
end
line([1,2],[mean(v1),mean(v2)], 'Color', [0 0 0], 'Linewidth',4);  hold on;
for i = 1: numel(v2)
    line([2,3],[v2(i),v3(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
end
line([2,3],[mean(v2),mean(v3)], 'Color', [0 0 0], 'Linewidth',4);  hold on;


v1_sem = std(v1) / sqrt(numel(v1));
v2_sem =  std(v2) / sqrt(numel(v2));
v3_sem =  std(v3) / sqrt(numel(v3));
line([1,1],[mean(v1)+v1_sem,mean(v1)-v1_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
line([2,2],[mean(v2)+v2_sem,mean(v2)-v2_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
line([3,3],[mean(v3)+v3_sem,mean(v3)-v3_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
xlim([0.5, 3.5]);
set(gca, 'XTick', []);
ylabel("#Replays/s");
xticks([1,2,3]);
xticklabels(location_seq);
end
%% calculating replay rate on the pooling level 
for hide = []
% loop through flight cluster 
cluster_idx_all = [cluster_idx numel(RP.replay_type)];
replay_types  = [{'remote'} {'takeoff'} {'landing'}];
clusters_occ = [];                                                          % time occupancy of replays 
clusters_replay_num = [];
distance_threshold = 0.5;
replay_rates =[];                                                           % replay rates by cluster {cluster} ('landing','remote','takeoff')
parameter = best_replay;
for i = 1: numel(cluster_idx)
    disp(i)
    if RP.cluster_id{cluster_idx(i)} == 1 || strcmp(RP.flight_type(cluster_idx(i)),'nonfeeder-nonfeeder')
        replay_rates (i,:) = [NaN NaN NaN];
        continue;
    end
    current_replay_idx = intersect(parameter, cluster_idx_all(i): cluster_idx_all (i+1)); % same cluster, and top replays 
    cellArray = sessions_info.curr_session_id;
    searchString = RP.session_id(cluster_idx(i));
    session_idx = find (cellfun(@(x) iscell(x) && ~isempty(x) && ischar(x{1}) && strcmp(x{1}, searchString), cellArray));
    r = sessions_info.r{session_idx};
    for ii = 1:numel(replay_types)
        current_type =  current_replay_idx(find(strcmp(RP.replay_type(current_replay_idx), {replay_types{ii}})));
        pos_idx=[];
        for iii = 1: numel(current_type) % mark all the positions near the replay for a specific cluster
            pos_idx = union (pos_idx, find(vecnorm(r{1}-RP.replay_position(:,current_type(iii))', 2,2)<distance_threshold));
        end
        Fs = cell2mat(RP.Fs_session(cluster_idx(i)));
        clusters_occ([ii]) = numel(pos_idx) /Fs;
        clusters_replay_num([ii]) = numel(current_type);
    end
    replay_rates (i,:) = clusters_replay_num./clusters_occ;
end
replay_rates
%% all replay
figure;
PlotDistr_AF_v1(replay_rates',turbo(size(replay_rates,2)),"Replay rate (#replays/second)",replay_types);
replay_rates_nonan = replay_rates(sum(~isnan(replay_rates),2)>1,:);
replay_rates_nonan (isnan(replay_rates_nonan)) = 0;
[p,tbl,stats]=friedman(replay_rates_nonan,1)
mtp_comparisons = multcompare(stats,[],'on');
mtp_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
%% separate by cluster type 
figure;
FN_cluster_idx = strcmp(cluster_type, "feeder-nonfeeder")
FN_replay_rates = replay_rates(FN_cluster_idx, :);
PlotDistr_AF_v1(FN_replay_rates',turbo(numel(FN_cluster_idx)),"Feeder-nonfeeder replay rate (#replays/second)",replay_types);
replay_rates_nonan = FN_replay_rates(sum(~isnan(FN_replay_rates),2)>1,:);
replay_rates_nonan (isnan(replay_rates_nonan)) = 0;
[p,tbl,stats]=friedman(replay_rates_nonan,1)
mtp_comparisons = multcompare(stats,[],'on');
mtp_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

figure;
NF_cluster_idx = strcmp(cluster_type, "nonfeeder-feeder")
NF_replay_rates = replay_rates(NF_cluster_idx, :);
PlotDistr_AF_v1(NF_replay_rates',turbo(numel(FN_cluster_idx)),"Nonfeeder-feeder replay rate (#replays/second)",replay_types);
replay_rates_nonan = NF_replay_rates(sum(~isnan(NF_replay_rates),2)>1,:);
replay_rates_nonan (isnan(replay_rates_nonan)) = 0;
[p,tbl,stats]=friedman(replay_rates_nonan,1)
mtp_comparisons = multcompare(stats,[],'on');
mtp_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

figure;
FN_NF_cluster_idx = FN_cluster_idx|NF_cluster_idx;
FN_NF_replay_rates = replay_rates(FN_NF_cluster_idx, :);
PlotDistr_AF_v1(FN_NF_replay_rates',turbo(numel(FN_NF_cluster_idx)),"Nonfeeder-feeder & feeder-nonfeedr replay rate (#replays/second)",replay_types);
replay_rates_nonan = FN_NF_replay_rates(sum(~isnan(FN_NF_replay_rates),2)>1,:);
replay_rates_nonan (isnan(replay_rates_nonan)) = 0;
[p,tbl,stats]=friedman(replay_rates_nonan,1)
mtp_comparisons = multcompare(stats,[],'on');
mtp_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
end
%% plot precalculated peri vs nonperi by cluster
y_lim_val = 0.1;
for hide = []% sessions = unique(RP.session_id);
curr_cluster_idx = cluster_idx(~strcmp(RP.flight_type(cluster_idx),  {'nonfeeder-nonfeeder'}));
peri_rates_all_cluster  = cell2mat(RP.peri_rate_cluster(curr_cluster_idx));
nonperi_rates_all_cluster = cell2mat(RP.nonperi_rate_cluster(curr_cluster_idx));

%=== boxplot 
% figure; 
% PlotDistr_AF_v1([peri_rates_all_cluster; nonperi_rates_all_cluster],turbo(size(flighttype_comp,2)),"# Replays per second",["replay near flight" "replay far from flight"]);
[~,p] = ttest(peri_rates_all_cluster,nonperi_rates_all_cluster);
ps = signrank(peri_rates_all_cluster,nonperi_rates_all_cluster);
disp(['ttest p: ' num2str(p) '   signrank p: ' num2str(ps)])
%=== line plot 
figure('units','normalized','outerposition',[0.3 .4 0.13 0.5])
ylim ([0 y_lim_val])
for i = 1: numel(peri_rates_all_cluster)
%     line([1,2],[peri_rates_all_cluster(i),nonperi_rates_all_cluster(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
end
line([1,2],[mean(peri_rates_all_cluster),mean(nonperi_rates_all_cluster)], 'Color', [0 0 0], 'Linewidth',4);  hold on;

peri_sem = std(peri_rates_all_cluster) / sqrt(numel(peri_rates_all_cluster));
nonperi_sem =  std(nonperi_rates_all_cluster) / sqrt(numel(nonperi_rates_all_cluster));
line([1,1],[mean(peri_rates_all_cluster)+peri_sem,mean(peri_rates_all_cluster)-peri_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
line([2,2],[mean(nonperi_rates_all_cluster)+nonperi_sem,mean(nonperi_rates_all_cluster)-nonperi_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
disp (['peri mean: ' num2str(mean(peri_rates_all_cluster)) '   peri sem: ' num2str(peri_sem)])
disp (['nonperi mean: ' num2str(mean(nonperi_rates_all_cluster)) '   peri sem: ' num2str(nonperi_sem)])
xlim([0.9, 2.1]);
ylim([0 y_lim_val]);
set(gca, 'XTick', []);
ylabel("#Replays/s");
xticks([1,2]);
xticklabels({"replay near flight", "replay far from flight"});

%% peri vs nonperi by cluster, seperated by flight type 

% figure; 
% PlotDistr_AF_v1([peri_rates_all_cluster; nonperi_rates_all_cluster],turbo(size(flighttype_comp,2)),"# Replays per second",["replay near flight" "replay far from flight"]);
% [~,p] = ttest(peri_rates_all_cluster,nonperi_rates_all_cluster)
% signrank(peri_rates_all_cluster,nonperi_rates_all_cluster)

%=== 3 plots, each one is a different flight type
% unique_flight_types = unique(RP.flight_type);
% unique_flight_types = unique_flight_types(~strcmp(unique(RP.flight_type), "nonfeeder-nonfeeder"));
unique_flight_types = [   {'feeder-nonfeeder'}  {'nonfeeder-feeder'} {'feeder-feeder'}   ];
max_len = 0;
replay_rate_pnp_ttest_p = [];
replay_rate_pnp_srank_p = [];
for clus = 1:numel(unique_flight_types)
    figure('units','normalized','outerposition',[0.3 .4 0.13 0.5])
    curr_cluster_idx = cluster_idx(strcmp(RP.flight_type(cluster_idx),unique_flight_types(clus)));
    max_len = max(max_len, numel(curr_cluster_idx));
    
    peri_rates_curr_cluster  = cell2mat(RP.peri_rate_cluster(curr_cluster_idx));
    nonperi_rates_curr_cluster = cell2mat(RP.nonperi_rate_cluster(curr_cluster_idx));
    %=== VISUALIZATION
    %line visualization 
%     PlotDistr_AF_v1([peri_rates_curr_cluster; nonperi_rates_curr_cluster],turbo(size(flighttype_comp,2)),"# Replays per second",["replay near flight" "replay far from flight"]);
    %=== line plot 
    v1 = peri_rates_curr_cluster;
    v2 = nonperi_rates_curr_cluster;

    for i = 1: numel(v1)
        line([1,2],[v1(i),v2(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
    end
    line([1,2],[mean(v1),mean(v2)], 'Color', [0 0 0], 'Linewidth',4);  hold on;
    v1_sem = std(v1) / sqrt(numel(v1));
    v2_sem =  std(v2) / sqrt(numel(v2));
    line([1,1],[mean(v1)+v1_sem,mean(v1)-v1_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
    line([2,2],[mean(v2)+v2_sem,mean(v2)-v2_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
    xlim([0.9, 2.1]);
    set(gca, 'XTick', []);
    ylabel("#Replays/s");
    xticks([1,2]);
    xticklabels({"near","far"});
    sgtitle(unique_flight_types(clus));
    ylim([0,y_lim_val])
    %significance tests 
    [~,p] = ttest(peri_rates_curr_cluster,nonperi_rates_curr_cluster);
    p_sr=signrank(peri_rates_curr_cluster,nonperi_rates_curr_cluster);
    replay_rate_pnp_ttest_p = [replay_rate_pnp_ttest_p,p];
    replay_rate_pnp_srank_p = [replay_rate_pnp_srank_p,p_sr];
end

%% one plot, separated by flight type and by peri-nonperi
% unique_flight_types = unique(RP.flight_type);
% unique_flight_types = unique_flight_types(~strcmp(unique(RP.flight_type), "nonfeeder-nonfeeder"));
for hide = []
replay_rate_pnp_cluster_all = nan(numel(unique_flight_types)*2, max_len);
x_labels_pnp = [];
curr_group = 1;
for clus = 1:numel(unique_flight_types)

    curr_cluster_idx = cluster_idx(strcmp(RP.flight_type(cluster_idx),unique_flight_types(clus)));
    
    peri_rates_curr_cluster  = cell2mat(RP.peri_rate_cluster(curr_cluster_idx));
    nonperi_rates_curr_cluster = cell2mat(RP.nonperi_rate_cluster(curr_cluster_idx));

    replay_rate_pnp_cluster_all(curr_group,1:numel(curr_cluster_idx))=peri_rates_curr_cluster;
    replay_rate_pnp_cluster_all(curr_group+1,1:numel(curr_cluster_idx))=nonperi_rates_curr_cluster;
    cluster_type_curr = unique_flight_types(clus);
    x_labels_pnp = [x_labels_pnp {[cluster_type_curr{1} ' P']} {[cluster_type_curr{1} ' NP']}]; 
    curr_group = curr_group+2;
end

% visualization 
figure('units','normalized','outerposition',[0.3 .2 0.5 0.7])
PlotDistr_AF_v1(replay_rate_pnp_cluster_all,turbo(size(flighttype_comp,2)),"# Replays per second",x_labels_pnp);
%     [~,p] = ttest(peri_rates_curr_cluster,nonperi_rates_curr_cluster)
%     signrank(peri_rates_curr_cluster,nonperi_rates_curr_cluster)
%     sgtitle(unique_flight_types(clus));
ylim([0,0.4])

% statistical test 
% peri only
[~,~,stats] = kruskalwallis(replay_rate_pnp_cluster_all([1,3,5],:)');
mtp_comparisons = multcompare(stats,[],'on');
fly_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

[~,~,stats] = kruskalwallis(replay_rate_pnp_cluster_all([2,4,6],:)');
mtp_comparisons = multcompare(stats,[],'on');
fly_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
end
end
%% CALCULATING  peri-flights and non-peri flights 
for hide = 1% loop through each session 
for hide = []
replays_used = good_replay;
session_ids_unique = unique(RP.session_id);
peri_rate_all = [];
nonperi_rate_all = [];
for session = 1: numel(session_ids_unique)
    current_session_str = session_ids_unique(session);
    % finding replays within 30s of flight
    current_session_replay = replays_used & strcmp(RP.session_id , current_session_str);
    peri_replay = current_session_replay & (RP.replay_flight_next_dist_any >=-30 | RP.replay_flight_prev_dist_any <=30);
    nonperi_replay = current_session_replay & (RP.replay_flight_next_dist_any <-30 & RP.replay_flight_prev_dist_any>30);
    % finding the amount of time that peri/non-peri occupies
    cellArray = sessions_info.curr_session_id;
    searchString = current_session_str;
    session_idx = find (cellfun(@(x) iscell(x) && ~isempty(x) && ischar(x{1}) && strcmp(x{1}, searchString), cellArray),1,'first');
    peri_dur = cell2mat(sessions_info.peri_dur{session_idx});
    nonperi_dur = cell2mat(sessions_info.nonperi_dur{session_idx});
    peri_rate_all = [peri_rate_all sum(peri_replay)/peri_dur];
    nonperi_rate_all = [nonperi_rate_all sum(nonperi_replay)/nonperi_dur];
end

[~,p] = ttest(peri_rate_all,nonperi_rate_all);
ps = signrank(peri_rate_all,nonperi_rate_all);
disp(['ttest p: ' num2str(p) '   signrank p: ' num2str(ps)])
%=== line plot 
figure('units','normalized','outerposition',[0.3 .4 0.13 0.5])
ylim ([0 y_lim_val])
for i = 1: numel(peri_rate_all)
    line([1,2],[peri_rate_all(i),nonperi_rate_all(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
end
line([1,2],[mean(peri_rate_all),mean(nonperi_rate_all)], 'Color', [0 0 0], 'Linewidth',4);  hold on;

peri_sem = std(peri_rate_all) / sqrt(numel(peri_rate_all));
nonperi_sem =  std(nonperi_rate_all) / sqrt(numel(nonperi_rate_all));
line([1,1],[mean(peri_rate_all)+peri_sem,mean(peri_rate_all)-peri_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
line([2,2],[mean(nonperi_rate_all)+nonperi_sem,mean(nonperi_rate_all)-nonperi_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
disp (['peri mean: ' num2str(mean(peri_rate_all)) '   peri sem: ' num2str(peri_sem)])
disp (['nonperi mean: ' num2str(mean(nonperi_rate_all)) '   peri sem: ' num2str(nonperi_sem)])
xlim([0.9, 2.1]);
ylim([0 y_lim_val]);
set(gca, 'XTick', []);
ylabel("#Replays/s");
xticks([1,2]);
xticklabels({"replay near flight", "replay far from flight"});
end
%  loop through each cluster
replays_used = good_replay;
session_ids_unique = unique(RP.session_id);
peri_rate_by_cluster = [];
nonperi_rate_by_cluster = [];
peri_sum = 0;
nonperi_sum = 0;
% loop through each cluster 
cluster_idx_all = [cluster_idx numel(RP.replay_type)];
for i = 1: numel(cluster_idx)
    if RP.cluster_id{cluster_idx(i)} == 1 || strcmp(RP.flight_type(cluster_idx(i)),'nonfeeder-nonfeeder')
        continue;
    end
%     current_local = sum(strcmp (RP_local_nonlocal (cluster_idx_all(i): cluster_idx_all (i+1)),'local'));
%     current_nonlocal = sum(strcmp (RP_local_nonlocal (cluster_idx_all(i): cluster_idx_all (i+1)),'nonlocal'));
    current_session_str = RP.session_id(cluster_idx(i));
    % finding replays within 30s of flight
    current_session_replay = replays_used & strcmp(RP.session_id , current_session_str);
    current_cluster_replay = intersect(find(current_session_replay), [(cluster_idx_all(i): cluster_idx_all (i+1))]);
    peri_replay = intersect(current_cluster_replay, find(RP.replay_flight_next_dist_any >=-30 | RP.replay_flight_prev_dist_any <=30));
    nonperi_replay = intersect(current_cluster_replay, find(RP.replay_flight_next_dist_any <-30 & RP.replay_flight_prev_dist_any>30));
    % finding the amount of time that peri/non-peri occupies
    cellArray = sessions_info.curr_session_id;
    searchString = current_session_str;
    session_idx = find (cellfun(@(x) iscell(x) && ~isempty(x) && ischar(x{1}) && strcmp(x{1}, searchString), cellArray),1,'first');
    peri_dur = cell2mat(sessions_info.peri_dur{session_idx});
    nonperi_dur = cell2mat(sessions_info.nonperi_dur{session_idx});
    peri_sum = peri_sum + numel(peri_replay);
    nonperi_sum = nonperi_sum + numel(nonperi_replay);
    peri_rate_by_cluster = [peri_rate_by_cluster numel(peri_replay)/peri_dur];
    nonperi_rate_by_cluster = [nonperi_rate_by_cluster numel(nonperi_replay)/nonperi_dur];
end

[~,p] = ttest(peri_rate_by_cluster,nonperi_rate_by_cluster);
ps = signrank(peri_rate_by_cluster,nonperi_rate_by_cluster);
disp(['ttest p: ' num2str(p) '   signrank p: ' num2str(ps)])
%=== line plot 
figure('units','normalized','outerposition',[0.3 .4 0.13 0.5])
ylim ([0 y_lim_val])
for i = 1: numel(peri_rate_by_cluster)
    line([1,2],[peri_rate_by_cluster(i),nonperi_rate_by_cluster(i)], 'Color', [0.8 0.8 0.8 ],'Linewidth',1);  hold on;
end
line([1,2],[mean(peri_rate_by_cluster),mean(nonperi_rate_by_cluster)], 'Color', [0 0 0], 'Linewidth',4);  hold on;

peri_sem = std(peri_rate_by_cluster) / sqrt(numel(peri_rate_by_cluster));
nonperi_sem =  std(nonperi_rate_by_cluster) / sqrt(numel(nonperi_rate_by_cluster));
line([1,1],[mean(peri_rate_by_cluster)+peri_sem,mean(peri_rate_by_cluster)-peri_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
line([2,2],[mean(nonperi_rate_by_cluster)+nonperi_sem,mean(nonperi_rate_by_cluster)-nonperi_sem], 'Color', [0 0 0], 'Linewidth',2);  hold on;
disp (['peri mean: ' num2str(mean(peri_rate_by_cluster)) '   peri sem: ' num2str(peri_sem)])
disp (['nonperi mean: ' num2str(mean(nonperi_rate_by_cluster)) '   peri sem: ' num2str(nonperi_sem)])
xlim([0.9, 2.1]);
ylim([0 0.015]);
set(gca, 'XTick', []);
ylabel("#Replays/s");
xticks([1,2]);
xticklabels({"replay near flight", "replay far from flight"});
disp(['numel of clusters: ', num2str(numel(peri_rate_by_cluster))])
disp(['percent of peri flights', num2str(peri_sum / (peri_sum + nonperi_sum))]);
end
%% total replay rate by 3 flight type
for hide = []
replay_rate_total_cluster_all = nan(numel(unique_flight_types), max_len);
for clus = 1:numel(unique_flight_types)
    curr_cluster_idx = cluster_idx(strcmp(RP.flight_type(cluster_idx),unique_flight_types(clus)));
    total_rates_curr_cluster  = cell2mat(RP.total_rate_cluster(curr_cluster_idx));
    replay_rate_total_cluster_all(clus,1:numel(curr_cluster_idx))=total_rates_curr_cluster;
end

figure('units','normalized','outerposition',[0.3 .2 0.25 0.5])
PlotDistr_AF_v1(replay_rate_total_cluster_all,turbo(size(replay_rate_total_cluster_all,2)),"# Replays per second",unique_flight_types);


v1 = replay_rate_total_cluster_all(1,:);
v2 = replay_rate_total_cluster_all(2,:);
v3 = replay_rate_total_cluster_all(3,:);
v1 = v1(~isnan(v1));
v2 = v2(~isnan(v2));
v3 = v3(~isnan(v3));


v1_sem = std(v1) / sqrt(numel(v1));
v2_sem =  std(v2) / sqrt(numel(v2));
v3_sem =  std(v3) / sqrt(numel(v3));

disp(['v1 mean: ' num2str(mean(v1)) '    v1 sem' num2str(v1_sem)])
disp(['v2 mean: ' num2str(mean(v2)) '    v2 sem' num2str(v2_sem)])
disp(['v3 mean: ' num2str(mean(v3)) '    v3 sem' num2str(v3_sem)])

[~,~,stats] = kruskalwallis(replay_rate_total_cluster_all');
mtp_comparisons = multcompare(stats,[],'on');
fly_tbl = array2table(mtp_comparisons,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
end
%% total replay rate by 2 flight type (merged)
for hide = []
RP_flight_to_feeder = ~strcmp(RP.flight_type,'feeder-nonfeeder');

cluster_idx_FN = cluster_idx(strcmp(RP.flight_type(cluster_idx),'feeder-nonfeeder'));
cluster_idx_NF = cluster_idx(strcmp(RP.flight_type(cluster_idx),'nonfeeder-feeder'));
cluster_idx_FF = cluster_idx(strcmp(RP.flight_type(cluster_idx),'feeder-feeder'));

cluster_idx_to_feeder =  cluster_idx(RP_flight_to_feeder(cluster_idx));
cluster_idx_to_nonfeeder =  cluster_idx(~RP_flight_to_feeder(cluster_idx));

cluster_idx_to_feeder = [cluster_idx_NF cluster_idx_FF];
cluster_idx_to_nonfeeder = cluster_idx_FN;


replay_rate_total_cluster_all = nan(2, max(numel(cluster_idx_to_feeder),numel(cluster_idx_to_nonfeeder)));

total_rates_to_feeder  = cell2mat(RP.total_rate_cluster(cluster_idx_to_feeder));
replay_rate_total_cluster_all(1,1:numel(total_rates_to_feeder))=total_rates_to_feeder;

total_rates_to_nonfeeder  = cell2mat(RP.total_rate_cluster(cluster_idx_to_nonfeeder));
replay_rate_total_cluster_all(2,1:numel(total_rates_to_nonfeeder))=total_rates_to_nonfeeder;

figure('units','normalized','outerposition',[0.3 .2 0.25 0.5])
PlotDistr_AF_v1(replay_rate_total_cluster_all,turbo(size(replay_rate_total_cluster_all,2)),"# Replays per second",{'to feeder', 'to nonfeeder'});

[~,p] = ttest(replay_rate_total_cluster_all(1,:),replay_rate_total_cluster_all(2,:))
end
%% FIGURE 2B total replay rate by 2 flight type 
for hide = []
RP_flight_to_feeder = ~strcmp(RP.flight_type,'feeder-nonfeeder');

cluster_idx_FN = cluster_idx(strcmp(RP.flight_type(cluster_idx),'feeder-nonfeeder'));
cluster_idx_NF = cluster_idx(strcmp(RP.flight_type(cluster_idx),'nonfeeder-feeder'));

replay_rate_total_cluster_all = nan(2, max(numel(cluster_idx_FN),numel(cluster_idx_NF)));

total_rates_to_feeder  = cell2mat(RP.total_rate_cluster(cluster_idx_NF));
replay_rate_total_cluster_all(1,1:numel(total_rates_to_feeder))=total_rates_to_feeder;

total_rates_to_nonfeeder  = cell2mat(RP.total_rate_cluster(cluster_idx_FN));
replay_rate_total_cluster_all(2,1:numel(total_rates_to_nonfeeder))=total_rates_to_nonfeeder;

figure('units','normalized','outerposition',[0.3 .2 0.25 0.5])
PlotDistr_AF_v1(replay_rate_total_cluster_all,turbo(size(replay_rate_total_cluster_all,2)),"# Replays per second",{'nonfeeder-feeder', 'feeder-nonfeeder'});

[~,p] = ttest(replay_rate_total_cluster_all(1,:),replay_rate_total_cluster_all(2,:))
end
%% flight type clustering/sigmoid fit, speed profile 
for hide = []
%% flight type clustering
for hide = 1
num_clus = numel(cluster_idx);
%=== loop across cluster 
med_curv = NaN(num_clus,1);
gaussR2 = NaN(num_clus,1);
for ii = 1: num_clus

    %=== Get mean path in 3D
    i = cluster_idx(ii);
    if RP.cluster_id{i}==1 % exclude the unclustered flights
        continue
    end
    mean_path3D = RP.mean_path3D{i}';  % Assume this is the mean position of the bat in 3D (sampling does not matter, curvature does not depend on that)
    
    %=== Calculate curvature
    norm_grad3D = diff(mean_path3D,1,1)./vecnorm(diff(mean_path3D,1,1),2,2);
    curv3D = vecnorm(diff(norm_grad3D,1,1),2,2);
    curv3D (isnan(curv3D)) = [];
    med_curv(ii) = median(curv3D);
    
    
    %=== Get mean paths in 3D
    mean_path1D =  RP.mean_path1D{i};  % Assume this is the mean lenght of the trajectory at each time sample
    mean_time1D = RP.mean_time1D{i};  % Assume this is the time at each time sample
%     NEED to create the time vector... 
    %=== Calculate fit R2 with 1 or 2 gaussians
    mean_v1D = diff(mean_path1D)./diff(mean_time1D);
    [fg1,gof1] = fit(mean_time1D(2:end)',mean_v1D','gauss1','StartPoint',[6 1 1]);
    [fg2,gof2] = fit(mean_time1D(2:end)',mean_v1D','gauss2','StartPoint',[6 1 1 6 1 1]);
    gaussR2(ii) = gof1.rsquare;
    
end
figure; scatter(med_curv,gaussR2)

%=== Cluster flight types
flight_type = kmeans([med_curv,gaussR2],2);
figure('units','normalized','outerposition',[.7 0.2 .13 0.25]);
gscatter(med_curv,gaussR2,flight_type); set(gca,'XScale','log');    axis square;
xlabel('Median Curvature'); ylabel('R square (fit)');

%=== assign flight type to each replay 
replay_flight_type_num = NaN(1,numel(RP.post))
cluster_idx_end = [cluster_idx numel(RP.post)]
for i = 2:numel(flight_type)
    replay_flight_type_num(cluster_idx(i-1):cluster_idx(i)-1) = flight_type(i-1)
end
end
%% sigmoid fit on the replays [slow]
for hide = 1
logisEqn = '1/(1+exp(-b*(x-c)))';
cut_size = 15;
new_post={};
for i = 1:numel(RP.post)
% for i = 1457
    posterior = RP_post_trimmed{i}';                % Get the posterior probability (cut the edges)
    if isempty(posterior)
        continue;
    end
    nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
    pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial'); % Smooth the posterior
    pstP = [nan(nS, cut_size), pstP, nan(nS, cut_size)];            % pad the posterior with nan
    [max_P,max_L] = max(pstP);                                      % Calculate location at max posterior
    invalidP = (max_P==0 | isnan(max_P));                           % Invalidate NaNs or zeros
    y1 = normalize(max_L,'range')';                                 % Normalize between 0 and 1 and assign x,y and weights
    x1 = [1:numel(y1)]';  w = max_P;                                 % Assign x and weights
    y1(1) = 0;   y1(end)=1;                                          % Anchor to start and stop (to help fitting)
    x1(invalidP)=[]; w(invalidP)=[]; y1(invalidP)=[];                % Remove invalid points
    [f,gof] = fit(x1,y1,logisEqn,'Start',[3 x1(end)/2],'Weights',w); % Fit with logistic function
    sigmoid_center = floor(f.c);
    if sigmoid_center-cut_size <=0 || sigmoid_center+cut_size > size(pstP,2)
        new_post{i} = nan(nS, cut_size*2+1);
        continue
    end
    new_post{i} = [pstP(:,sigmoid_center-cut_size:sigmoid_center+cut_size)];
end
end
%% calculating the movement of replay 
post = RP_post_trimmed;
movement_all_replay = zeros(numel(RP_post_trimmed),size(new_post{i},2)-1);
for i = 1:numel(RP_post_trimmed)
    event = new_post{i}';
    if isempty(event)
        continue
    end
    event = imgaussfilt(event,[.01 1],'FilterDomain','spatial');        % Smooth the posterior
    nS = size(event,2);
    nT = size(event,1);
    uniform_prob = 1/nS; 
    event_position_vector = [1:size(event,2)];
    event_pos_vector_rescaled = rescale (event_position_vector);        % rescaling y
    event(find(isnan(event)))=1/nS;                                     % assign NaN to uniform prob
    [position_prob, decoded_position] = max(event');                    % finding decoded position with max prob
    decoded_position = event_pos_vector_rescaled(decoded_position);      % rescaling decoded pos
    decoded_position(find(position_prob < uniform_prob*3)) = NaN;       % filter out position with low probability 
    event_time_vector = [1:numel(decoded_position)];                    % create time vector 
    event_time_mid = event_time_vector(1:end-1)+0.5;                    % create time vector for difference
%     decoded_pos_smooth =  fillmissing(decoded_position,"linear");       % interpolate the NaN positions
%     movement = diff(decoded_pos_smooth);                                % find movement 
    movement = diff (decoded_position);                                 % do not smooth event 
%     speed = abs(movement);                                              % (%trajectory covered/bin)
    movement_all_replay (i,:) = movement;
    % visualization 
%     figure; subplot(2,1,1);
%     % imagesc(event');
%     imagesc(event_time_vector, event_pos_vector_rescaled,event');
%     axis xy;
%     hold on
%     scatter(event_time_vector, decoded_position)
%     xlim([0.5,nT+0.5])
%     subplot(2,1,2);plot(event_time_mid, speed);
%     xlim([0.5,nT+0.5])
end
%% plotting the speed-time profile 
% results vary greatly based on different values of replay_scores_replay,
% as well as sign of weighted_corr
for hide = []
speed_1 = zeros(1,nT-1);
speed_2 = zeros(1,nT-1);
figure('units','normalized','outerposition',[.2 0.2 .7 0.7]);
tiledlayout(3,4);
for ii =  [0:0.1:0.8]
    nexttile;
    speed_1 = zeros(1,nT-1);
    speed_2 = zeros(1,nT-1);
    condition =  weighted_corr_replay >0.5  & segment_len_frac >0.5 & replay_scores_replay >ii & replay_horiz == 0  & RP_clus_id ~=1 & shuffle_sig;
%     condition = weighted_corr_replay>0.5;
%     for i = find(replay_flight_type_num==1 & condition)
%         current = movement_all_replay{i};
%         current(isnan(current))=0;
%         speed_1 = speed_1 + current;
%     end
%     for i = find(replay_flight_type_num==2  & condition)
%         current = movement_all_replay{i};
%         current(isnan(current))=0;
%         speed_2 = speed_2 + current;
%     end
    %   loop through each position 
    for i = 1:size(movement_all_replay,2)
        test = movement_all_replay(:,i);
        speed_1(i) = mean(test(~isnan(test') & condition & replay_flight_type_num ==1));
    end
    for i = 1:size(movement_all_replay,2)
        test = movement_all_replay(:,i);
        speed_2(i) = mean(test(~isnan(test') & condition & replay_flight_type_num ==2));
    end
    % subplot(1,2,1);plot(rescale(speed_1));xlabel("Type 1")
    % subplot(1,2,2);plot(rescale(speed_2));xlabel("Type 2")
    plot(normalize(speed_1));hold on;
    plot(normalize(speed_2))
    title(["replay score > ", ii])
    ylabel("movement")
    xlabel("time bins")
    legend("Type 1","Type 2","Location", "northwest")
    end
end
%% cropping posteriors to around the center of the sigmoid 
for hide = []
center_space = 10; % 5 bins on either side of the center of sigmoid 
center_idx = ceil(size(new_post{1},2)/2);
center_post = {};
center_slope = [];
center_R2 = [];
for i = 1:numel(RP.post)
    event = new_post{i};
    if isempty(event)
        continue
    end
    center_post{i} = event(:,center_idx-center_space:center_idx+center_space); % saving only the center of events 
    % save both slope and goodness of fit  
    event = center_post{i};
    [position_prob, decoded_position] = max(event');
    event_time_vector = [1:numel(decoded_position)];
    event_position_vector = [1:size(event,2)];
    event_pos_vector_rescaled = rescale (event_position_vector);
    event_pos_scaled = event_pos_vector_rescaled(decoded_position);
%     line_slope_intercept = polyfit(event_time_vector, event_pos_scaled,1); 
%     center_slope (i) =  line_slope_intercept(1);
%     fitted_y = polyval(line_slope_intercept,event_time_vector);
    mdl = fitlm(event_time_vector, event_pos_scaled);
    center_R2 (i) = mdl.Rsquared.Ordinary;
    center_slope (i) =  mdl.Coefficients.Estimate(2);
end
end
%% plotting the heatmap, slope, R^2 of replay (whole or centered on sigmoid), by flight type
for hide = []
test=new_post; % plot the entire replay
% test=center_post; % plot only the center of replays 

new_post_sum = zeros(20,size(test{1},2));
% condition = weighted_corr_replay >0.8& segment_len_frac >0.6 & replay_scores_replay >0.4 & replay_horiz == 0  & RP_clus_id ~=1 & shuffle_sig;
 condition = weighted_corr_replay >0.2 & shuffle_sig& RP_dur >50  & segment_len_frac >0.7  & RP_clus_id ~=1 & spatial_coverage > 0.5;
for i = find(replay_flight_type_num==1)
    if condition(i)
        current = test{i};
        current(isnan(current))=0;
        current(current<0) = 0;
        current = imresize(current,[20, size(test{1},2)]);
        current=current./sum(current);
        current(isnan(current))=0;
        new_post_sum = new_post_sum + current;
    end
end
new_post_sum_2 = zeros(20,size(test{1},2));
for i = find(replay_flight_type_num==2)
    if condition(i)
        current = test{i};
        current(isnan(current))=0;
        current(current<0) = 0;
        current = imresize(current,[20, size(test{1},2)]);
        current=current./sum(current);
        current(isnan(current))=0;
        new_post_sum_2 = new_post_sum_2 + current;
    end
end
figure('units','normalized','outerposition',[.7 0.2 .3 0.3]);
subplot(1,2,1);imagesc(new_post_sum);title("Type 1"); xlabel("time bins"); ylabel("normalized position")
set(gca,'YTick', [])
axis xy;
subplot(1,2,2);imagesc(new_post_sum_2);title("Type 2");
set(gca, 'YTick', [])
axis xy;

% visualizing slope and R2 of the center region of replays 
% setting the condition 
% condition =  weighted_corr_replay >0  & segment_len_frac >0.5 & replay_scores_replay >0.4 & replay_horiz == 0  & RP_clus_id ~=1 & shuffle_sig;

% R2
figure('units','normalized','outerposition',[.7 0.4 .15 0.3]);
center_R2_flight1 = center_R2(find( ~isnan(center_R2) & replay_flight_type_num==1 & condition));
center_R2_flight2 = center_R2(find( ~isnan(center_R2) & replay_flight_type_num==2 & condition));
temp_length =  max(numel(center_R2_flight1),numel(center_R2_flight2));
center_R2_flight1 = [center_R2_flight1 nan(1,temp_length - numel(center_R2_flight1))];
center_R2_flight2 = [center_R2_flight2 nan(1,temp_length - numel(center_R2_flight2))];

PlotDistr_AF_v1([center_R2_flight1; center_R2_flight2],turbo(size(flighttype_comp,2)),"R^2",["Type 1" "Type2"]);
[~,p] = ttest(center_R2_flight1,center_R2_flight2)
signrank(center_R2_flight1,center_R2_flight2)
ranksum(center_R2_flight1,center_R2_flight2)
% slope
figure('units','normalized','outerposition',[.56 0.4 .15 0.3]);
center_slope_flight1 = center_slope(find( ~isnan(center_slope) & replay_flight_type_num==1 & condition));
center_slope_flight2 = center_slope(find( ~isnan(center_slope) & replay_flight_type_num==2 & condition));
temp_length =  max(numel(center_slope_flight1),numel(center_slope_flight2));
center_slope_flight1 = [center_slope_flight1 nan(1,temp_length - numel(center_slope_flight1))];
center_slope_flight2 = [center_slope_flight2 nan(1,temp_length - numel(center_slope_flight2))];

PlotDistr_AF_v1([center_slope_flight1; center_slope_flight2],turbo(size(flighttype_comp,2)),"slope",["Type 1" "Type2"]);
[~,p] = ttest2(center_slope_flight1,center_slope_flight2)
signrank(center_slope_flight1,center_slope_flight2)

ranksum(center_slope_flight1,center_slope_flight2)
end
%% plotting the maximum posterior line of diff flight types on the same graph 
for hide = []
%=== calculating the trajectory of each summed trajectory 
nS = size(new_post_sum,1);
nT = size(new_post_sum,2);

new_post_sum_1_norm = new_post_sum./sum(new_post_sum)
[position_prob_1, decoded_position_1] = max(new_post_sum_1_norm);
event_time_vector_1 = [1:numel(decoded_position_1)];


new_post_sum_2_norm = new_post_sum_2./sum(new_post_sum_2)
[position_prob_2, decoded_position_2] = max(new_post_sum_2_norm);
event_time_vector_2 = [1:numel(decoded_position_2)];

%=== calculating the width of the line 

% condition is 2 times the uniform probability 
for hide = []
broad_line_matrix_1 = new_post_sum_1_norm>uni_prob*2;
strfind(broad_line_matrix_1(:,3)', [0 1])
% lower_bound_1 = nan(1,nT);
% upper_bound_1 = nan(1,nT);
lower_bound_1 = decoded_position_1;
upper_bound_1 = decoded_position_1;
for i = 1: nT
    current_vector = broad_line_matrix_1(:,i)';
    current_boundary =  strfind(current_vector, [0 1]);
    if ~isempty(current_boundary)
        lower_bound_1(i) = current_boundary;
    end
    current_boundary =  strfind(current_vector, [1 0]);
    if ~isempty(current_boundary)
        upper_bound_1(i) = current_boundary;
    end
end
broad_line_matrix_2 = new_post_sum_2_norm>uni_prob*2;
strfind(broad_line_matrix_2(:,3)', [0 1])
lower_bound_2 = decoded_position_2;
upper_bound_2 = decoded_position_2;
for i = 1: nT
    current_vector = broad_line_matrix_2(:,i)';
    current_boundary =  strfind(current_vector, [0 1]);
    if ~isempty(current_boundary)
        lower_bound_2(i) = current_boundary(1);
    end
    current_boundary =  strfind(current_vector, [1 0]);
    if ~isempty(current_boundary)
        upper_bound_2(i) = current_boundary(1);
    end
end
end 
% condition is percent of maximum
percent_max = 0.9;
threshold = max(new_post_sum_1_norm)*percent_max;
lower_bound_1 = decoded_position_1;
upper_bound_1 = decoded_position_1;
for i = 1: nT
   current_vector =new_post_sum_1_norm(:,i)';
   compare = current_vector > threshold(i);
   if isempty(strfind(compare,[0 1]))
       lower_bound_1(i)=0.5;
   else
       lower_bound_1(i) = strfind(compare,[0 1])+0.5;
   end
   if isempty(strfind(compare,[1 0]))
       upper_bound_1(i) = nS-0.5;
   else
       upper_bound_1(i) = strfind(compare,[1 0])+0.5;
   end
end

threshold = max(new_post_sum_2_norm)*percent_max;
lower_bound_2 = decoded_position_2;
upper_bound_2 = decoded_position_2;
for i = 1: nT
   current_vector =new_post_sum_2_norm(:,i)';
   compare = current_vector > threshold(i);
   if isempty(strfind(compare,[0 1]))
       lower_bound_2(i)=0.5;
   else
       lower_bound_2(i) = strfind(compare,[0 1])+0.5;
   end
   if isempty(strfind(compare,[1 0]))
       upper_bound_2(i) = nS-0.5;
   else
       upper_bound_2(i) = strfind(compare,[1 0])+0.5;
   end
end



figure('units','normalized','outerposition',[.7 0.2 .3 0.3]);
subplot(1,3,1);imagesc(new_post_sum);title("Type 1"); xlabel("time bins"); ylabel("normalized position")
hold on; 
plot(event_time_vector_1,decoded_position_1)
set(gca,'YTick', [])
axis xy;
subplot(1,3,2);imagesc(new_post_sum_2);title("Type 2");
hold on; 
plot(event_time_vector_2,decoded_position_2)
set(gca, 'YTick', [])
axis xy;
subplot(1,3,3);
plot(event_time_vector_1,decoded_position_1); hold on;
plot(event_time_vector_2,decoded_position_2); hold on;
fill_y = [lower_bound_1 flip(upper_bound_1)];
% nonnan_idx =  ~isnan( fill_y);
fill_x = [event_time_vector_1 flip(event_time_vector_1)];
fill ( fill_x, fill_y,[0 0 0],'FaceAlpha',0.3);hold on;
fill ( [event_time_vector_2 flip(event_time_vector_2)], [lower_bound_2 flip(upper_bound_2)],'r','FaceAlpha',0.3)

ylim([0,nS+0.5])
xlim([0,nT+0.5])
axis xy;

legend("Type 1","Type 2","Type 1","Type 2");
end

%% plotting the heatmap, slope, R^2 of replay (whole or centered on sigmoid), by FF/NF/FN
for hide = []
test=new_post; % plot the entire replay
% test=center_post; % plot only the center of replays 

condition = weighted_corr_replay <-0.4& segment_len_frac >0.8 & replay_scores_replay >0.4  & RP_clus_id ~=1 & shuffle_sig;
flight_FF = strcmp(RP.flight_type, "feeder-feeder");
flight_FN = strcmp(RP.flight_type, "feeder-nonfeeder");
flight_NF = strcmp(RP.flight_type, "nonfeeder-feeder");

nS = 30;
nT = size(test{1},2);
% condition = weighted_corr_replay >0   & segment_len_frac >0.5
new_post_sum_FF = zeros(nT,nS);
new_post_sum_FN = zeros(nT,nS);
new_post_sum_NF = zeros(nT,nS);

for i = find(flight_FF)
    if condition(i)
        current = test{i};
        current = imresize(current,[nT, nS]);
        current=current./sum(current);
        
        current(isnan(current))=0;
        new_post_sum_FF = new_post_sum_FF + current;
    end
end

for i = find(flight_FN)
    if condition(i)
        current = test{i};
        current(isnan(current))=0;
        current = imresize(current,[nT, nS]);
        current=current./sum(current);
        
        current(isnan(current))=0;
        new_post_sum_FN = new_post_sum_FN + current;
    end
end

for i = find(flight_NF)
    if condition(i)
        current = test{i};
        current(isnan(current))=0;
        current = imresize(current,[nT, nS]);
        current=current./sum(current);
        
        current(isnan(current))=0;
        new_post_sum_NF = new_post_sum_NF + current;
    end
end


new_post_sum_FF_norm = new_post_sum_FF./sum(new_post_sum_FF,2)
new_post_sum_FN_norm = new_post_sum_FN./sum(new_post_sum_FN,2)
new_post_sum_NF_norm = new_post_sum_NF./sum(new_post_sum_NF,2)

[position_prob_FF, decoded_position_FF] = max(new_post_sum_FF');
[position_prob_FN, decoded_position_FN] = max(new_post_sum_FN');
[position_prob_NF, decoded_position_NF] = max(new_post_sum_NF');

figure('units','normalized','outerposition',[.7 0.2 .3 0.5]);
subplot(2,3,1);imagesc(new_post_sum_FF');title("Feeder-feeder"); xlabel("time bins"); ylabel("normalized position")
set(gca,'YTick', [])
axis xy;
subplot(2,3,2);imagesc(new_post_sum_FN');title("Feeder-nonfeeder");
set(gca, 'YTick', [])
axis xy;
subplot(2,3,3);imagesc(new_post_sum_NF');title("Nonfeeder-feeder");
set(gca, 'YTick', [])
axis xy;
% subplot(2,3,4); plot([1:nS],sum(new_post_sum_FF_norm))
% subplot(2,3,5); plot([1:nS],sum(new_post_sum_FN_norm))
% subplot(2,3,6); plot([1:nS],sum(new_post_sum_NF_norm))
subplot(2,3,4); plot([1:nT],position_prob_FF); xlabel ("time bins"); ylabel("maximum probability (AU)"); set(gca,'YTick', [])
subplot(2,3,5); plot([1:nT],position_prob_FN);set(gca,'YTick', [])
subplot(2,3,6); plot([1:nT],position_prob_NF);set(gca,'YTick', [])


end
%% plotting the heatmap, slope, R^2 of replay (whole or centered on sigmoid), by to feeder/not to feeder
for hide = []
test=new_post; % plot the entire replay
% test=center_post; % plot only the center of replays 

condition = weighted_corr_replay >0.3& segment_len_frac >0.7 & replay_scores_replay >0.2  & RP_clus_id ~=1 & shuffle_sig;
flight_to_feeder = strcmp(RP.flight_type,'nonfeeder-feeder')+strcmp(RP.flight_type,'feeder-feeder');
flight_to_nonfeeder = strcmp(RP.flight_type,'feeder-nonfeeder');

nS = 30;
nT = size(test{1},2);
% condition = weighted_corr_replay >0   & segment_len_frac >0.5
new_post_sum_to_feeder = zeros(nT,nS);
new_post_sum_to_nonfeeder = zeros(nT,nS);
new_post_sum_NF = zeros(nT,nS);

for i = find(flight_to_feeder)
    if condition(i)
        current = test{i};
        current = imresize(current,[nT, nS]);
        current=current./sum(current);
        
        current(isnan(current))=0;
        new_post_sum_to_feeder = new_post_sum_to_feeder + current;
    end
end

for i = find(flight_to_nonfeeder)
    if condition(i)
        current = test{i};
        current(isnan(current))=0;
        current = imresize(current,[nT, nS]);
        current=current./sum(current);
        
        current(isnan(current))=0;
        new_post_sum_to_nonfeeder = new_post_sum_to_nonfeeder + current;
    end
end


[position_prob_FF, decoded_position_FF] = max(new_post_sum_to_feeder');
[position_prob_FN, decoded_position_FN] = max(new_post_sum_to_nonfeeder');

figure('units','normalized','outerposition',[.7 0.2 .3 0.5]);
subplot(2,3,1);imagesc(new_post_sum_to_feeder');title("To feeder"); xlabel("time bins"); ylabel("normalized position")
set(gca,'YTick', [])
axis xy;
subplot(2,3,2);imagesc(new_post_sum_to_nonfeeder');title("To nonfeeder");
set(gca, 'YTick', [])
axis xy;
% subplot(2,3,4); plot([1:nS],sum(new_post_sum_FF_norm))
% subplot(2,3,5); plot([1:nS],sum(new_post_sum_FN_norm))
% subplot(2,3,6); plot([1:nS],sum(new_post_sum_NF_norm))
subplot(2,3,4); plot([1:nT],position_prob_FF); xlabel ("time bins"); ylabel("maximum probability (AU)"); set(gca,'YTick', [])
subplot(2,3,5); plot([1:nT],position_prob_FN);set(gca,'YTick', [])


end
end
%% histogram of replay lengths by flight type
for hide = []
flight_FF = strcmp(RP.flight_type, "feeder-feeder");
flight_FN = strcmp(RP.flight_type, "feeder-nonfeeder");
flight_NF = strcmp(RP.flight_type, "nonfeeder-feeder");
NF_replay = find(shuffle_sig & segment_len_frac > 0.5 & top_replay & flight_NF)
FN_replay = find(shuffle_sig & segment_len_frac > 0.5 & top_replay & flight_FN)
FF_replay = find(shuffle_sig & segment_len_frac > 0.5 & top_replay & flight_FF)

NF_flight_len = RP.flight_length(NF_replay);
FN_flight_len = RP.flight_length(FN_replay);
FF_flight_len = RP.flight_length(FF_replay);
bins = [0:0.5:15]
figure;
histogram(cell2mat(NF_flight_len),bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');hold on;
histogram(cell2mat(FN_flight_len),bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');hold on;
histogram(cell2mat(FF_flight_len),bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');
xlabel("flight length");
ylabel("normalized count of replays");
legend("Nonfeeder-feeder","Feeder-nonfeeder","feeder-feeder")
end
%% histograms of distance to prev/next flight (any flight cluster)
for hide = []
%===(normalized)
xlim_size = max([max(abs(RP.replay_flight_prev_dist_any)),max(abs(RP.replay_flight_next_dist_any))]);
xlim_size = 30;
hist_bins_num = 30;
hist_bins_gap = xlim_size*2/hist_bins_num;
hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
% hist_bins =  [-xlim_size: 60: xlim_size];
% condition = abs(weighted_corr_replay)>0  & replay_scores_replay >0.4 & replay_horiz == 0 & RP_clus_id ~=1 & segment_len_frac>0.5 ;
condition = top_replay;

figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);
% xlim_size = max(abs(RP.replay_flight_prev_dist));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist_any(find(RP.corr>0 & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0  & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.2])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")
% xlim_size = max(abs(RP.replay_flight_next_dist));
subplot(1,2,1); histogram (RP.replay_flight_next_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.18])

legend("forward replays","reverse replays");
xlabel("time to next flight (s)")


sgtitle("any cluster flights");
%=== (not normalized)

for hide = []
% xlim_size = max([max(abs(RP.replay_flight_st_dist)), max(abs(RP.replay_flight_en_dist)),max(abs(RP.replay_flight_prev_dist_any)),max(abs(RP.replay_flight_next_dist_any))]);
% xlim_size = 10;
% hist_bins_num = 25;
% hist_bins_gap = xlim_size*2/hist_bins_num;
% hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);

% xlim_size = max(abs(RP.replay_flight_prev_dist_any));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.15])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")

% xlim_size = max(abs(RP.replay_flight_next_dist_any));
subplot(1,2,1); histogram (RP.replay_flight_next_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.15])
legend("forward replays","reverse replays");
xlabel("time to next flight (s)")

end
end
%%(replicating caitlin's results) histograms of distance to prev/next flight BY flight type (same cluster) 
for hide = []
%===(normalized)
xlim_size = max([max(abs(RP.replay_flight_prev_dist)),max(abs(RP.replay_flight_next_dist))]);
xlim_size = 30;
hist_bins_num = 30;
hist_bins_gap = xlim_size*2/hist_bins_num;
hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
% hist_bins =  [-xlim_size: 60: xlim_size];

figure('units','normalized','outerposition',[0 0.05 0.4 0.91]);
subplot(4,2,1); histogram (RP.replay_flight_next_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');title("all flights")
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
xlabel("time to next flight (s)")
subplot(4,2,2); histogram (RP.replay_flight_prev_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")

subplot(4,2,3); histogram (RP.replay_flight_next_dist(find(flight_FF & RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
title("feeder-feeder")
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
subplot(4,2,4); histogram (RP.replay_flight_prev_dist(find(flight_FF & RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])

subplot(4,2,5); histogram (RP.replay_flight_next_dist(find(flight_NF & RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
title("nonfeeder-feeder")
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
subplot(4,2,6); histogram (RP.replay_flight_prev_dist(find(flight_NF & RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])


subplot(4,2,7); histogram (RP.replay_flight_next_dist(find(flight_FN & RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
title("feeder-nonfeeder")
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
subplot(4,2,8); histogram (RP.replay_flight_prev_dist(find(flight_FN & RP.corr>0.8 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
sgtitle("same cluster flights");






%=== (not normalized)

for hide = []
% xlim_size = max([max(abs(RP.replay_flight_st_dist)), max(abs(RP.replay_flight_en_dist)),max(abs(RP.replay_flight_prev_dist)),max(abs(RP.replay_flight_next_dist))]);
% xlim_size = 10;
% hist_bins_num = 25;
% hist_bins_gap = xlim_size*2/hist_bins_num;
% hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);

% xlim_size = max(abs(RP.replay_flight_prev_dist));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")

% xlim_size = max(abs(RP.replay_flight_next_dist));
subplot(1,2,1); histogram (RP.replay_flight_next_dist(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
legend("forward replays","reverse replays");
xlabel("time to next flight (s)")
end
end
%% histograms of distance to prev/next flight (any flight cluster), by flight type
for hide = []
%===(normalized)
xlim_size = max([max(abs(RP.replay_flight_prev_dist_any)),max(abs(RP.replay_flight_next_dist_any))]);
xlim_size = 30;
hist_bins_num = 30;
hist_bins_gap = xlim_size*2/hist_bins_num;
hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
% hist_bins =  [-xlim_size: 60: xlim_size];
condition = top_replay ;


figure('units','normalized','outerposition',[0.25 0.05 0.4 0.91]);
subplot(4,2,1); histogram (RP.replay_flight_next_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.13])
title("all flights");
subplot(4,2,2); histogram (RP.replay_flight_prev_dist_any(find(RP.corr>0 & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0  & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.13])
legend("forward replays","reverse replays");

subplot(4,2,3); histogram (RP.replay_flight_next_dist_any(find(flight_FF & RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.13])
title("Feeder-feeder")
subplot(4,2,4); histogram (RP.replay_flight_prev_dist_any(find(flight_FF & RP.corr>0 & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0  & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.13])

subplot(4,2,5); histogram (RP.replay_flight_next_dist_any(find(flight_NF &RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.13])
title("Nonfeeder-feeder")
subplot(4,2,6); histogram (RP.replay_flight_prev_dist_any(find(flight_NF &RP.corr>0 & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0  & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.13])

subplot(4,2,7); histogram (RP.replay_flight_next_dist_any(find(flight_FN & RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.13])
xlabel("time to next flight (s)")
title("Feeder-nonfeeder")
subplot(4,2,8); histogram (RP.replay_flight_prev_dist_any(find(flight_FN & RP.corr>0 & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0  & condition & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.13])
xlabel("time to previous flight (s)")

sgtitle("any cluster flights");


%=== (not normalized)

for hide = []
% xlim_size = max([max(abs(RP.replay_flight_st_dist)), max(abs(RP.replay_flight_en_dist)),max(abs(RP.replay_flight_prev_dist_any)),max(abs(RP.replay_flight_next_dist_any))]);
% xlim_size = 10;
% hist_bins_num = 25;
% hist_bins_gap = xlim_size*2/hist_bins_num;
% hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);

% xlim_size = max(abs(RP.replay_flight_prev_dist_any));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_prev_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim([0 0.15])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")

% xlim_size = max(abs(RP.replay_flight_next_dist_any));
subplot(1,2,1); histogram (RP.replay_flight_next_dist_any(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist_any(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist_any) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim([0 0.15])
legend("forward replays","reverse replays");
xlabel("time to next flight (s)")

end
end
%% histograms of distance to prev/next flight (same cluster),by flight type
for hide = []
%===(normalized)
xlim_size = max([max(abs(RP.replay_flight_prev_dist)),max(abs(RP.replay_flight_next_dist))]);
xlim_size = 30;
hist_bins_num = 20;
hist_bins_gap = xlim_size*2/hist_bins_num;
hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
% hist_bins =  [-xlim_size: 60: xlim_size];

figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);
% xlim_size = max(abs(RP.replay_flight_prev_dist));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim ([0,0.25])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")
% xlim_size = max(abs(RP.replay_flight_next_dist));
subplot(1,2,1); histogram (RP.replay_flight_next_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim ([0,0.25])

legend("forward replays","reverse replays");
xlabel("time to next flight (s)")


sgtitle("same cluster flights");
%=== (not normalized)

for hide = 1
% xlim_size = max([max(abs(RP.replay_flight_st_dist)), max(abs(RP.replay_flight_en_dist)),max(abs(RP.replay_flight_prev_dist)),max(abs(RP.replay_flight_next_dist))]);
% xlim_size = 10;
% hist_bins_num = 25;
% hist_bins_gap = xlim_size*2/hist_bins_num;
% hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);

% xlim_size = max(abs(RP.replay_flight_prev_dist));
subplot(1,2,2); histogram (RP.replay_flight_prev_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_prev_dist(find(RP.corr<0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")

% xlim_size = max(abs(RP.replay_flight_next_dist));
subplot(1,2,1); histogram (RP.replay_flight_next_dist(find(RP.corr>0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; histogram (RP.replay_flight_next_dist(find(RP.corr<-0.8 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins, 'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
legend("forward replays","reverse replays");
xlabel("time to next flight (s)")
end
end
%% histograms of distance to prev/next flight (same cluster, normalize by cluster and then average)
for hide = []
%===(normalized)
xlim_size = max([max(abs(RP.replay_flight_prev_dist)),max(abs(RP.replay_flight_next_dist))]);
xlim_size = 30;
hist_bins_num = 20;
hist_bins_gap = xlim_size*2/hist_bins_num;
hist_bins = [-xlim_size: hist_bins_gap: xlim_size];
hist_bins_c = [-xlim_size+hist_bins_gap/2: hist_bins_gap: xlim_size-hist_bins_gap/2];

replays_current = RP.corr>0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size;
cluster_ids_current = flight_ids(find(replays_current_fwd));
cluster_ids_curr_unique = unique(cluster_ids_current);
replay_count_fwd_prev = zeros(1,numel(hist_bins)-1);
for i = cluster_ids_curr_unique
    current_replays = RP.replay_flight_prev_dist(replays_current & flight_ids == i);
    replay_count_fwd_prev  = replay_count_fwd_prev + histcounts(current_replays,hist_bins,'Normalization','probability');
end
replay_count_fwd_prev = replay_count_fwd_prev ./ numel(cluster_ids_curr_unique);

replays_current = RP.corr<0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size;
cluster_ids_current = flight_ids(find(replays_current));
cluster_ids_curr_unique = unique(cluster_ids_current);
replay_count_rev_prev = zeros(1,numel(hist_bins)-1);
for i = cluster_ids_curr_unique
    current_replays = RP.replay_flight_prev_dist(replays_current & flight_ids == i);
    replay_count_rev_prev  = replay_count_rev_prev + histcounts(current_replays,hist_bins,'Normalization','probability');
end
replay_count_rev_prev = replay_count_rev_prev ./ numel(cluster_ids_curr_unique);


replays_current = RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size;
cluster_ids_current = flight_ids(find(replays_current));
cluster_ids_curr_unique = unique(cluster_ids_current);
replay_count_fwd_next = zeros(1,numel(hist_bins)-1);
for i = cluster_ids_curr_unique
    current_replays = RP.replay_flight_next_dist(replays_current & flight_ids == i);
    replay_count_fwd_next  = replay_count_fwd_next + histcounts(current_replays,hist_bins,'Normalization','probability');
end
replay_count_fwd_next = replay_count_fwd_next ./ numel(cluster_ids_curr_unique);

replays_current = RP.corr<0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size;
cluster_ids_current = flight_ids(find(replays_current));
cluster_ids_curr_unique = unique(cluster_ids_current);
replay_count_rev_next = zeros(1,numel(hist_bins)-1);
for i = cluster_ids_curr_unique
    current_replays = RP.replay_flight_next_dist(replays_current & flight_ids == i);
    replay_count_rev_next  = replay_count_rev_next + histcounts(current_replays,hist_bins,'Normalization','probability');
end
replay_count_rev_next = replay_count_rev_next ./ numel(cluster_ids_curr_unique);


figure('units','normalized','outerposition',[0.25 .25 0.5 0.5]);
% xlim_size = max(abs(RP.replay_flight_prev_dist));
subplot(1,2,2); bar (hist_bins_c,replay_count_fwd_prev,1,'k','edgecolor','none'); alpha(0.2);
% xlim ([-xlim_size xlim_size])
hold on; bar (hist_bins_c, replay_count_rev_prev, 1,'r','edgecolor','none'); alpha(0.2);
% histogram (RP.replay_flight_prev_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_prev_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([0 xlim_size])
ylim ([0,0.25])
legend("forward replays","reverse replays");
xlabel("time to previous flight (s)")
% xlim_size = max(abs(RP.replay_flight_next_dist));
subplot(1,2,1);
bar (hist_bins_c,replay_count_fwd_next,1,'k','edgecolor','none'); alpha(0.2);
% histogram (RP.replay_flight_next_dist(find(RP.corr>0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
% xlim ([-xlim_size xlim_size])
hold on; 
bar (hist_bins_c, replay_count_rev_next, 1,'r','edgecolor','none');alpha(0.2)
% histogram (RP.replay_flight_next_dist(find(RP.corr<-0 & top_replay & abs(RP.replay_flight_next_dist) <xlim_size)),hist_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
xlim ([-xlim_size 0])
ylim ([0,0.25])

legend("forward replays","reverse replays");
xlabel("time to next flight (s)")


sgtitle("same cluster flights");

end

%% plotting replays (other conditions)
for hide = []
    n_row=10;
    n_col=15;
    resize = 0;
    % defining best replay
    % best_replay = find(weighted_corr_replay<-0.7 & replay_scores_replay >0.5 & replay_horiz == 0 & RP_clus_id ~=1 );
%     best_replay = find(top_replay & RP_bat_id' == 14445);
%     best_replay = find(shuffle_sig & RP.replay_flight_prev_dist_any<2 & segment_len_frac > 0.5 & top_replay);
    condition = weighted_corr_replay >0 & RP_dur >50  & segment_len_frac >0.7  & RP_clus_id ~=1 & spatial_coverage > 0.5;
    best_replay = find(replay_flight_type_num==1 &  condition );
%     best_replay = find(shuffle_sig & segment_len_frac > 0.5 & top_replay & flight_FN & weighted_corr_replay >0);
    % randomly sampling a certain number of best replays
    if numel(best_replay) > n_row * n_col
        best_replay = datasample(best_replay, n_row*n_col, 'Replace',false);
    end
    % filter by graph order by flight length
    sort_by_variable = RP_flight_length;
    [~, sorted_index]=sort(abs(sort_by_variable), 'descend'); % ascend or descend
    top_idx = best_replay;
    top_idx=sorted_index(ismember(sorted_index,top_idx));
    
    % finding what to set as the maximum x axis value
    max_x = 0;
    for i = top_idx
        max_x = max(size(RP.post{i},1), max_x);
    end
    
    %visualization
    figure('units','normalized','outerposition',[.2 0 .6 1]);
    tiledlayout(n_row,n_col,'TileSpacing','compact');
    for ii = 1:n_row*n_col
        if ii>numel(top_idx)
            break;
        end
        
        i=top_idx(ii);
        nexttile;
        colormap(hot);
        %     current_post = imresize(RP.post{i},[size(RP.post{i},1), 10]); % make y time bins uniform
        current_post = RP_post_trimmed{i};
        if resize == 1
            current_post = imresize(current_post,[20, 20]);
            current_post = current_post./sum(current_post)
        end
        uniform_prob = 1/size(current_post,2);
        [position_prob, decoded_position] = max(current_post');
        current_post(find(position_prob' < 3*uniform_prob),:) = NaN;
        
        
        imagesc(current_post',prctile(current_post',[1 99],"all")'); hold on;
        xlabel([num2str(round(sort_by_variable(i),1)), ' m']);
        set(gca, 'YTick', []);
        set(gca, 'XTick', []);
        %     xlabel([num2str(size(RP.post{i},1)*5) ' ms']);
        axis xy;
        %     xlim([1, max_x-1]); % making all x axis consistent
        set(gca,'Color','black') % setting the extra volume as black
    end
end
%% plotting local/nonlocal replays 
for hide = 1
RP_replay_locations = string(num2str(nan(1,numel(RP.replay_type))));
RP_replay_locations(find(contains(RP.replay_type,"takeoff"))) = "takeoff";
RP_replay_locations(find(contains(RP.replay_type,"landing"))) = "landing";
RP_replay_locations(find(contains(RP.replay_type,"remote"))) = "remote";
RP_replay_locations(find(contains(RP.replay_type,"takeoff landing"))) = NaN;

RP_local_nonlocal =  string(num2str(nan(1,numel(RP.replay_type))));
RP_local_nonlocal(find(contains(RP_replay_locations,"takeoff") & weighted_corr_replay>0)) = "local"; % take off + forward
RP_local_nonlocal(find(contains(RP_replay_locations,"landing") & weighted_corr_replay<0)) = "local"; % landing + reverse 
RP_local_nonlocal(find(contains(RP_replay_locations,"takeoff") & weighted_corr_replay<0)) = "nonlocal"; % take off + forward
RP_local_nonlocal(find(contains(RP_replay_locations,"landing") & weighted_corr_replay>0)) = "nonlocal"; % landing + reverse 
RP_local_nonlocal(find(contains(RP_replay_locations,"remote"))) = "nonlocal"; % landing + reverse 



% loop through each cluster 
for hide = []
cluster_idx_all = [cluster_idx numel(RP.replay_type)];
local_frac_clusters = zeros(1, numel(cluster_idx));
nonlocal_frac_clusters = zeros(1, numel(cluster_idx));
for i = 1: numel(cluster_idx)
    if RP.cluster_id{cluster_idx(i)} == 1 || strcmp(RP.flight_type(cluster_idx(i)),'nonfeeder-nonfeeder')
        local_frac_clusters(i) = NaN;
        nonlocal_frac_clusters(i) = NaN;
        continue;
    end
    current_local = sum(strcmp (RP_local_nonlocal (cluster_idx_all(i): cluster_idx_all (i+1)),'local'));
    current_nonlocal = sum(strcmp (RP_local_nonlocal (cluster_idx_all(i): cluster_idx_all (i+1)),'nonlocal'));
    total_count = current_local + current_nonlocal;
    local_frac_clusters(i) = current_local / total_count; 
    nonlocal_frac_clusters(i) = current_nonlocal / total_count; 
end


figure;
PlotDistr_AF_v1([local_frac_clusters; nonlocal_frac_clusters],turbo(2),"Fraction replays",["Local" "Nonlocal"]);
[~,p] = ttest(local_frac_clusters,nonlocal_frac_clusters);
ps = signrank(local_frac_clusters,nonlocal_frac_clusters);
disp(['ttest p: ' num2str(p) '   signrank p: ' num2str(ps)])
end
% loop through each session
[~, session_idx] = unique(RP.session_id);
% loop through each session 
session_idx = session_idx';
session_idx_all = [session_idx numel(RP.replay_type)];
local_frac_sessions = zeros(1, numel(session_idx));
nonlocal_frac_sessions = zeros(1, numel(session_idx));
for i = 1: numel(session_idx)
    current_local = sum(strcmp (RP_local_nonlocal (session_idx_all(i): session_idx_all (i+1)),'local'));
    current_nonlocal = sum(strcmp (RP_local_nonlocal (session_idx_all(i): session_idx_all (i+1)),'nonlocal'));
    total_count = current_local + current_nonlocal;
    local_frac_sessions(i) = current_local / total_count; 
    nonlocal_frac_sessions(i) = current_nonlocal / total_count; 
end

% 
% PlotDistr_AF_v1([local_frac_sessions; nonlocal_frac_sessions],turbo(2),"Fraction replays",["Local" "Nonlocal"]);
% [~,p] = ttest(local_frac_sessions,nonlocal_frac_sessions);
% ps = signrank(local_frac_sessions,nonlocal_frac_sessions);
% disp(['ttest p: ' num2str(p) '   signrank p: ' num2str(ps)])

% actual analysis used
replays_to_use = good_replay;
% replays_to_use = top_replay;
local_frac_sessions_fw = zeros(1, numel(session_idx));
nonlocal_frac_sessions_fw = zeros(1, numel(session_idx));
for i = 1: numel(session_idx)
    condition = replays_to_use & weighted_corr_replay > 0;
    current_local = sum(condition(session_idx_all(i): session_idx_all (i+1)) & strcmp (RP_local_nonlocal(session_idx_all(i): session_idx_all (i+1)),'local'));
    current_nonlocal = sum(condition(session_idx_all(i): session_idx_all (i+1))  & strcmp (RP_local_nonlocal (session_idx_all(i): session_idx_all (i+1)),'nonlocal'));
    total_count = current_local + current_nonlocal;
    local_frac_sessions_fw(i) = current_local / total_count; 
    nonlocal_frac_sessions_fw(i) = current_nonlocal / total_count; 
end

local_frac_sessions_rev = zeros(1, numel(session_idx));
nonlocal_frac_sessions_rev = zeros(1, numel(session_idx));

ttest(local_frac_sessions_fw-0.5);
signrank(local_frac_sessions_fw-0.5);

for i = 1: numel(session_idx)
    condition = replays_to_use & weighted_corr_replay < 0;
    current_local = sum(condition(session_idx_all(i): session_idx_all (i+1)) & strcmp (RP_local_nonlocal(session_idx_all(i): session_idx_all (i+1)),'local'));
    current_nonlocal = sum(condition(session_idx_all(i): session_idx_all (i+1))  & strcmp (RP_local_nonlocal (session_idx_all(i): session_idx_all (i+1)),'nonlocal'));
    total_count = current_local + current_nonlocal;
    local_frac_sessions_rev(i) = current_local / total_count; 
    nonlocal_frac_sessions_rev(i) = current_nonlocal / total_count; 
end

ttest(local_frac_sessions_rev-0.5)

disp(['local forward fraction: ' num2str(mean(local_frac_sessions_fw))]);
disp(['nonlocal forward fraction: ' num2str(mean(nonlocal_frac_sessions_fw))]);
disp(['local rev fraction: ' num2str(mean(local_frac_sessions_rev))]);
disp(['nonlocal rev fraction: ' num2str(mean(nonlocal_frac_sessions_rev))]);
disp(['signrank rev: ', num2str(signrank(local_frac_sessions_rev-0.5))]);
disp(['signrank fw: ', num2str(signrank(local_frac_sessions_fw-0.5))]);


PlotDistr_AF_v1([local_frac_sessions_fw; local_frac_sessions_rev],turbo(2),"Fraction replays",["forward" "reverse"]);


% %%
% [~,~,RP.sessionID] =  unique(string(RP.session_id));                          % Id for each cluster from each session
% RP_local_nonlocal_binary = NaN(1,numel(RP.session_id));                   % 1 if local 
% RP_local_nonlocal_binary (find(strcmp(RP_local_nonlocal, 'local'))) =1;
% RP_local_nonlocal_binary (find(strcmp(RP_local_nonlocal, 'nonlocal'))) =0;
% RP.local_nonlocal = RP_local_nonlocal_binary';    
% % RP_tbl = table(RP);
% Rpl_sessionG = groupsummary(RP,'sessionID','mean',{'local_nonlocal'}); %Group Summary
% % Rpl_sessionG = groupsummary(RP,'sessionID','mean','local_nonlocal'); %Group Summary
end
%% short flights vs long flights violin plots 
top_replay = RP_clus_id~=1 & abs(weighted_corr_replay) > 0.4 & replay_scores_replay>0.4 & segment_len_frac>0.7 & shuffle_sig;
disp(['Number of top replays:   '  num2str(numel(find(top_replay)))])
disp(['Percent forward replay:   ' num2str(sum(weighted_corr_replay(top_replay)>0)/numel(find(top_replay))) ])
figure;
flight_len = cell2mat(RP.flight_length);
replay_short_idx = find(flight_len>3 & flight_len<7 );
replay_long_idx = find(flight_len>9 & flight_len<13);
flight_len_short = (flight_len(intersect(replay_short_idx, cluster_idx)));
flight_len_long = (flight_len(intersect(replay_long_idx, cluster_idx)));
top_replay_idx = find(top_replay);
% top_replay_idx = find(good_replay);
% loop through each cluster 
cluster_idx_all = [cluster_idx numel(RP.replay_type)];
flight_duration_med_long = [];
flight_duration_med_short =[];
wc_short = [];
wc_long = [];
score_short = [];
score_long = [];
for i = 1: numel(cluster_idx)
    if RP.cluster_id{cluster_idx(i)} == 1 || strcmp(RP.flight_type(cluster_idx(i)),'nonfeeder-nonfeeder')
        continue;
    end
    current_flight_len = flight_len(cluster_idx(i));
    if current_flight_len > 3 & current_flight_len <7
        replay_idx = (intersect(top_replay_idx,(cluster_idx_all(i): cluster_idx_all (i+1))));
        all_replay_dur = RP_dur (replay_idx);
        flight_duration_med_short = [flight_duration_med_short median(all_replay_dur)];
        all_replay_wc = weighted_corr_replay(replay_idx);
        wc_short = [wc_short median(abs(all_replay_wc))];
        all_replay_score = replay_scores_replay(replay_idx);
        score_short = [score_short median(all_replay_score)];
        continue;
    end
    if current_flight_len > 9 & current_flight_len <13
        replay_idx = (intersect(top_replay_idx,(cluster_idx_all(i): cluster_idx_all (i+1))));
        all_replay_dur = RP_dur (replay_idx);
        flight_duration_med_long = [flight_duration_med_long median(all_replay_dur)];
        all_replay_wc = weighted_corr_replay(replay_idx);
        wc_long = [wc_long median(abs(all_replay_wc))];
        all_replay_score = replay_scores_replay(replay_idx);
        score_long = [score_long median(all_replay_score)];
    end
end
plot_type = 'scatter'; % scatter or violin
subplot(1,4,1);plot_distr_AF_violin(flight_len_short, flight_len_long, {'Short Flights', 'Long Flights'}, 'SEM', 'Flight Length (m)', plot_type);      
currentYLim = ylim;  
% ylim([0, 14]);

subplot(1,4,2); plot_distr_AF_violin(flight_duration_med_short, flight_duration_med_long, {'Short Flights', 'Long Flights'}, 'SEM', 'Replay Duration (s)', plot_type);      
currentYLim = ylim;  
ylim([0, 400]);

subplot(1,4,3); plot_distr_AF_violin(abs(wc_short), abs(wc_long), {'Short Flights', 'Long Flights'}, 'SEM', 'Weighted Corr', plot_type);      
currentYLim = ylim;  
ylim([0, 1]);
subplot(1,4,4); plot_distr_AF_violin(abs(score_short), abs(score_long), {'Short Flights', 'Long Flights'}, 'SEM', 'Replay score', plot_type);      
currentYLim = ylim;  ylim([0, 1]);

%% quantifications 
disp(['Number of top replays:   '  num2str(numel(find(top_replay)))])
disp(['Percent forward replay:   ' num2str(sum(weighted_corr_replay(top_replay)>0)/numel(find(top_replay))) ])
disp(['Number of good replays:   '  num2str(numel(find(good_replay)))])
disp(['Percent forward replay:   ' num2str(sum(weighted_corr_replay(good_replay)>0)/numel(find(good_replay))) ])
disp(['good forward replay:   ' num2str(sum(weighted_corr_replay(good_replay)>0)) ])
disp(['good reverse replay:   ' num2str(sum(weighted_corr_replay(good_replay)<0)) ])

%% plotting RMS error 
% count = 0;
same_clus_error_all = [];
for i =  1:numel(sessions_info.same_clus_error)
    if isempty(sessions_info.same_clus_error{i})
        continue;
    end
    same_clus_error_all = [same_clus_error_all sessions_info.same_clus_error{i}{1}{1}]
%     count = count +1;
end
diff_clus_error_all = [];
for i =  1:numel(sessions_info.diff_clus_error)
    if isempty(sessions_info.diff_clus_error{i})
        continue;
    end
    diff_clus_error_all = [diff_clus_error_all sessions_info.diff_clus_error{i}{1}{1}]
end

figure;
histogram(same_clus_error_all,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    hold on;
histogram(diff_clus_error_all,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
xlabel("RMS error");
ylabel("fraction");
legend("same cluster decoding error flight","different cluster decoding error flight");
% count

save("same_clus_error_all", "same_clus_error_all");
save("diff_clus_error_all", "diff_clus_error_all");







