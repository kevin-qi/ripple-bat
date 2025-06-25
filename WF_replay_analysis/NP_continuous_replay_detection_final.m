%%=== LOAD DATA
for hide=1
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);
    nrnFile = dir('SU_kilosort*');
    NP_unit = table();
    MUA_unit = table();
    n_probes = size(nrnFile,1);    
    merge_cells = 1;
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        %=== If merging potential duplicates
        if merge_cells
%             load('merges.mat'); 
            mrgFile = dir('merges_*');              
            load([mrgFile.folder '/' mrgFile.name]);    % This contains the ids of potential merges
            if ~isempty(merges{i,1})
                for jj=1:size(merges{i,1},1)
                    toBeMerged = find(any(out.good_units.cluster_id == merges{i,1}{jj,1},2));
                    spikeTimes_usec = vertcat(out.good_units.spikeTimes_usec{toBeMerged,1});
                    temp_unit = out.good_units(end,:);
%                     temp_unit.cluster_id = temp_unit.cluster_id+1;
                    temp_unit.cluster_id = out.good_units.cluster_id(toBeMerged(1));
                    temp_unit.spikeTimes_usec = {spikeTimes_usec};
                    out.good_units = [out.good_units;temp_unit];
                end
                toBeDeleted = find(any(out.good_units.cluster_id == horzcat(merges{i,1}{:}),2));
                out.good_units(toBeDeleted,:) = [];
            end
        end
        NP_unit = [NP_unit; out.good_units];                                            % Single units are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);
    %=== Load LFP from Probe1
    load('LFP_probe1.mat');                                                             % This file contains LFP data from probe 1
    unique_ID = options.unique_ID;      
    
    %=== Create analysis folder for storing the results
    options.savefigures = 1;                                   % If creating folder for saving the data or not
    options.figure_show = 'off';
    current_version = 'v37';
    if options.savefigures == 1
        analysis_directory=fullfile(pwd,[current_version,'_merge_test_Decoding_analysis_figures_',datestr(now, 'yymmdd_HHMM')]);
        if ~exist(analysis_directory,'dir')
            mkdir(analysis_directory);
        end
    end
end 
%% loading RPL data 
for hide =1 
FileList = dir(fullfile(cd,"RPL*"));
all_RPL = [];
RPL_table = [];
RPL_n = [];
for nc = 1:length(FileList)
    all_RPL = [all_RPL load(FileList(nc).name)];
    RPL_table = [RPL_table;all_RPL(nc).RPL_out.table];
    RPL_n(nc) = size(all_RPL(nc).RPL_out.table,1);
end
% [~,RPL_max_idx] = max(RPL_n);
RPL_max_idx = 2;

% avoid duplicates of SWR
%=== Sort by time, find burst initiators and keep largest amplitude events
RPL = table();
candidate_RPL = sortrows(RPL_table,'t','ascend');
candidate_RPL.brst = [0; diff(candidate_RPL.t)<0.05];
initiators = find(~candidate_RPL.brst);
for ii = 1:numel(initiators)-1
    subtable = candidate_RPL(initiators(ii):initiators(ii+1)-1,:);
    [~,tmp_idx] = max(subtable.amp);
    RPL = [RPL; subtable(tmp_idx,:)];
end

%=== get all RPL times and duration  
RPL_t_all =  RPL.t;
LFP_rp = bandpass(all_RPL(RPL_max_idx).RPL_out.LFP_trace,[100 200],all_RPL(RPL_max_idx).RPL_out.Fs);                                           % Ripple trace



end
%% RECLUSTER FLIGHTS
for hide=1
    disp('Reclustering flights...');
    fig_count = 1;                                          % Id of the first figure
    alpha_clus = .8;   clear f_clus    % Default is 0.7
    f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);

end
%% BASIC PARAMS, ASSIGNMENTS AND SANITY CHECKS
for hide=1
    %=== Params
    if isrow(t),t=t';end                                    % Make sure t is a column vector
    NP_unit = NP_unit(NP_unit.fr<1,:);                      % Exclude neurons with high-firing rate
    %NP_unit = NP_unit(NP_unit.Amplitude>1e3,:);            % Exclude neurons with low amplitude
    n_cells = size(NP_unit,1);                              % Number of units
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 100;                                            % Number of repetitions for shuffling (Place Cells)
    times_to_shift = t(randi([0.5*Fs T-0.5*Fs],1,n_rep));   % Time intervals for circshift
%     fig_count = 1;                                          % Id of the first figure
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    rstr_int = [-3 3];                                      % Interval for looking at rasters
    t_tko = t(f_smp(:,1));                                  % Times of takeoff
    t_lnd = t(f_smp(:,2));                                  % Times of landing
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
    plot_field1D = 0;                                       % If plotting the 1D fields
    clus_id = 3;                                            % Identity of the flight cluster for place-cell analysis
    n2show = 5;                                             % Number of place cells to show
    Fs_imu = NP_imu.Fs;                                     % Sampling frequency of the IMU data
    %=== Populate the s cell with the single units and calculate firing rate
    s = cell(n_cells,n_rep);            % Single units' spikes
    Rate = zeros(length(t),n_cells);    % Smoothed Firing Rate
    for nc = 1:n_cells
        s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
        %=== Clean up duplicated spikes
        too_close = find(diff(s{nc,1})<0.001);
        if ~isempty(too_close)
            disp([num2str(numel(too_close)), ' cleaned spikes']);
        s{nc,1}(too_close+1) = [];
        end
        Rate(2:end,nc) = histcounts(s{nc,1},t)*Fs;                      % Calculate firing rate at behavioral samples
        Rate(:,nc) = smoothdata(Rate(:,nc),1,'gaussian',round(Fs*0.1)); % Smooth firing rate
        for n = 1:n_rep;    s{nc,n+1} = mod(s{nc,1}+times_to_shift(n),t(T)); end        % Shuffling
    end
    Rate_norm = normalize(Rate,1,'range');
    NP_unit = table2struct(NP_unit);                %=== Convert NP unit to structure
    %=== Populate the mua cell with the single units and calculate firing rate
    n_mua = size(MUA_unit,1);              % Number of mua units
    s_mua = cell(n_mua,1);                    % Single units' spikes
    for nc = 1:n_mua
        s_mua{nc,1} = MUA_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
        too_close = find(diff(s_mua{nc,1})<0.001);                          % Clean up duplicated spikes
        if ~isempty(too_close)
            %disp([num2str(numel(too_close)), ' cleaned spikes']);
            s_mua{nc,1}(too_close+1) = [];
        end
    end
    MUA_unit = table2struct(MUA_unit);                    % Convert MUA unit to structure
    %=== Transform the LFP into double and interpolate at 500 Hz sampling
    LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,NP_imu.t,'linear','extrap');     % Interpolate LFP at the accelerometer time samples
    LFP_mn = normalize(mean(LFP,2),'range');                                                                 % Calculate average LFP
    LFP_th = bandpass(LFP_mn,[4 12],Fs_imu);                                                                 % Filtered in the theta band
     %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_imu);                 % Filtered at the wing-beat frequency
    [up,lo] = envelope(a_flt_NP,round(0.06*Fs_imu),'peak');     % Upper and lower envelopes
    env = normalize(up - lo,'range');                           % Amplitude of the envelope
    env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
    wBeats = movsum(env>env_th,2*Fs_imu)>Fs_imu/5;              % Euristic criterion for flight detection
    %=== Calculate Population Spike Density
    all_s = sort([vertcat(s{:,1});vertcat(s_mua{:,1})]);
    t_rate = [t(1):0.001:t(end)];
    all_rate = kernel_rate_AF_v1(all_s,0.05,t_rate)/n_cells;
end
%% flight cluster calculations 
for hide = 1
imp_bat_clusters = unique(f_clus.id); % get all clusters 
n_surv_clusters = numel(unique(f_clus.id));
col_clus = hsv(n_surv_clusters);
% clusters_of_interest = setdiff(imp_bat_clusters,1); % exclude the flights that are unclustered
clusters_of_interest = imp_bat_clusters; % including analysis on cluster 1, as a control and bc im curious
r_lim = [-4.6 5.5; -2.3 1.8; 0 2.8]*1.1;           % Override Room boundaries

count_subplot = 1;

figure('units','normalized','outerposition',[.05 .1 .9 .4],'visible',options.figure_show);
for clus_id = clusters_of_interest
    id = find(f_clus.id==clus_id); % find the index of all the trajectories with a certain cluster id
%     plotting trajectory of all flights
    subplot(1,numel(clusters_of_interest),count_subplot);
    count_subplot = count_subplot+1;
%     plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
    xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
    xlabel('x');    ylabel('y');    hold on;
    avg_take_off = [];

    for ii=1:size(id,2)
        title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
        plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(clus_id,:));
        avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
    end
    take_off = mean(avg_take_off,2);    
    textscatter(take_off(1),take_off(2),"Take-off");     
    hold off;   axis equal;
end

% figure;

mean_path_3d_all_clusters = {};
mean_len_1d_all_clusters = {};
mean_time_1d = {};
for clus_id = clusters_of_interest
    id = find(f_clus.id==clus_id); % find the index of all the trajectories with a certain cluster id
%     subplot(1,numel(clusters_of_interest),count_subplot);
%     count_subplot = count_subplot+1;
%     xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
%     xlabel('x');    ylabel('y');    hold on;
    
    
    avg_take_off = [];
    %=== initiating xyz positions 
    pos_x_sum = f_clus.pos(1,:,id(1));
    pos_y_sum = f_clus.pos(2,:,id(1));
    pos_z_sum = f_clus.pos(3,:,id(1));
    
    %=== using max len for 1d length 
%     max_len = 0;
%     for ii=1:numel(id)
%        max_len = max(max_len, numel(f_clus.lin_tr{id(ii)})); 
%     end
%     path_1d_sum = zeros(1, max_len);
    %=== using min len for 1d length 
    min_len = numel(f_clus.lin_tr{id(1)});
    for ii=2:numel(id)
       min_len = min(min_len, numel(f_clus.lin_tr{id(ii)})); 
    end
    path_1d_sum = zeros(1, min_len);

    current_path_1d = (f_clus.length(id(1))) * f_clus.lin_tr{id(1)}';
%     path_1d_sum(1:numel(current_path_1d)) = path_1d_sum
%     (1:numel(current_path_1d)) +current_path_1d; % for max length
    path_1d_sum= path_1d_sum +current_path_1d(1:min_len); % for min length
    for ii=2:numel(id)
        %=== finding the mean 3d path
        pos_x_sum = pos_x_sum + f_clus.pos(1,:,id(ii));
        pos_y_sum = pos_y_sum + f_clus.pos(2,:,id(ii));
        pos_z_sum = pos_z_sum + f_clus.pos(3,:,id(ii));
        %=== finding the mean 1d length
        current_path_1d = (f_clus.length(id(ii))) * f_clus.lin_tr{id(ii)}';
%         path_1d_sum(1:numel(current_path_1d)) = path_1d_sum (1:numel(current_path_1d)) +current_path_1d;
        path_1d_sum= path_1d_sum +current_path_1d(1:min_len); % for min length
    end
%     plot3(pos_x_sum/numel(id),pos_y_sum/numel(id),pos_z_sum/numel(id),'-','LineWidth',1,'Color', col_clus(clus_id,:));

    mean_path_3d_all_clusters {clus_id} =  [pos_x_sum/numel(id);pos_y_sum/numel(id);pos_z_sum/numel(id)];
    mean_len_1d_all_clusters {clus_id} = path_1d_sum/numel(id);
    mean_time_1d {clus_id} = [0:1/Fs:min_len/Fs-1/Fs];
%     hold off;   axis equal;
end

% sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
% fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);

end
%% === concatenating spike data spk(time of spike, cell number)
for hide = 1
    if merge_cells == 0
        if isfile("spike_data.mat")
            disp('Loading spike data...');
            load spike_data;
        else
            disp('Concatenating spike data...');
            spk=[];
            for i=1:n_cells %this code takes a couple minutes to run
                for ii=1:size(s{i,1})
                    spk=vertcat(spk, [s{i,1}(ii) i]);
                end
            end
            save('spike_data','spk')
        end
    else 
        if isfile("spike_data_merged.mat")
            disp('Loading spike data...');
            load spike_data_merged;
        else
            disp('Concatenating spike data...');
            spk=[];
            for i=1:n_cells %this code takes a couple minutes to run
                for ii=1:size(s{i,1})
                    spk=vertcat(spk, [s{i,1}(ii) i]);
                end
            end
            save('spike_data_merged','spk')
        end
    end
end
%% === cluster-specific analysis 
disp('Starting replay analysis...');
% declaring variables
for hide = 1
posteriors_all_clusters = {};
posteriers_all_flights = {};
temporal_compression_all_clusters = {};
replay_velocity_all_clusters = {};
weighted_corr_all_clusters = {};
avg_jump_all_clusters = {};
max_jump_all_clusters = {};
replay_score_all_clusters = {};
slope_all_clusters = {};
replay_post_all_clusters = {};
replay_keep_all_clusters = {};


%=== v6 new variables 
flight_type_all_clusters = {};
flight_error_all_clusters = {};
replay_type_all_clusters = {};
flight_active_cells_all_clusters = {};
replay_position_occupancy_all_clusters = {};
replay_rate_all_clusters = {};
replay_times_all_clusters = {};
replay_dist_flight_start = {};
replay_position_all_clusters = {};
replay_segment_len_all_clusters = {};

wc_p_1_all_clusters = {};
wc_p_2_all_clusters = {};

end
for clus_id = clusters_of_interest
disp(['Processing cluster ', num2str(clus_id)])
%% generating place fields
for hide = 1
%=== choosing a flight cluster to analyze 
ids = find(f_clus.id==clus_id);

%=== finding length of all flights
flight_length_mean = mean(f_clus.length(ids));
%=== generate time_position matrix of flights using time and normalized position
time_position = [];
for ii = 1:numel(ids) % loop through each flight within cluster
    distance_from_prev = [];
    st_time = f_clus.strt_frame(ids(ii)); % start time of trajectory
    en_time = f_clus.stop_frame(ids(ii)); % end time of trajectory
    time_stamps = [st_time:en_time];
    time_actual = t(time_stamps);
    positon_normalized = f_clus.lin_tr{ids(ii)};
    time_position = vertcat(time_position, [time_actual positon_normalized]);
end 
% figure; plot(time_position(:,1),time_position(:,2)) % super quick visualization

%=== calculating spike-position data during flight (making matrix of time spike position)
Position_Data = sortrows(time_position,1); % sort by time; note that this is only the position data of the points during the trajectory flight
Spike_Data = sortrows(spk,1); % sort by time; this is the spike data for all spikes

% restricting to time intervals of selected flights 
time_intervals = [];
for ii = 1:numel(ids) % loop through each trajectory within cluster
    st_time = f_clus.strt_frame(ids(ii)); % start time of trajectory
    en_time = f_clus.stop_frame(ids(ii)); % end time of trajectory
    time_intervals = vertcat(time_intervals, [t(st_time), t(en_time)]);
end

Time_Spike_Position=[];
% loop through Spike_Data and assign position data to spike data
for N = 1:size(Spike_Data,1)
    spike_time = Spike_Data(N,1);
    spike_include = 0;
    for i = 1:size(time_intervals)
        if (spike_time >= time_intervals(i,1) && spike_time<=time_intervals(i,2))
            spike_include = 1;
            break
        end
    end
    if (spike_include ==1)
        Index=find(abs(Position_Data(:,1)-Spike_Data(N,1))==min(abs(Position_Data(:,1)-Spike_Data(N,1))),1,'first');
        Time_Spike_Position=vertcat(Time_Spike_Position, [Position_Data(Index,1) Spike_Data(N,2) Position_Data(Index,2)]); % time, spikes, position
        clear Index;
    end

end
%=== splitting the flights into testing and training data 
% training_flights_time_intervals = time_intervals(1:2:end,:);
training_flights_time_intervals = time_intervals; % using all flight to calculate error
testing_flights_time_intervals = time_intervals(2:2:end,:);

%=== getting the flight time, spike, and position data for the training data
time_spike_position_training = [];
for i = 1: size(training_flights_time_intervals,1)
    st_time = training_flights_time_intervals(i,1);
    en_time = training_flights_time_intervals(i,2);
    time_spike_position_training = vertcat(time_spike_position_training, Time_Spike_Position(find(Time_Spike_Position(:,1) >= st_time & Time_Spike_Position(:,1) <= en_time),:));
end
time_spike_position_training = Time_Spike_Position; % when using all the training data 

%=== calculating position distribution
position_training = time_spike_position_training(:,[1 3]);

% generating position bins using average flight distance 
mean_flight_length = mean(f_clus.length(find(f_clus.id==clus_id))); 

%=== using a set number of position bins 
number_position_bin = 30; 
position_bins = [min(position_training(:,2)):1/number_position_bin:max(position_training(:,2))];

%=== using a set length
% position_bin_size = 0.15; % in meters 
% position_bin_size_norm = position_bin_size/mean_flight_length;
% % position_bins = [min(position_training(:,2)):position_bin_size_norm:max(position_training(:,2))];
% position_bins = [0:position_bin_size_norm:1];


if position_bins(end) <1
    position_bins = [position_bins 1];
end
position_distribution = hist(position_training(:,2),position_bins);
position_distribution_smoothed = conv2(position_distribution,normpdf([-10:1:10],0,2),"same");
% figure;plot(position_distribution); xlabel ("Position bins"); ylabel("Occupancy"); title("Position Distribution");hold on;
% plot(position_distribution_smoothed);

%=== calculating number of spikes at a given position
cells_detected = unique(time_spike_position_training(:,2));
spikes_map = zeros(size(s,1), size(position_distribution,2));
for i = cells_detected'
    spikes_map(i,:) = hist(time_spike_position_training(find(time_spike_position_training(:,2)==i),3),position_bins);
end
map_smoothed = conv2(spikes_map,normpdf([-10:1:10],0,2),"same");

%=== calculating average firing rate, number of spikes at a location adjusted by time spent at location
average_firing_rate_training = 100*spikes_map./position_distribution_smoothed;

%=== plotting average firing rate
field_smoothed = conv2(average_firing_rate_training,normpdf([-10:1:10],0,2),"same");
% field_smoothed = average_firing_rate_training; % unsmoothed
[~, field_peak] = max (field_smoothed'); % finding the peak intensity of field
field_peak_indexed = [[1:size(field_smoothed,1)]',field_peak']; % indexing
b = sortrows(field_peak_indexed,2); % sort based on max peak
zfields = zscore(field_smoothed,1,2);
sorted_fields = zfields(b(:,1),:);
fired_cells = sorted_fields(find(b(:,2)>1),:);
figure('visible',options.figure_show);
imagesc(fired_cells)
flight_active_cells_all_clusters{clus_id} = {size(fired_cells,1)};
% title("Flight aligned activity");
ylabel("Cell (sorted by peak)");
xlabel("Position");
colorbar;
field_w = field_smoothed;
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%using only place cells 
% load cell_place_index;
% place_field=zeros(size(field_smoothed)); %using only place cells 
% place_field(cell_place_index,:) = field_smoothed(cell_place_index,:);
% imagesc(zscore(place_field,1,2))
end
%% decoding flight trajectory 
for hide = 1
errors = [];
flight_velocity = [];
field_smoothed=field_w;

% field_smoothed = place_field; % use if using only place cells 

% if using only test trials 
% for current_flight = 1: size(testing_flights_time_intervals,1)
%     flight_st = testing_flights_time_intervals(current_flight,1);
%     flight_en = testing_flights_time_intervals(current_flight,2);

% if using all flights
max_posteriors_all = {};
posterior_all = {};

figure('visible',options.figure_show);set(gcf, 'Position',  [200, 100, 1400, 800]);
plot_index = 1;
for current_flight = 1: size(time_intervals,1)
    if plot_index > 3*4
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('visible',options.figure_show);set(gcf, 'Position',  [200, 100, 1400, 800]);
        plot_index =1;
    end
    flight_st = time_intervals(current_flight,1);
    flight_en = time_intervals(current_flight,2);
    
    tau = 0.05; % each time bin is 100 ms 
    step = tau; % there is no overlap between time bins 
    bin_st = [flight_st:step:flight_en-tau];

    posterior = [];
    for i = 1: numel(bin_st)
%         time_spike_position_flight = Time_Spike_Position(find(Time_Spike_Position >= bin_st(i) & Time_Spike_Position <= bin_st(i)+tau),:);
%         n_i = hist(time_spike_position_flight(:,2), [1:n_cells]);%number of spikes per cell
        spk_current = spk(find(spk(:,1) >= bin_st(i) & spk(:,1) <= bin_st(i)+tau),:);
        n_i = hist(spk_current(:,2), [1:n_cells]);%number of spikes per cell

%     figure;plot(n_i);xlabel("cell id"); ylabel("number of firing during flight");
        prod_term = prod(field_smoothed.^(n_i'));
        exp_term = exp(-tau*(sum(field_smoothed)));
        posterior(i,:) = exp_term.*prod_term;
    end    
    
    posterior_norm = posterior./sum(posterior,2);
    posterior_all{current_flight} = posterior_norm';
    
    % visualization
    subplot(3,4,plot_index);
    imagesc(bin_st,position_bins,posterior_norm'); 
    set(gca,'YDir','normal')
    if plot_index == 1
        title("Decoded position using spike activity");
        xlabel("Flight Time (s)"); ylabel("Flight Position (m)");
    end
    xlim([bin_st(1) bin_st(end)]);
    xticks(bin_st(end));
    xticklabels([num2str((flight_en-flight_st)) ' s']);
    axis xy;
    plot_index =1+plot_index;
end

sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);



% calculating the error from actual position 
figure('visible',options.figure_show);set(gcf, 'Position',  [200, 100, 1400, 800]);
plot_index = 1;
for current_flight = 1: size(time_intervals,1)
    if plot_index > 3*4
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('visible',options.figure_show);set(gcf, 'Position',  [200, 100, 1400, 800]);
        plot_index =1;
    end
    posterior_norm = posterior_all{current_flight};
    flight_st = time_intervals(current_flight,1);
    flight_en = time_intervals(current_flight,2);
    bin_st = [flight_st:step:flight_en-tau];
    
    [~,max_posteriors] = max(posterior_norm, [],1);
    % plotting predicted position (max value)
%     figure; imagesc(posterior_norm);xlabel("time bins"); ylabel("position"); title("Decoded position using spike activity");
%     hold on; scatter ([1:numel(max_posteriors)], max_posteriors)

    % plotting actual position 
    actual_positions = [];
    for i = 1: numel(bin_st)
         time_position_flight = time_position(find(time_position(:,1) >= bin_st(i) & time_position(:,1) <= bin_st(i)+tau),:);
         actual_positions = [actual_positions mean(time_position_flight(:,2))];
%         time_spike_position_flight = Time_Spike_Position(find(Time_Spike_Position(:,1) >= bin_st(i) & Time_Spike_Position(:,1) <= bin_st(i)+tau),:);
%         actual_positions = [actual_positions mean(time_spike_position_flight(:,3))];
    end
    % scale actual position and position bins 
    max_posteriors_scaled = max_posteriors/numel(position_bins);
    % calculate root mean error 
    errors(current_flight) = sqrt(sum((max_posteriors_scaled-actual_positions).^2)/numel(max_posteriors_scaled));

    % visualization
    subplot(3,4,plot_index);
    colormap(hot);

    % visualizing actual vs predicted position 
    scatter (bin_st, actual_positions, "filled")
    hold on; scatter (bin_st, max_posteriors_scaled, "filled");

    xlim([bin_st(1) bin_st(end)]);
    xticks(bin_st(end));
    xticklabels([num2str((flight_en-flight_st)) ' s']);
    axis xy;
    if plot_index == 1
        xlabel("Flight time (s)"); ylabel("Flight Position (m)");
        legend({'Actual position','Decoded position'}, 'Location', "east");
    end
    title(['RMSE: ', num2str(errors(current_flight))]);
%     text(bin_st(1),position_bins(end),['RMSE: ', num2str(errors(current_flight))]);

    %calculating flight velocity
    flight_velocity(current_flight) = (max(actual_positions) - min(actual_positions))/(flight_en-flight_st);
    plot_index =1+plot_index;
end
flight_error_all_clusters{clus_id} = mean(errors);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

end
%% decoding entire session (flights)
for hide = 1
    
disp('Decoding entire session...');
f = field_smoothed;

st=min(spk(:,1));
en=max(spk(:,1));
post_test = [];
f(find(f==0)) = 0.00001;
quantum = 0.05;  % overlap 0.005
tau = 0.05; %decoding window 0.02

nspks=[];
term2=[];
stepnum=0;
count=0;
for i=unique(spk(:,2))'
    count=count+1;
    nspks(count,:)=histc(spk(find(spk(:,2)==i),1),[st:quantum:(en-quantum)]);
end
% nspksC = nspks;
nspksC = conv2(nspks, [1 1 1], "same");
for pos=1:size(f,2)
    term2(pos,:)=prod((f(:,pos)*...
        ones(1,size(nspksC,2))).^nspksC,1);
end
term3=exp(-tau*sum(f,1))' * ones(1,size(term2,2));
num=term2.*term3;
post_test=(num./(ones(size(num,1),1)*sum(num,1)))';
end
posteriers_all_flights {clus_id} = post_test;
%% decoding entire session (replays)
for hide = 1
    
disp('Decoding entire session...');
f = field_smoothed;

st=min(spk(:,1));
en=max(spk(:,1));
post_test = [];
f(find(f==0)) = 0.00001;
quantum = 0.005;  % overlap 0.005
tau = 0.02; %decoding window 0.02

nspks=[];
term2=[];
stepnum=0;
count=0;
for i=unique(spk(:,2))'
    count=count+1;
    nspks(count,:)=histc(spk(find(spk(:,2)==i),1),[st:quantum:(en-quantum)]);
end
% nspksC = nspks;
nspksC = conv2(nspks, [1 1 1], "same");
for pos=1:size(f,2)
    term2(pos,:)=prod((f(:,pos)*...
        ones(1,size(nspksC,2))).^nspksC,1);
end
term3=exp(-tau*sum(f,1))' * ones(1,size(term2,2));
num=term2.*term3;
post_test=(num./(ones(size(num,1),1)*sum(num,1)))';
end
post = post_test;
%% visualization of entire session decoding
for hide = 1
%=== spike density
spk_density_t = [min(spk(:,1)):0.001:max(spk(:,1))];
spk_density = histc(spk(:,1), spk_density_t);
spk_density_smoothed=conv(spk_density,normpdf([-150:150],0,15), 'same');
spk_density_std=std(spk_density_smoothed);
%=== maximum posterior 
[max_prob_all,predicted_position] = max(post'); 
predicted_position=rescale(predicted_position); % converting from position bin to position [0 1]
uniform_prob = 1/size(post_test,2);
predicted_position(find(max_prob_all<uniform_prob*3))=0;

t_NP_max = max(spk(:,1));
t_NP_min = min(spk(:,1));
t_NP_spacing_for_position = (t_NP_max-t_NP_min)/numel(predicted_position);
predicted_position_t = [t_NP_min:t_NP_spacing_for_position:t_NP_max-t_NP_spacing_for_position];
predicted_position_t_less = [t_NP_min:t_NP_spacing_for_position:t_NP_max-2*t_NP_spacing_for_position];
end
%% automatically detecting replay: using posterior spread and posterior center of mass 
for hide = 1
disp('Detecting replay...');
post_spread = []; 
post_COM = [];
post_t = predicted_position_t; % defining the time vector associated with each posterior
position_vector = rescale([1:size(post,2)]);
max_post_spread = 1; 
current_post = [];
%=== calculating posterior spread 
for i =  1:size(post,1) % loop through every time point 
    current_post = post(i,:);
    [~,predicted_pos] = max(current_post'); 
    predicted_pos = position_vector(predicted_pos);
    post_spread(i) = sqrt(sum(((position_vector- predicted_pos).^2).*current_post));
end 
% max_pos_diff = diff(predicted_position);
% max_pos_diff_filtered = max_pos_diff;
% max_pos_diff_filtered(find(max_pos_diff_filtered==0)) = NaN;
subsequent_events = ones(size(post_spread));
subsequent_events(find(post_spread>0.3))=0;
subsequent_events(find(predicted_position==0))=0;
subsequent_events(find(max_prob_all<uniform_prob*3))=0;


%=== filter out events that are horizontal lines 
% predicted_pos_diff = diff(predicted_position);
% predicted_pos_diff = [predicted_pos_diff predicted_pos_diff(end)];
% 
% % calculating the different between the closest value before it while skipping the nans in between; could be helpful to implement eventually  
% % pos_potential_events = predicted_position;
% % pos_potential_events(find(subsequent_events==0))= NaN;
% % max_pos_diff_potential_events = nan(size(predicted_position));
% % max_pos_diff_potential_events(find(isnan(pos_potential_events)==0))=[diff(pos_potential_events(find(isnan(pos_potential_events)==0))),NaN];
% % subsequent_events_strict = subsequent_events;
% % subsequent_events_strict(find(max_pos_diff_potential_events>0.2))=0;
% horiz_perc = 0.010;% if position moved less than this percent in the selected timebins, then erase
% subsequent_events_horiz = zeros(size(post_spread));
% subsequent_events_horiz(find(isnan(predicted_pos_diff)))=1;
% subsequent_events_horiz(abs(predicted_pos_diff)<=horiz_perc)=1;
% % subsequent_events_horiz(strfind(subsequent_events_horiz, [1 0])+1) =1; 
% event_st = strfind(subsequent_events_horiz, [0 1]);
% event_en = strfind(subsequent_events_horiz, [1 0]);
% if subsequent_events_horiz(1)==1
%     event_st = [1 event_st];
% end
% if subsequent_events_horiz(end)==1
%     event_en = [event_en numel(subsequent_events_horiz)];
% end
% event_dur = event_en - event_st;
% max_dur = 4; % erase the events whose diff more than n bins, which means at least n+1 bins, 5*6 = 30ms 
% event_erase_idx = find(event_dur>max_dur); 
% for i = event_erase_idx
%     if (abs(mean(predicted_pos_diff([event_st(i)+1: event_en(i)])))<=horiz_perc*2)
%         subsequent_events([event_st(i): event_en(i)])=0;
%     end
% end
subsequent_events_ori = subsequent_events;

%=== filter out events that are too short
event_st = strfind(subsequent_events, [0 1]);
event_en = strfind(subsequent_events, [1 0]);
if subsequent_events(1)==1
    event_st = [1 event_st];
end
if subsequent_events(end)==1
    event_en = [event_en numel(subsequent_events)];
end
event_dur = event_en - event_st; %hist(event_dur,1000) % visualize event_dur
min_dur = 5; %25 ms
event_erase_idx = find(event_dur<min_dur); 
for i = event_erase_idx
    subsequent_events([event_st(i): event_en(i)])=0;
end

%=== fill in gaps
event_st = strfind(subsequent_events, [0 1]);
event_en = strfind(subsequent_events, [1 0]);
if subsequent_events(1)==1
    event_st = [1 event_st];
end
if subsequent_events(end)==1
    event_en = [event_en numel(subsequent_events)];
end
event_st_en = sort([event_st event_en]);
all_gap = diff(event_st_en);
event_gap = all_gap(2:2:end);
event_dur = all_gap(1:2:end);
gap_max = 15;% fill in gap less than 75ms
event_fill_idx = find(event_gap<gap_max);
for i = event_fill_idx
    subsequent_events([event_en(i):event_en(i)+event_gap(i)+1])=1;
end

%=== erase events that are more than 70% non-noise 
event_st = strfind(subsequent_events, [0 1]);
event_en = strfind(subsequent_events, [1 0]);
if subsequent_events(1)==1
    event_st = [1 event_st];
end
if subsequent_events(end)==1
    event_en = [event_en numel(subsequent_events)];
end
num_bins = event_en-event_st;
event_st_en = sort([event_st event_en]);
all_gap = diff(event_st_en);
event_gap = all_gap(2:2:end);
event_dur = all_gap(1:2:end);
num_zeros = [];
event_erase_idx = [];
for i = 1:numel(event_st)
   num_zeros(i) = numel(strfind(subsequent_events_ori(event_st(i):event_en(i)),0))/num_bins(i);
   if num_zeros(i) >= 0.3
       event_erase_idx = [event_erase_idx i];
   end
end
for i = event_erase_idx
    subsequent_events([event_st(i): event_en(i)])=0;
end

%=== filter out events that are too short  
event_st = strfind(subsequent_events, [0 1]);
event_en = strfind(subsequent_events, [1 0]);
if subsequent_events(1)==1
    event_st = [1 event_st];
end
if subsequent_events(end)==1
    event_en = [event_en numel(subsequent_events)];
end
event_dur = event_en - event_st;
min_dur = 15; %75 ms
event_erase_idx = find(event_dur<min_dur); 
for i = event_erase_idx
    subsequent_events([event_st(i): event_en(i)])=0;
end

%=== restrict replay detection to times when the bat isn't flying 
flight_st = f_clus.strt_frame;
flight_en = f_clus.stop_frame;
flight_st_t = (t(flight_st))';
flight_en_t = (t(flight_en))';
bat_flying = zeros(size(subsequent_events));
for i = 1:numel(flight_st_t)
    bat_flying(find(predicted_position_t>flight_st_t(i) & predicted_position_t<flight_en_t(i)))=1;
    subsequent_events (find(predicted_position_t>flight_st_t(i) & predicted_position_t<flight_en_t(i)))=0;
end
end
%% segmenting candidate events
segment_event = 1;
if segment_event
    for hide = 1
    disp('Selecting candidate events...');
    warning('off');
    start_times = predicted_position_t(strfind(subsequent_events, [0 1]));
    end_times = predicted_position_t(strfind(subsequent_events, [1 0]));
    if subsequent_events(1)==1
        start_times = [predicted_position_t(1) start_times];
    end
    if subsequent_events(end)==1
        end_times = [end_times predicted_position(end)];
    end
    new_st_times = [];
    new_en_times = [];
    for i = 1:numel(start_times)
        %=== not cropping the replays that don't satisfy basic parameters 
%         if abs(weighted_corr_replay(i)) < 0.1 || replay_scores_replay(i) <0.2 || segment_len_frac (i) <0.4
%             new_st_times =  [new_st_times start_times];
%             new_en_times =  [new_en_times end_times];
%             continue;
%         end
%         
        current_times = (find(predicted_position_t >= start_times(i) & predicted_position_t<=(end_times(i))));
        current_posteriors = post(current_times,:);
        [new_st_pos,new_en_pos]= segment_replay_v1(current_posteriors);
        if new_st_pos == 0
            new_st_times =  [new_st_times start_times];
            new_en_times =  [new_en_times end_times];
           continue;
        end
        new_st_times = [new_st_times predicted_position_t(current_times(new_st_pos))];
        new_en_times = [new_en_times predicted_position_t(current_times(new_en_pos))];
    end
    end
else 
    start_times = predicted_position_t(strfind(subsequent_events, [0 1]));
    end_times = predicted_position_t(strfind(subsequent_events, [1 0]));
    if subsequent_events(1)==1
        start_times = [predicted_position_t(1) start_times];
    end
    if subsequent_events(end)==1
        end_times = [end_times predicted_position(end)];
    end
    new_st_times = start_times;
    new_en_times = end_times;
end
%% trimming horizontal parts of candidate events 
for hide = 1
horizontal_replay = zeros(1,numel(new_st_times));
%===for testing purposes
for hide = []
    warning('off');
    for i = 1:numel(new_st_times)
        current_times = (find(predicted_position_t >= new_st_times(i) & predicted_position_t<=(new_en_times(i))));
        if numel(current_times)==1
            continue
        end
        current_posteriors = post(current_times,:);
        [new_st_pos,new_en_pos]= horiz_segment_v0(current_posteriors);
        if isempty(new_st_pos) || isempty(new_en_pos)
            horizontal_replay (i) = 1;
            continue;
        end
        new_st_times (i) = predicted_position_t(current_times(new_st_pos));
        new_en_times (i) = predicted_position_t(current_times(new_en_pos));
    end
end
replay_duration = new_en_times - new_st_times;
end
%% pre/post segmentation visualization 
for hide = []
% top_idx = 1:numel(new_st_times);
% best_idx = abs(weighted_corr_replay)>0.7 & replay_scores_replay > 0.6 ; % index of replays passing certain criteria 
% best_idx = segment_len_frac > 0.1 &abs(weighted_corr_replay)> 0.1 & replay_scores_replay >0.1 ;
best_idx = replay_scores_replay >= 0 &horizontal_replay ~=1;
top_idx = find(best_idx); % only plot events that are valid after segmenting
num_plots = ceil(numel(top_idx)/7/15);

%===sorting top_idx by selected variable
sort_by_variable = weighted_corr_replay;
[~, sorted_index]=sort(abs(sort_by_variable), 'descend');
top_idx=sorted_index(find(ismember(sorted_index,top_idx)));
top_weighted_corr = weighted_corr_replay(top_idx);
top_st_times_ori = start_times(top_idx);
top_en_times_ori = end_times(top_idx);
top_st_times = new_st_times(top_idx);
top_en_times = new_en_times(top_idx);
top_fitted_y = fitted_y_replay(top_idx);
replay_velocity = [];
replay_post_current = {};
nrow = 10;
ncol = 15;
plot_index = 1;
figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
tiledlayout(nrow,ncol,'TileSpacing','tight');
ii=1;
post_t = predicted_position_t;
for i = 1:numel(top_st_times)
    if plot_index > nrow*ncol
%         break
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
        plot_index =1;
        ii=1;
        tiledlayout(nrow,ncol,'TileSpacing','tight');
    end
    if top_st_times(i) <=0 || top_en_times(i) <=0
        continue
    end
    temp_st_times = max(1,find(post_t >= top_st_times_ori(i),1));
    temp_en_times = min(find(post_t <= top_en_times_ori(i),1,'last'), numel(post_t));
    current_times = [temp_st_times:temp_en_times];
    
    event_time_vector = [1:numel(current_times)];
    current_positions = predicted_position(current_times);
    current_posteriors = post_test(current_times,:);
    current_maxprob = max_prob_all(current_times);
    
    ax = nexttile;
    imagesc([post_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
  
    colormap(ax,hot);
    title([i]);
%     title([weighted_corr_replay(top_idx(i))]);
    if plot_index ==1
        ylabel("Position (m)");
        title("candidate replay events");
%         xlabel("Time (s)"); 
    else 
        set(gca, 'YTick', []);
    end
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(top_en_times_ori(i)-top_st_times_ori(i)))) ' ms']);
    axis xy;
    
%=== plotting after segmentation  
    ax = nexttile;  
    temp_st_times = max(1,find(post_t >= top_st_times(i),1));
    temp_en_times = min(find(post_t <= top_en_times(i),1,'last'), numel(post_t));
    current_times = [temp_st_times:temp_en_times];
    
    event_time_vector = [1:numel(current_times)];
    current_positions = predicted_position(current_times);
    current_posteriors = post_test(current_times,:);
    current_maxprob = max_prob_all(current_times);
    imagesc([post_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
  
    
    yticks([]); 
    xticks([]);
%     hold on; area(predicted_position_t,detected_RPL,0,'FaceColor',[0 0 0],'FaceAlpha',0.5,'LineStyle','none');  xlim([top_st_times(i),top_en_times(i)]);  yticks([]);

    ii=ii+1;
    colormap(hot)
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
    axis xy;
    plot_index = plot_index +2;
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% calculating metrics using cropped replay times
for hide = 1
% variables for storing replay metrics 
weighted_corr_replay = [];
weighted_corr_p_replay = [];
avg_jump_distance_replay  = []; 
max_jump_distance_replay=[];
replay_scores_replay =  [];
slope_replay = [];
num_peaks_replay = [];
fitted_y_replay = {};
posterior_spread_replay = [];
com_jump_replay = [];
segment_len_frac = [];
percent_horizontal = [];

for i = 1:numel(new_st_times)
    if new_st_times <=0
        % variables for storing replay metrics 
        weighted_corr_replay(i) = NaN;
        weighted_corr_p_replay (i,:) = NaN;
        avg_jump_distance_replay(i) = NaN;
        max_jump_distance_replay(i) = NaN;
        replay_scores_replay(i) = NaN;
        slope_replay(i) = NaN;
        num_peaks_replay(i) = NaN;
        fitted_y_replay{1}= NaN ;
        posterior_spread_replay(i) = NaN;
        com_jump_replay(i) = NaN;
        segment_len_frac(i) = NaN;
        percent_horizontal(i) = NaN;
    else
    current_posteriors = post((find(predicted_position_t >= new_st_times(i) & predicted_position_t<=(new_en_times(i)))),:);
    [weighted_corr_replay(i),max_jump_distance_replay(i), avg_jump_distance_replay(i), num_peaks_replay(i), replay_scores_replay(i), slope_replay(i),fitted_y_replay{i}, posterior_spread_replay(i), com_jump_replay(i),segment_len_frac(i),percent_horizontal(i)] = evaluate_candidate_event_v6(current_posteriors);
    [weighted_corr_p_replay(i,:), ~, ~, ~, ~, ~, ~,~] = shuffle_validation_v2(current_posteriors);
end
end 
end
%% selecting for good events 
for hide = 1
wc_p_1 = weighted_corr_p_replay(:,1)';
wc_p_2 = weighted_corr_p_replay(:,2)';
replay_padding_bins = 40; %add 20 bins on each side - total 200 ms
selected_replays = find(new_st_times>quantum*replay_padding_bins & new_en_times<max(spk(:,1))-replay_padding_bins & replay_duration > 0.05);
new_st_times = new_st_times(selected_replays);
new_en_times = new_en_times(selected_replays);
weighted_corr_replay = weighted_corr_replay(selected_replays);
avg_jump_distance_replay = avg_jump_distance_replay(selected_replays);
max_jump_distance_replay = max_jump_distance_replay(selected_replays);
replay_scores_replay = replay_scores_replay(selected_replays);
slope_replay = slope_replay(selected_replays);
num_peaks_replay = num_peaks_replay(selected_replays);
fitted_y_replay= fitted_y_replay(selected_replays);
posterior_spread_replay = posterior_spread_replay(selected_replays);
com_jump_replay = com_jump_replay(selected_replays);
segment_len_frac = segment_len_frac(selected_replays);
percent_horizontal = percent_horizontal(selected_replays);
wc_p_1 = wc_p_1(selected_replays);
wc_p_2 = wc_p_2(selected_replays);
end
%% visualizing detected events
for hide = []
    
%=== selecting index of all significant replay
top_idx = find(new_st_times);
num_plots = ceil(numel(top_idx)/7/15);

%===sorting top_idx by selected variable
sort_by_variable = weighted_corr_replay;
[~, sorted_index]=sort(abs(sort_by_variable), 'descend');
top_idx=sorted_index(find(ismember(sorted_index,top_idx)));
top_weighted_corr = weighted_corr_replay(top_idx);
top_st_times = new_st_times(top_idx);
top_en_times = new_en_times(top_idx);
top_fitted_y = fitted_y_replay(top_idx);
replay_velocity = [];
replay_post_current = {};
plot_index = 1;

figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
tiledlayout(15,20,'TileSpacing','tight');
ii=1;
for i = 1:numel(top_st_times)
    if plot_index > 20*15
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
        plot_index =1;
        ii=1;
        tiledlayout(15,20,'TileSpacing','tight');
    end
    nexttile(40*floor((ii-1)/20)+ii,[2 1]);
    current_times = find(predicted_position_t >= top_st_times(i) & predicted_position_t <= top_en_times(i));  
    event_time_vector = [1:numel(current_times)];
    current_positions = predicted_position(current_times);
    current_posteriors = post(current_times,:);
    current_maxprob = max_prob_all(current_times);

    
    
    
%     calculating replay speed 
    current_times_s = predicted_position_t(current_times); % converting from time bins to seconds
    max_posterior_time_pos =[];
    max_posterior_time_pos(:,1) = current_times_s;
    max_posterior_time_pos(:,2) = current_positions;
    max_posterior_time_pos_filtered = max_posterior_time_pos;
    
    % plotting line used in evaluation of replay 
    colormap(hot);
    imagesc([predicted_position_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;

    % calculating replay velocity 
    robust_fit = robustfit(max_posterior_time_pos_filtered(:,1),max_posterior_time_pos_filtered(:,2));
    line_slope_intercept = [robust_fit(2) robust_fit(1)]; 
    line_f = polyval(line_slope_intercept,current_times_s);
    replay_velocity(i) = (line_f(end)-line_f(1))*flight_length_mean/(current_times_s(end)-current_times_s(1));
    replay_post_current{i} = current_posteriors;
    hold on; plot(current_times_s,line_f,'-','Color',[1 1 1],'LineWidth',1)
    title([i]);
    if plot_index ==1
        ylabel("Position (m)");
        title("candidate replay events");
    else 
        set(gca, 'YTick', []);
    end
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
    axis xy;
    
%=== plotting LFP trace 
    nexttile(40*(1+floor((ii-1)/20))+ii);  
    current_times_lfp = find(all_RPL(RPL_max_idx).RPL_out.time_vector  >= top_st_times(i) &all_RPL(RPL_max_idx).RPL_out.time_vector  <= top_en_times(i)) ;
    plot(all_RPL(RPL_max_idx).RPL_out.time_vector(current_times_lfp),LFP_rp(current_times_lfp),'b'); 
    yticks([]); 
    xticks([]);
    ii=ii+1;  
    plot_index = plot_index +3;
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% visualizing detected events vs include/exclude (w/ LFP - very slow)
for hide = []
    
%=== selecting index of all significant replay
% top_idx = 1:numel(new_st_times);
top_idx = find(new_st_times); % only plot events that are valid after segmenting
best_idx = abs(weighted_corr_replay)>0.8  & replay_scores_replay >0.4 & new_st_times; % index of replays passing certain criteria 
num_plots = ceil(numel(top_idx)/7/15);

%===sorting top_idx by selected variable
sort_by_variable = weighted_corr_replay;
[~, sorted_index]=sort(abs(sort_by_variable), 'descend');
top_idx=sorted_index(find(ismember(sorted_index,top_idx)));
top_weighted_corr = weighted_corr_replay(top_idx);
top_st_times = new_st_times(top_idx);
top_en_times = new_en_times(top_idx);
top_fitted_y = fitted_y_replay(top_idx);
replay_velocity = [];
replay_post_current = {};
plot_index = 1;
figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
tiledlayout(15,20,'TileSpacing','tight');
ii=1;
for i = 1:numel(top_st_times)
% for i = 500:599
    if plot_index > 20*15
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
        plot_index =1;
        ii=1;
        tiledlayout(15,20,'TileSpacing','tight');
    end
    ax = nexttile(40*floor((ii-1)/20)+ii,[2 1]);
    
    
    current_times = find(predicted_position_t >= top_st_times(i) & predicted_position_t <= top_en_times(i));  
    event_time_vector = [1:numel(current_times)];
    current_positions = predicted_position(current_times);
    current_posteriors = post(current_times,:);
    current_maxprob = max_prob_all(current_times);

%     calculating replay speed 
    current_times_s = predicted_position_t(current_times); % converting from time bins to seconds
    max_posterior_time_pos =[];
    max_posterior_time_pos(:,1) = current_times_s;
    max_posterior_time_pos(:,2) = current_positions;
    max_posterior_time_pos_filtered = max_posterior_time_pos;
    
    % plotting line used in evaluation of replay 

    imagesc([predicted_position_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
  
    if  best_idx(top_idx(i))==1
        colormap(ax,hot);
    else
        colormap(ax, parula);
    end
    % calculating replay velocity 
    robust_fit = robustfit(max_posterior_time_pos_filtered(:,1),max_posterior_time_pos_filtered(:,2));
    line_slope_intercept = [robust_fit(2) robust_fit(1)]; 
    line_f = polyval(line_slope_intercept,current_times_s);
    replay_velocity(i) = (line_f(end)-line_f(1))*flight_length_mean/(current_times_s(end)-current_times_s(1));
    replay_post_current{i} = current_posteriors;
    hold on; plot(current_times_s,line_f,'-','Color',[1 1 1],'LineWidth',1)
 
%     title(["slope: ",slope_replay(find(start_times== top_st_times(i)))]);
%     title([i]);
    title([weighted_corr_replay(top_idx(i))]);
    if plot_index ==1
        ylabel("Position (m)");
        title("candidate replay events");
%         xlabel("Time (s)"); 
    else 
        set(gca, 'YTick', []);
    end
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
    axis xy;
    
%=== plotting LFP trace 
    nexttile(40*(1+floor((ii-1)/20))+ii);  
    current_times_lfp = find(all_RPL(RPL_max_idx).RPL_out.time_vector  >= top_st_times(i) &all_RPL(RPL_max_idx).RPL_out.time_vector  <= top_en_times(i)) ;
    if  best_idx(top_idx(i))==1
        plot(all_RPL(RPL_max_idx).RPL_out.time_vector(current_times_lfp),LFP_rp(current_times_lfp),'r'); 
    else
        plot(all_RPL(RPL_max_idx).RPL_out.time_vector(current_times_lfp),LFP_rp(current_times_lfp),'b'); 
    end

    yticks([]); 
    xticks([]);
%     hold on; area(predicted_position_t,detected_RPL,0,'FaceColor',[0 0 0],'FaceAlpha',0.5,'LineStyle','none');  xlim([top_st_times(i),top_en_times(i)]);  yticks([]);

    ii=ii+1;

    
    plot_index = plot_index +3;
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% visualizing detected replays (w/ LFP) & calculating replay speed
for hide = 1
% % adjusting time vector for the visualization of test decoding 
% [max_prob_all,predicted_position] = max(post_test'); 
% predicted_position=rescale(predicted_position); % converting from position bin to position [0 1]
% uniform_prob = 1/size(post_test,2);
% predicted_position(find(max_prob_all<uniform_prob*3))=0;
% t_NP_max = max(spk(:,1));
% t_NP_min = min(spk(:,1));
% t_NP_spacing_for_position = (t_NP_max-t_NP_min)/numel(predicted_position);
% predicted_position_t_test = [t_NP_min:t_NP_spacing_for_position:t_NP_max-t_NP_spacing_for_position];
% predicted_position_t_test = predicted_position_t;
%=== TESTING selecting index of all significant replay
% top_idx = 1:numel(new_st_times);
%=== test
% best_idx = replay_scores_replay > 0.5 & segment_len_frac >0.5 & abs(weighted_corr_replay) > 0.4 & percent_horizontal<0.5 ; % index of replays passing certain criteria 
%=== actual run 
best_idx = abs(weighted_corr_replay) >-1; % for saving all replays

top_idx = find(best_idx); % only plot events that are valid after segmenting
num_plots = ceil(numel(top_idx)/7/15);

%===sorting top_idx by selected variable
% sort_by_variable = weighted_corr_replay;
% [~, sorted_index]=sort(abs(sort_by_variable), 'descend');
% top_idx=sorted_index(find(ismember(sorted_index,top_idx)));
% top_weighted_corr = weighted_corr_replay(top_idx);
top_st_times = new_st_times(top_idx);
top_en_times = new_en_times(top_idx);
% top_fitted_y = fitted_y_replay(top_idx);
plot_index = 1;
figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
tiledlayout(15,20,'TileSpacing','tight');
ii=1;
post_t = predicted_position_t;
replay_velocity = [];
replay_post_current = {};
% post_t = predicted_position_t_test;
for i = 1:numel(top_st_times)
% for i = 1:20*15
    if plot_index > 20*15
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'visible',options.figure_show);
        plot_index =1;
        ii=1;
        tiledlayout(15,20,'TileSpacing','tight');
    end
    ax = nexttile(40*floor((ii-1)/20)+ii,[2 1]);
    if top_st_times(i) <=0+replay_padding_bins || top_en_times(i) >=max(spk(:,1))-replay_padding_bins
        continue
    end
%     temp_st_times = max(1,find(post_t >= top_st_times(i),1));
%     temp_en_times = min(find(post_t <= top_en_times(i),1,'last'), numel(post_t));
%     

    temp_st_times = max(1,find(post_t >= top_st_times(i),1))-replay_padding_bins;
    temp_en_times = min(find(post_t <= top_en_times(i),1,'last'), numel(post_t))+replay_padding_bins;
    current_times = [temp_st_times:temp_en_times];
    
    event_time_vector = [1:numel(current_times)];
    current_positions = predicted_position(current_times);
    current_posteriors = post_test(current_times,:);
    current_maxprob = max_prob_all(current_times);

%     calculating replay speed 
    current_times_s = post_t(current_times); % converting from time bins to seconds
    max_posterior_time_pos =[];
    max_posterior_time_pos(:,1) = current_times_s;
    max_posterior_time_pos(:,2) = current_positions;
    max_posterior_time_pos_filtered = max_posterior_time_pos;
    
    % plotting line used in evaluation of replay 

    imagesc([post_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
  
    if  best_idx(top_idx(i))==1
        colormap(ax,hot);
    else
        colormap(ax, parula);
    end
    % calculating replay velocity 
    robust_fit = robustfit(max_posterior_time_pos_filtered(:,1),max_posterior_time_pos_filtered(:,2));
    line_slope_intercept = [robust_fit(2) robust_fit(1)]; 
    line_f = polyval(line_slope_intercept,current_times_s);
    replay_velocity(i) = (line_f(end)-line_f(1))*flight_length_mean/(current_times_s(end)-current_times_s(1));
    replay_post_current{i} = current_posteriors;
    hold on; plot(current_times_s,line_f,'-','Color',[1 1 1],'LineWidth',1)
 
    title([i]);
%     title([weighted_corr_replay(top_idx(i))]);
    if plot_index ==1
        ylabel("Position (m)");
        title("candidate replay events");
%         xlabel("Time (s)"); 
    else 
        set(gca, 'YTick', []);
    end
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
    axis xy;
    
%=== plotting LFP trace 
    nexttile(40*(1+floor((ii-1)/20))+ii);  
    current_times_lfp = find(all_RPL(RPL_max_idx).RPL_out.time_vector  >= top_st_times(i) &all_RPL(RPL_max_idx).RPL_out.time_vector  <= top_en_times(i)) ;
    if  best_idx(top_idx(i))==1
        plot(all_RPL(RPL_max_idx).RPL_out.time_vector(current_times_lfp),LFP_rp(current_times_lfp),'r'); 
    else
        plot(all_RPL(RPL_max_idx).RPL_out.time_vector(current_times_lfp),LFP_rp(current_times_lfp),'b'); 
    end

    yticks([]); 
    xticks([]);
%     hold on; area(predicted_position_t,detected_RPL,0,'FaceColor',[0 0 0],'FaceAlpha',0.5,'LineStyle','none');  xlim([top_st_times(i),top_en_times(i)]);  yticks([]);

    ii=ii+1;

    
    plot_index = plot_index +3;
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% visualize entire session with detected replay
%=== visualization
for hide = []
    % variable to visualize replay
    detected_replays = zeros(size(post_spread));
    for i = 1:numel(top_idx)
      detected_replays (find(predicted_position_t>top_st_times(i) & predicted_position_t<top_en_times(i)))=1;
    end
    % variable to visualize ripple times 
    detected_RPL = zeros(size(post_spread));
    
    for i = 1:numel(RPL_t_all)
        [~, min_idx] = min(abs(RPL_t_all (i)-predicted_position_t));
        detected_RPL(min_idx) = max(LFP_rp);
    end
    
    
    figure('visible',options.figure_show);
    colormap(hot)
    tiledlayout(4,1)
    ax1 = nexttile;
%     imagesc([t_NP_min, t_NP_max],[0,max(Position_Data(:,2))],post');xlabel("Time"); ylabel("Position bins"); title("posterior_raw"); hold on;
%     scatter(predicted_position_t, predicted_position, 5,"filled","w");
    imagesc([t_NP_min, t_NP_max],[0,max(Position_Data(:,2))],post',prctile(post,[5 99.5],"all")');xlabel("Time"); ylabel("Position bins"); title("posteior");
    axis xy;
    hold on; area(predicted_position_t,detected_replays,0,'FaceColor',[1 1 1],'FaceAlpha',0.5,'LineStyle','none'); 
    axis xy;
    ax2 = nexttile;
    plot(all_RPL(RPL_max_idx).RPL_out.time_vector,LFP_rp)
    hold on; area(predicted_position_t,detected_RPL,0,'FaceColor',[0 0 0],'FaceAlpha',0.5,'LineStyle','none'); 

    ax3 = nexttile;
    plot(t, r(:,1)); title ("position")
    ax4 = nexttile;
    plot(spk_density_t,spk_density_smoothed); title ("spike density");

    linkaxes([ax1,ax2, ax3, ax4],'x')
end
%% calculating & visualizing temporal compression of replay
for hide = 1
bin_size = 25;
flight_velocity_avg = mean(flight_velocity)*flight_length_mean; 
velocity_compression = replay_velocity/flight_velocity_avg;
% figure;hist(abs(replay_velocity),38); ylabel("Number of events"); xlabel("Replay velocity (m/s)"); % title("Replay speed"); 
% sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');

figure('visible',options.figure_show);hist(abs(velocity_compression),bin_size); ylabel("Number of events"); xlabel("Temporal Compression amount (fold)"); %title("Replay speed / flight speed"); 
xlim([0 50]);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

% plotting time compression of forward vs reverse replay 
figure('visible',options.figure_show);
fwd_replay=find(velocity_compression>0);
rev_replay=find(velocity_compression<0);
% scatter(weighted_corr_replay(top_idx), abs(velocity_compression))
% scatter(abs(weighted_corr_replay(top_idx)), abs(velocity_compression))
scatter(abs(weighted_corr_replay(top_idx(fwd_replay))), abs(velocity_compression(fwd_replay)),"filled","b"); hold on;
scatter(abs(weighted_corr_replay(top_idx(rev_replay))), abs(velocity_compression(rev_replay)),"filled","r");
ylabel("time compression (absolute value)");
xlabel("weighted correlation (absolute value)");
legend("forward replay", "reverse replay");
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
figure('visible',options.figure_show);
histogram(abs(velocity_compression(fwd_replay)),bin_size,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');    hold on;
histogram(abs(velocity_compression(rev_replay)),bin_size,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
legend("Forward Replay", "Reverse Replay")
ylabel ("Number of events (normalized)")
xlabel("Temporal compression (fold)")
xlim([0 50]);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% plotting replays with specific values of velocity_compression 
for hide = []
% current_groups = {};
% % current_groups{1} = top_idx(find(velocity_compression >0));
% % current_groups{2}= top_idx(find(velocity_compression <0));
% 
% current_groups{1} = top_idx(find(abs(velocity_compression)>0 & abs(velocity_compression)<5 ));
% current_groups{2}= top_idx(find(abs(velocity_compression)>5 & abs(velocity_compression)<15 ));
% current_groups{3}= top_idx(find(abs(velocity_compression)>15 ));
% 
% for ii = 1:numel(current_groups)
% current_group = current_groups {ii};
% %=== selecting index of all significant replay
% num_plots = ceil(numel(current_group)/7/15);
% 
% % top_idx = find(abs(slope_replay)>0.01 & abs(slope_replay)<0.2);
% top_st_times = start_times(current_group);
% top_en_times = end_times(current_group);
% 
% top_fitted_y = fitted_y_replay(current_group);
% % sorted_posteriors = post(index);
% 
% figure('units','normalized','outerposition',[0 .05 1 0.95]);
% sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id),', group ',num2str(ii)],'Interpreter','none');
% 
% plot_index = 1;
% for i = 1:numel(top_st_times)
%     if plot_index > 7*15
%         sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id),', group ',num2str(ii)],'Interpreter','none');
%         fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
%         figure('units','normalized','outerposition',[0 .05 1 0.95]);
%         plot_index =1;
%     end
%     subplot(7,15,plot_index);
%     current_times = find(predicted_position_t >= top_st_times(i) & predicted_position_t <= top_en_times(i));  
%     event_time_vector = [1:numel(current_times)];
%     current_positions = predicted_position(current_times);
%     current_posteriors = post(current_times,:);
%     current_maxprob = max_prob_all(current_times);
%     
% %     calculating replay speed 
%     current_times_s = predicted_position_t(current_times); % converting from time bins to seconds
% 
%     max_posterior_time_pos =[];
%     max_posterior_time_pos(:,1) = current_times_s;
%     max_posterior_time_pos(:,2) = current_positions;
%     max_posterior_time_pos_filtered = max_posterior_time_pos(find(current_maxprob' > (3*uniform_prob) & max_posterior_time_pos(:,2) > 0 & max_posterior_time_pos(:,2) < max(max_posterior_time_pos(:,2))),:);
% 
%     colormap(hot);
%     imagesc([predicted_position_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
% 
%     % plotting line with polyfit
% %     hold on; plot([predicted_position_t(current_times)],top_fitted_y{i},'-','Color',[1 1 1])
% %     hold on; plot([predicted_position_t(current_times)],line_f,'-','Color',[1 1 1])
%     % plotting line with lsline     
% %     hold on; scatter(predicted_position_t(current_times), current_positions,7,"filled",'w');
% %     hold on; lsline;
% 
% %     title(["slope: ",slope_replay(find(start_times== top_st_times(i)))]);
% %     title(["i: ",i]);
% %     title(["wc: ", weighted_corr_replay(find(start_times== top_st_times(i)))]);
%     if plot_index ==1
%         ylabel("Position (m)");
%         title("candidate replay events");
%         xlabel("Time (s)"); 
%     else 
%         set(gca, 'YTick', []);
%     end
%     xticks(predicted_position_t(current_times(end)));
%     xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
%     axis xy;
%     plot_index = plot_index +1;
% end
% end
end
%% plotting replays in order of velocity_compression 
for hide = []
% %===sorting top_idx by selected variable
% sort_by_variable = abs(velocity_compression);
% % sort_by_variable = slope_replay;
% [~, sorted_index]=sort(abs(sort_by_variable), 'descend');
% top_idx=sorted_index(find(ismember(sorted_index,top_idx)));
% 
% 
% 
% % top_idx = find(abs(slope_replay)>0.01 & abs(slope_replay)<0.2);
% top_st_times = start_times(top_idx);
% top_en_times = end_times(top_idx);
% 
% top_fitted_y = fitted_y_replay(top_idx);
% % sorted_posteriors = post(index);
% 
% figure('units','normalized','outerposition',[0 .05 1 0.95]);
% plot_index = 1;
% for i = 1:numel(top_st_times)
%     if plot_index > 7*15
%         sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
%         fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
%         figure('units','normalized','outerposition',[0 .05 1 0.95]);
%         plot_index =1;
%     end
%     subplot(7,15,plot_index);
%     current_times = find(predicted_position_t >= top_st_times(i) & predicted_position_t <= top_en_times(i));  
%     event_time_vector = [1:numel(current_times)];
%     current_positions = predicted_position(current_times);
%     current_posteriors = post(current_times,:);
%     
% %     calculating replay speed 
%     current_times_s = predicted_position_t(current_times); % converting from time bins to seconds
% 
%     max_posterior_time_pos =[];
%     max_posterior_time_pos(:,1) = current_times_s;
%     max_posterior_time_pos(:,2) = current_positions;
%     max_posterior_time_pos_filtered = max_posterior_time_pos(find(max_posterior_time_pos(:,2) > 0 & max_posterior_time_pos(:,2) < max(max_posterior_time_pos(:,2))),:);
%     colormap(hot);
%     imagesc([predicted_position_t(current_times)],[0,max(Position_Data(:,2))],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
% 
% 
%     hold on; plot([predicted_position_t(current_times)],top_fitted_y{i},'-','Color',[1 1 1])
%     
%     % plotting line with lsline     
% %     hold on; scatter(predicted_position_t(current_times), current_positions,7,"filled",'w');
% %     hold on; lsline;
% 
% %     title(["slope: ",slope_replay(find(start_times== top_st_times(i)))]);
% %     title(["i: ",i]);
% %     title(["wc: ", weighted_corr_replay(find(start_times== top_st_times(i)))]);
%     if plot_index ==1
%         ylabel("Position (m)");
%         title("candidate replay events");
%         xlabel("Time (s)"); 
%     else 
%         set(gca, 'YTick', []);
%     end
%     xticks(predicted_position_t(current_times(end)));
%     xticklabels([num2str(round(1000*(top_en_times(i)-top_st_times(i)))) ' ms']);
%     axis xy;
%     plot_index = plot_index +1;
% end
% sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
% fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
% 
end
%% saving whole session variables
for hide = 1
posteriors_all_clusters{clus_id} = post;
replay_times_all_clusters{clus_id} = [top_st_times' top_en_times'];
temporal_compression_all_clusters{clus_id} = velocity_compression;
replay_velocity_all_clusters{clus_id} = replay_velocity;
weighted_corr_all_clusters{clus_id} = weighted_corr_replay(top_idx);
avg_jump_all_clusters{clus_id} = avg_jump_distance_replay(top_idx); 
max_jump_all_clusters{clus_id} = max_jump_distance_replay(top_idx);
replay_score_all_clusters{clus_id} = replay_scores_replay(top_idx);
slope_all_clusters{clus_id} = slope_replay(top_idx);
replay_post_all_clusters{clus_id} = replay_post_current;
end
%% analyzing the location and timing of replay 
for hide = 1
%% === Analysis of type of replay
for hide = 1
figure('visible',options.figure_show);
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64];                                       % Feeder coordinates
if strcmp(unique_ID{1,1},'Dataset_2')
    r_fd = [-3,0.82,2; 100,100,100];                                            % Correct coordinates for Field Station
end
X = r(~bflying & v_abs<0.5,:);                                                  % Consider only stationary epochs
[~,centroid,~,~,~] = Cluster3D_AF_v1(X,30*Fs,0.3,100);                          % Perform clustering
for i=1:2                                                                       % Correct feeder position(s)
    [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid,2,2));
    if feeder_distance<0.3,r_fd(i,:) =  centroid(fd_idx,:);end                  % Do not correct if further than 30 cm
end
bat_on_feeder = vecnorm(r_fd(1,:)-r,2,2)<0.3 | vecnorm(r_fd(2,:)-r,2,2)<0.3;    % Bat on feeder
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== finding fraction of replays on feeder vs not 
bat_on_feeder_t = t(bat_on_feeder);
replay_on_feeder =  any (abs(top_st_times-bat_on_feeder_t)<0.01);
figure('visible',options.figure_show);histogram(replay_on_feeder, 2,'Normalization','probability'); xlabel("On Feeder"); ylabel("Fraction");
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== finding bat's position during replay 
pos_during_replay =  [];
for i = 1: numel(top_st_times)
    [~, min_idx] = min(abs(abs(top_st_times(i)-t)));
    pos_during_replay = vertcat(pos_during_replay, r(min_idx,:));     
%     disp([i size(pos_during_replay)])
end
end
%% find the average start/end position of flight 
for hide = 1 
pos_st = [];
flight_st = f_clus.strt_frame(ids);
flight_en = f_clus.stop_frame(ids);
for i = 1: numel(flight_st)
    pos_st = vertcat(pos_st, r(flight_st(i),:)); 
end
mean_takeoff = mean(pos_st);

pos_en = [];
for i = 1: numel(flight_st)
    pos_en = vertcat(pos_en, r(flight_en(i),:)); 
end
mean_landing = mean(pos_en);

% visualizing take off, landing positions, and mean takeoff/landing positions 
% figure;
% scatter3(pos_st(:,1),pos_st(:,2),pos_st(:,3),'filled'); 
% hold on; scatter3(pos_en(:,1),pos_en(:,2),pos_en(:,3),'filled'); alpha(0.5);
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled');
% xlabel ("x"); ylabel ("y"); zlabel ("z");
% title ('take off and landing positions');
% legend (["take off", "landing", "mean take off", "mean landing"])

%=== visualizing all replay positions
% figure;
% scatter3(pos_during_replay(:,1),pos_during_replay(:,2),pos_during_replay(:,3),5);
% hold on; 
% scatter3(pos_st(:,1),pos_st(:,2),pos_st(:,3),'filled'); 
% hold on; scatter3(pos_en(:,1),pos_en(:,2),pos_en(:,3),'filled');
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled');
% xlabel ("x"); ylabel ("y"); zlabel ("z");
% title('position of bat during replay');
% legend (["replay position","take off", "landing"], 'Location','best');
if isempty(pos_during_replay)
    continue
end
end
%% clustering replays based on their positions 
for hide = []
%=== putting replays into clusters 
X = pos_during_replay;
replay_clus = dbscan(X, 0.5,1)';
replay_clus_uni = unique(replay_clus);
replay_clus_centroid = [];

%=== visualizing the clusters
% figure;
for i = replay_clus_uni
    current_clus = find(replay_clus==i);
    if i ~= -1 
        if numel(current_clus)==1
            current_centroid = pos_during_replay(current_clus,:);
        else
            current_centroid = mean(pos_during_replay(current_clus,:));
        end
        replay_clus_centroid ([i],:) = current_centroid;
%         hold on;scatter3(current_centroid(1),current_centroid(2),current_centroid(3),'filled','DisplayName',['center ', num2str(i)]);

    end
%     hold on;scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'DisplayName',['cluster ', num2str(i)]);
%     max(vecnorm((pos_during_replay(current_clus,:)-current_centroid)')); % calculating
%     the max distance between each point and centroid 
end
% legend;
% xlabel ("x"); ylabel ("y"); zlabel ("z");
% title('replays colored by position');
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing');
% xlabel ("x"); ylabel ("y"); zlabel ("z");

%%=== find out if the replays are happening at feeders, landing, takeoff, or remote locations 
%=== defining important coordinates 
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64];                                       % Feeder coordinates
if strcmp(unique_ID{1,1},'Dataset_2')
    r_fd = [-3,0.82,2; 100,100,100];                                            % Correct coordinates for Field Station
end
replay_clus_centroid;
replay_clus_labels = cell(1,size(replay_clus_centroid,1));

%=== labelling each centroid as takeoff/landing/remote/feeders 
distance_tolerance = 0.6;
for i = 1:size(replay_clus_centroid,1)
    current_centroid = replay_clus_centroid(i,:);
    for ii = 1: size(r_fd,1)
        if vecnorm(current_centroid-r_fd(ii,:),2,2)<distance_tolerance
            replay_clus_labels{i} = [replay_clus_labels{i} 'feeder '];
        end
    end
    if vecnorm(current_centroid-mean_takeoff,2,2)<distance_tolerance 
        replay_clus_labels{i} = [replay_clus_labels{i} 'takeoff '];
    end
    if vecnorm(current_centroid-mean_landing,2,2)<distance_tolerance 
        replay_clus_labels{i} = [replay_clus_labels{i} 'landing '];
    end
    if isempty(replay_clus_labels{i})
        replay_clus_labels{i} = ['remote'];
    end


%=== visualizing the clusters with new labels 
figure('visible',options.figure_show);
for i = replay_clus_uni
    current_clus = find(replay_clus==i);
    if i == -1 
        hold on; scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'k','DisplayName',['remote']);

    else
        hold on;scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'DisplayName',[replay_clus_labels{i}]);
    end
end
legend;
xlabel ("x"); ylabel ("y"); zlabel ("z");
title('replays colored by position');
hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off');
hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing');
end

%%=== reassign clusters & visualize

replay_clus_char = cell(size(current_clus));
for i = replay_clus_uni
    current_clus =find(replay_clus==i);
    if i == -1
        for ii = current_clus
            replay_clus_char{ii} = 'remote';
        end
    else
       for ii = current_clus
            replay_clus_char{ii} =  replay_clus_labels{i};
       end
    end
end

figure('visible',options.figure_show);
unique_clusters = unique(replay_clus_char);
for i = 1:numel(unique_clusters)
    current_clus =  find(strcmpi(replay_clus_char, {unique_clusters{i}}));
    hold on;scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'DisplayName',[unique_clusters{i}]);
end

legend;
xlabel ("x"); ylabel ("y"); zlabel ("z");
title('replays colored by position');
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing');

id = find(f_clus.id==clus_id);
plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001,'HandleVisibility', 'off');
xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
xlabel('x');    ylabel('y');    hold on;
avg_take_off = [];

for ii=1:size(id,2)
    title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
    plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(clus_id,:),'HandleVisibility', 'off');
    avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
end
take_off = mean(avg_take_off,2);    

textscatter(take_off(1),take_off(2),"Take-off",'HandleVisibility', 'off');     
hold off;   axis equal;
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

replay_type_all_clusters{clus_id} = replay_clus_char;
end
%% determining type of replay based on distance to feeders, take-off, and landing 
%=== defining important coordinates 
% mean_takeoff, mean_landing 
% pos_during_replay(current_clus,:)
r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64];                                       % Feeder coordinates
if strcmp(unique_ID{1,1},'Dataset_2')
    r_fd = [-3,0.82,2; 100,100,100];                                            % Correct coordinates for Field Station
end

replay_clus_char = cell(size(top_st_times));

%=== labelling each replay as takeoff/landing/remote/feeders 
distance_tolerance = 0.5;
for i = 1:numel(top_st_times)
    current_position = pos_during_replay(i,:);
%     for ii = 1: size(r_fd,1)
%         if vecnorm(current_position-r_fd(ii,:),2,2)<distance_tolerance
%             replay_clus_char{i} = [replay_clus_char{i} 'feeder '];
%         end
%     end
    if vecnorm(current_position-mean_takeoff,2,2)<distance_tolerance 
        replay_clus_char{i} = ['takeoff'];
    end
    if vecnorm(current_position-mean_landing,2,2)<distance_tolerance 
        replay_clus_char{i} = ['landing'];
    end
    if isempty(replay_clus_char{i})
        replay_clus_char{i} = ['remote'];
    end
end

replay_type_all_clusters{clus_id} = replay_clus_char;

%=== visualizing 
figure('visible',options.figure_show);
replay_clus_uni = unique(replay_clus_char);
for i = 1:numel(replay_clus_uni)
    current_clus = find(strcmp(replay_clus_char,replay_clus_uni{i}));
    hold on;
    scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),10,'DisplayName',[replay_clus_uni{i}]);

end
legend;
xlabel ("x"); ylabel ("y"); zlabel ("z");
title('replays colored by position');
hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off location');
hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing location');

%% determining type of flight 
for hide = 1
curr_flight_type = [];
max_dist_same_clus = distance_tolerance;
if size(r_fd,1)>1
% takeoff 
    if (vecnorm(mean_takeoff-r_fd(1,:), 2,2)<max_dist_same_clus | vecnorm(mean_takeoff-r_fd(2,:), 2,2)< max_dist_same_clus)
        curr_flight_type = ['feeder-'];
    else
        curr_flight_type = ['nonfeeder-'];
    end
% landing 
    if (vecnorm(mean_landing-r_fd(1,:), 2,2)<max_dist_same_clus | vecnorm(mean_landing-r_fd(2,:), 2,2)< max_dist_same_clus)
        curr_flight_type = [curr_flight_type, 'feeder'];
    else
        curr_flight_type = [curr_flight_type, 'nonfeeder'];
    end
else
    if (vecnorm(mean_takeoff-r_fd, 2,2)<max_dist_same_clus)
        curr_flight_type = ['feeder-'];
    else
        curr_flight_type = ['nonfeeder-'];
    end
% landing 
    if (vecnorm(mean_landing-r_fd, 2,2)<max_dist_same_clus)
        curr_flight_type = [curr_flight_type, 'feeder'];
    else
        curr_flight_type = [curr_flight_type, 'nonfeeder'];
    end
end
flight_type_all_clusters{clus_id} = curr_flight_type;

end
%% figure 3a visualize
for hide = 1
figure('visible',options.figure_show);
xlabel ("x"); ylabel ("y"); zlabel ("z");
title(curr_flight_type);
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing');

id = find(f_clus.id==clus_id);
% plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001,'HandleVisibility', 'off');
% xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
xlabel('x');    ylabel('y');    hold on;
avg_take_off = [];


for ii=1:size(id,2)
%     title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
    plot(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),'-','LineWidth',1,'Color', 'k','HandleVisibility', 'off');
    avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
end
take_off = mean(avg_take_off,2);    
xlim([-6,6])
xlim([-4,4])
axis equal;
textscatter(take_off(1),take_off(2),"Take-off",'HandleVisibility', 'off');     
hold off;   
figure('visible',options.figure_show);
unique_clusters = unique(replay_clus_char);
for i = 1:numel(unique_clusters)
    current_clus =  find(strcmpi(replay_clus_char, {unique_clusters{i}}));
    hold on;scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'DisplayName',[unique_clusters{i}]);
end

% legend;
xlabel ("x"); ylabel ("y"); zlabel ("z");
% title('replays colored by position');
% hold on; scatter3 (mean_takeoff(1),mean_takeoff(2),mean_takeoff(3),'filled','DisplayName','take off');
% hold on; scatter3 (mean_landing(1),mean_landing(2),mean_landing(3),'filled','DisplayName','landing');

id = find(f_clus.id==clus_id);
% plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001,'HandleVisibility', 'off');
% xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
xlabel('x');    ylabel('y');    hold on;
avg_take_off = [];


for ii=1:size(id,2)
%     title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
    plot(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),'-','LineWidth',1,'Color', [0.85 0.85 0.85],'HandleVisibility', 'off');
    avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
end
take_off = mean(avg_take_off,2);    

% textscatter(take_off(1),take_off(2),"Take-off",'HandleVisibility', 'off');     
hold off;   axis equal;
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

replay_type_all_clusters{clus_id} = replay_clus_char;
end
%% calculating occupancy of clusters 
for hide = 1
clusters_occ = [];% the amount of seconds spent in the vicinity of replay location
clusters_all_pos = {};
clusters_replay_num = [];
for i = 1:numel(unique_clusters)
    current_clus =  find(strcmpi(replay_clus_char, {unique_clusters{i}}));
    pos_idx=[];
    for ii = 1: numel(current_clus) % mark all the positions near the replay for a specific cluster
        pos_idx = union (pos_idx, find(vecnorm(r-pos_during_replay(current_clus(ii),:), 2,2)<0.3));
    end
	clusters_all_pos{i} = pos_idx;
    clusters_occ([i]) = numel(pos_idx) /Fs;
    clusters_replay_num([i]) = numel(current_clus);
end

%     hold on;scatter3(pos_during_replay(current_clus,1),pos_during_replay(current_clus,2),pos_during_replay(current_clus,3),5,'DisplayName',[unique_clusters{i}]);

% === sanity check that the algorithm is selecting the right points in r that correspond to the replay clusters  (slow)
% for i = 1:numel(unique_clusters)
%     figure;
%     scatter3(r(:,1),r(:,2),r(:,3), "filled","k"); alpha(0.2); 
%     hold on; scatter3(r(clusters_all_pos{i},1),r(clusters_all_pos{i},2),r(clusters_all_pos{i},3), "filled","r"); alpha(0.2); title (unique_clusters{i})
% end
end
%% histogram of types of replay
for hide = 1
replay_rate = (clusters_replay_num ./clusters_occ);
% replay_rate = replay_rate./sum(replay_rate);
unique_clusters;
figure('visible',options.figure_show); subplot(1,2,1); bar(reordercats(categorical(unique_clusters),unique_clusters), clusters_replay_num); ylabel("Count");
subplot(1,2,2); bar(reordercats(categorical(unique_clusters),unique_clusters), replay_rate);ylabel("#Replay /s");
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== saving some variables 
replay_occ = [];
replay_fraction = [];
for i = 1:numel(unique_clusters)
    current_clus =  find(strcmpi(replay_clus_char, {unique_clusters{i}}));
    replay_occ([current_clus]) = clusters_occ(i);
    replay_fraction([current_clus]) = replay_rate(i);
end 
replay_position_occupancy_all_clusters{clus_id} = replay_occ;
replay_rate_all_clusters{clus_id}= replay_fraction;
end
%% saving the position of replay 
replay_position_all_clusters{clus_id}=pos_during_replay;
wc_p_1_all_clusters {clus_id} = wc_p_1;
wc_p_2_all_clusters {clus_id} = wc_p_2;
replay_segment_len_all_clusters{clus_id} = segment_len_frac;
end
end
%% === session analysis (pool together all clusters) 
%% plotting what happens at replay times when using different place maps
for hide = 1
% current_replay_clus = 2;
% current_replays = zeros(size(subsequent_events));
% for i = 1: size(replay_times_all_clusters{current_replay_clus},1)
%     current_times = replay_times_all_clusters{current_replay_clus}(i,:);
%     current_replays(find(predicted_position_t>current_times(1) & predicted_position_t < current_times(2)))=1;
% end
% figure;
% colormap(hot)
% tiledlayout(4,1)
% ax1 = nexttile;
% area(predicted_position_t,current_replays,0,'FaceColor',[0 0 0],'FaceAlpha',0.5,'LineStyle','none'); 
% ax2 = nexttile;
% imagesc([t_NP_min, t_NP_max],[0,max(Position_Data(:,2))],posteriors_all_clusters{1}',prctile(posteriors_all_clusters{2},[5 99.5],"all")');xlabel("Time"); ylabel("Position bins"); title("posteior cluster 1");
% axis xy;
% ax3 = nexttile;
% imagesc([t_NP_min, t_NP_max],[0,max(Position_Data(:,2))],posteriors_all_clusters{2}',prctile(posteriors_all_clusters{3},[5 99.5],"all")');xlabel("Time"); ylabel("Position bins"); title("posteior cluster 2");
% axis xy;
% ax4 = nexttile;
% imagesc([t_NP_min, t_NP_max],[0,max(Position_Data(:,2))],posteriors_all_clusters{3}',prctile(posteriors_all_clusters{3},[5 99.5],"all")');xlabel("Time"); ylabel("Position bins"); title("posteior cluster 3");
% axis xy;
% linkaxes([ax1,ax2, ax3, ax4],'x')
end
%% plotting replay times at different place maps 
for hide = []
num_col = numel(imp_bat_clusters)-1;
plot_index =1;
figure('units','normalized','outerposition',[0 .05 0.05*num_col 0.95],'visible',options.figure_show);
replay_duration 
% for i = 1:size(replay_times_all_clusters{current_replay_clus},1)
for current_replay_clus = setdiff(imp_bat_clusters,1) 
    padding = 0.01;
% for current_replay_clus = 5
    replay_index_all = find(weighted_corr_all_clusters{current_replay_clus}>0.6 &  replay_score_all_clusters {current_replay_clus} >0.3);
    if numel(replay_index_all > 10)
        replay_index_all = replay_index_all(randperm(numel(replay_index_all),10));
    end
    for replay_index = replay_index_all
        if replay_index >  size(replay_times_all_clusters{current_replay_clus},1)
            break;
        end
        
        
        if plot_index > 10*num_col
            sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
            fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
            figure('units','normalized','outerposition',[0 .05 0.05*num_col 0.95],'visible',options.figure_show);
            plot_index =1;
        end

        current_times_vector = replay_times_all_clusters{current_replay_clus}(replay_index,:);
        current_times = find(predicted_position_t >= current_times_vector(1)-padding & predicted_position_t <= current_times_vector(2)+padding);  
        event_time_vector = [1:numel(current_times)];
        for current_field = setdiff(imp_bat_clusters,1)
            subplot(10,num_col,plot_index);
            current_posteriors = posteriors_all_clusters{current_field}(current_times,:);
            uniform_prob = 1/size(current_posteriors,2);
            [position_prob, decoded_position] = max(current_posteriors');
            current_posteriors(find(position_prob' < 3*uniform_prob),:) = NaN;
            event_time_vector = [1:numel(current_times)];
            colormap(hot);
            imagesc([predicted_position_t(current_times)],[0,1],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;

%             imagesc([predicted_position_t(current_times)],[0,1],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
            if plot_index ==1
                ylabel(["replay cluster " num2str(current_replay_clus)]);
                xlabel("Time (s)"); 
                title(["field " num2str(current_field)])
            elseif plot_index <=num_col
                title(["field " num2str(current_field)])
            else 
                set(gca, 'YTick', []);
            end
            xticks(predicted_position_t(current_times(end)));
            xticklabels([num2str(round(1000*(current_times_vector(2)-current_times_vector(1)))) ' ms']);
            axis xy;
            plot_index = plot_index +1;
        end
    end
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% plotting bat trajectory at the point of replay - picking replays 
for hide = []
%=== finding the replays of interest 
current_replay_clus = 1;
current_location = 'remote';
replay_index_all = find(weighted_corr_all_clusters{current_replay_clus}>0.2 &  replay_score_all_clusters {current_replay_clus} >0.3 & contains(replay_type_all_clusters{current_replay_clus},current_location));

padding = 0.01
figure('units','normalized','outerposition',[0.25 .05 0.7 0.95],'visible',options.figure_show);
tiledlayout(7,10);
if numel(replay_index_all) > 70
    replay_index_all = replay_index_all(randperm(numel(replay_index_all),70));
end
for replay_index = replay_index_all
    nexttile;
    current_times_vector = replay_times_all_clusters{current_replay_clus}(replay_index,:);
    current_times = find(predicted_position_t >= current_times_vector(1)-padding & predicted_position_t <= current_times_vector(2)+padding);  
    event_time_vector = [1:numel(current_times)];
    current_posteriors = posteriors_all_clusters{current_replay_clus}(current_times,:);
    uniform_prob = 1/size(current_posteriors,2);
    [position_prob, decoded_position] = max(current_posteriors');
    current_posteriors(find(position_prob' < 3*uniform_prob),:) = NaN;
    event_time_vector = [1:numel(current_times)];
    colormap(hot);
    imagesc([predicted_position_t(current_times)],[0,1],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
    title(replay_type_all_clusters{current_replay_clus}(replay_index));
    title(replay_index);
    xticks(predicted_position_t(current_times(end)));
    xticklabels([num2str(round(1000*(current_times_vector(2)-current_times_vector(1)))) ' ms']);
    axis xy;
end 

%% plotting bat trajectory at the point of replay - plotting
replay_index = 28;
for hide = []
%=== plotting all flights of cluster of interest 
for hide = 1
    figure('units','normalized','outerposition',[.05 .1 .9 .4],'visible',options.figure_show);
clus_id = current_replay_clus;
id = find(f_clus.id==clus_id); % find the index of all the trajectories with a certain cluster id
xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
xlabel('x');    ylabel('y');    hold on;
avg_take_off = [];

for ii=1:size(id,2)
    title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
    plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(clus_id,:));
    avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
end
take_off = mean(avg_take_off,2);    
textscatter(take_off(1),take_off(2),"Take-off");     
hold off;   axis equal;
end
%=== plotting mean flight trajectory of cluster of interest 
figure; subplot(1,2,1);
plot3(mean_path_3d_all_clusters{clus_id}(1,:),mean_path_3d_all_clusters{clus_id}(2,:),mean_path_3d_all_clusters{clus_id}(3,:))
xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
xlabel('x');    ylabel('y');    hold on;
axis equal;

%=== plotting bat position at the time of replay 

current_times_vector = replay_times_all_clusters{current_replay_clus}(replay_index,:);
% current_times = find(t>current_times_vector(1) & t<current_times_vector(2));
% current_pos = mean(r(current_times,:,:));
current_time = find(t>current_times_vector(1), 1, 'first');
current_pos = r(current_time,:,:);

hold on; scatter3(current_pos(1), current_pos(2), current_pos(3), 10)
%=== plotting 
subplot(1,2,2);
current_times_vector = replay_times_all_clusters{current_replay_clus}(replay_index,:);
current_times = find(predicted_position_t >= current_times_vector(1)-padding & predicted_position_t <= current_times_vector(2)+padding);  
event_time_vector = [1:numel(current_times)];
current_posteriors = posteriors_all_clusters{current_replay_clus}(current_times,:);
uniform_prob = 1/size(current_posteriors,2);
[position_prob, decoded_position] = max(current_posteriors');
current_posteriors(find(position_prob' < 3*uniform_prob),:) = NaN;
event_time_vector = [1:numel(current_times)];
colormap(hot);
imagesc([predicted_position_t(current_times)],[0,1],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
title(replay_type_all_clusters{current_replay_clus}(replay_index));
title(replay_index);
xticks(predicted_position_t(current_times(end)));
xticklabels([num2str(round(1000*(current_times_vector(2)-current_times_vector(1)))) ' ms']);
axis xy;

end
end
%% plotting decoding of flights
for hide = 1

% % adjusting time vector for the visualization of test decoding 
[max_prob_all,predicted_position] = max( posteriers_all_flights{1}'); 
predicted_position=rescale(predicted_position); % converting from position bin to position [0 1]
uniform_prob = 1/size(post_test,2);
predicted_position(find(max_prob_all<uniform_prob*3))=0;
t_NP_max = max(spk(:,1));
t_NP_min = min(spk(:,1));
t_NP_spacing_for_position = (t_NP_max-t_NP_min)/numel(predicted_position);
predicted_position_t_test = [t_NP_min:t_NP_spacing_for_position:t_NP_max-t_NP_spacing_for_position];
% predicted_position_t_test = predicted_position_t;

current_bat_clusters = 2:numel(imp_bat_clusters);
num_col = numel(current_bat_clusters);
plot_index =1;
figure('units','normalized','outerposition',[0.25 .05 0.5 0.95],'visible',options.figure_show);
for current_flight_clus = 1:numel(current_bat_clusters)
        flight_idx = find(f_clus.id==current_flight_clus);
        flight_idx = flight_idx (3);
        current_times_vector = [t(f_clus.strt_frame(flight_idx)),t(f_clus.stop_frame(flight_idx))];
        current_times = find(predicted_position_t_test >= current_times_vector(1) & predicted_position_t_test <= current_times_vector(2));  
        event_time_vector = [1:numel(current_times)];
        for current_field = 1:numel(current_bat_clusters)
            subplot(numel(current_bat_clusters),numel(current_bat_clusters),plot_index);
%             current_posteriors = posteriors_all_clusters{current_field}(current_times,:);
            current_posteriors =  posteriers_all_flights{current_field} (current_times,:);
%             colormap(hot);
            imagesc([predicted_position_t_test(current_times)],[0,1],current_posteriors',prctile(current_posteriors',[1 99],"all")'); hold on;
            if plot_index ==1
                xlabel("Time (s)"); 
                title(["field " num2str(current_field)])
            elseif plot_index <=num_col
                title(["field " num2str(current_field)])
            else 
                set(gca, 'YTick', []);
            end
            if current_field ==1
                ylabel(["Flight cluster " num2str(current_flight_clus)]);
                xlabel([num2str(round((current_times_vector(2)-current_times_vector(1)),2)) ' s']);

            end
            set(gca, 'XTick', []);
            axis xy;
            plot_index = plot_index+1;
    end
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% calculating the same cluster vs diff cluster decoding error 
for hide = 1
% looping through each cluster 
same_clus_error = [];
diff_clus_error = [];
session_st = predicted_position_t(find(predicted_position_t>=0,1,"first"));
session_st_idx = find(predicted_position_t>=0,1,"first");
session_en = min([max(NP_imu.t), max(RPL.t)]);
session_en_idx = find(predicted_position_t<session_en,1,"last");
for current_flight_clus = 2:numel(imp_bat_clusters)
    flight_all = find(f_clus.id==current_flight_clus);
    for ii = 1:numel(flight_all) % loop through each flight within cluster
        distance_from_prev = [];
        st_time = f_clus.strt_frame(flight_all(ii)); % start time of trajectory
        en_time = f_clus.stop_frame(flight_all(ii)); % end time of trajectory
        time_stamps = [st_time:en_time];
        time_actual = t(time_stamps);
        positon_normalized = f_clus.lin_tr{flight_all(ii)};
        time_position = vertcat(time_position, [time_actual positon_normalized]);
    end 
    for flight_current = 1:numel(flight_all)
        flight_idx = flight_all (flight_current);
        current_times_vector = [t(f_clus.strt_frame(flight_idx)),t(f_clus.stop_frame(flight_idx))];
        flight_st = current_times_vector (1);
        flight_en = current_times_vector (2);
        if flight_st < session_st | flight_en > session_en
            continue;
        end
        step = 0.05;
        tau = step;
        
        current_times = find(predicted_position_t_test >= current_times_vector(1) & predicted_position_t_test <= current_times_vector(2));  
%         bin_st = [flight_st:step:flight_en - tau];
        bin_st = predicted_position_t_test(current_times);
        event_time_vector = [1:numel(current_times)];
%         diff_cluster_field = find(imp_bat_clusters~= 1 & imp_bat_clusters ~= current_flight_clus);
%         same_cluster_field = find(imp_bat_clusters~= 1 & imp_bat_clusters == current_flight_clus);
        for current_field = find(imp_bat_clusters~= 1)
            current_posteriors = posteriers_all_flights{current_field} (current_times,:);
            [~,max_posteriors] = max(current_posteriors', [],1);
            max_posteriors_scaled = max_posteriors/size(current_posteriors,2);
            actual_positions = [];
            for i = 1: numel(bin_st)
                 time_position_flight = time_position(find(time_position(:,1) >= bin_st(i) & time_position(:,1) <= bin_st(i)+tau),:);
                 actual_positions = [actual_positions mean(time_position_flight(:,2))];
            end
            if current_field == current_flight_clus
                same_clus_error = [same_clus_error sqrt(sum((max_posteriors_scaled-actual_positions).^2)/numel(max_posteriors_scaled))];
            else
                diff_clus_error = [diff_clus_error sqrt(sum((max_posteriors_scaled-actual_positions).^2)/numel(max_posteriors_scaled))];
            end
        end
    end
end

figure;
histogram(same_clus_error,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    hold on;
histogram(diff_clus_error,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
xlabel("RMS error");
ylabel("fraction");
legend("same cluster decoding error flight","different cluster decoding error flight");
end
%% confusion matrix of replay metrics
for hide = 1
weighted_corr_mat = [];
avg_jump_dist_mat = [];
replay_scores_mat = [];
post_spread_mat =[];
warning("off");
for current_replay_clus = imp_bat_clusters
    for current_field = imp_bat_clusters
        current_replay_times = replay_times_all_clusters {current_replay_clus};
        weighted_corr_replay = [];
        avg_jump_distance_replay = [];
        replay_scores_replay = [];
        posterior_spread_replay =[];
        for i = 1:size(current_replay_times, 1)
            current_times_vector = replay_times_all_clusters{current_replay_clus}(i,:);
            current_times = find(predicted_position_t >= current_times_vector(1) & predicted_position_t <= current_times_vector(2));  
            current_posteriors = posteriors_all_clusters{current_field}(current_times,:);
            [weighted_corr_replay(i),~, avg_jump_distance_replay(i), ~, replay_scores_replay(i), ~,~, posterior_spread_replay(i), ~] = evaluate_candidate_event_v3(current_posteriors);
        end
%         weighted_corr_replay(find(isnan(weighted_corr_replay)))=0;
        weighted_corr_mat(current_replay_clus, current_field) = mean(abs(weighted_corr_replay(find(~isnan(weighted_corr_replay)))));
        avg_jump_dist_mat(current_replay_clus, current_field) = mean(abs(avg_jump_distance_replay(find(~isnan(avg_jump_distance_replay)))));
        replay_scores_mat(current_replay_clus, current_field) = mean(abs(replay_scores_replay(find(~isnan(replay_scores_replay)))));
        post_spread_mat(current_replay_clus, current_field) = mean(abs(posterior_spread_replay(find(~isnan(posterior_spread_replay)))));
    end
end
warning("on");
figure('visible',options.figure_show);
subplot(2,2,1);heatmap(weighted_corr_mat); xlabel ("Place field"); ylabel("Replay cluster"); title ("average weighted correlation");
subplot(2,2,2);heatmap(avg_jump_dist_mat); title ("average jump distance")
subplot(2,2,3);heatmap(replay_scores_mat); title ("average replay score");
subplot(2,2,4);heatmap(post_spread_mat);  title ("average posterior spread")
end
%% plotting flight duration vs compression times 
for hide = 1
flight_dur_clus = [];
top_replays_num = 5;
%=== plotting temporal compression of all significant events 
figure('visible',options.figure_show);
for clus_id = setdiff(unique(f_clus.id),1)
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));% getting the duration of flight
%     rand_amt = rand(size(temporal_compression_all_clusters{clus_id}))/100;
    rand_amt = zeros(size(temporal_compression_all_clusters{clus_id}));
    flight_dur_plot = flight_dur_clus + rand_amt;
    hold on; scatter (abs(temporal_compression_all_clusters{clus_id}),flight_dur_plot,'filled','DisplayName',['cluster ' num2str(clus_id)]);xlabel('temporal compression of top replays'); ylabel('flight duration');
end
legend();
alpha(0.2); 
xlim([0 50]);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== plotting mean temporal compression of replay
figure('visible',options.figure_show);
for clus_id = setdiff(unique(f_clus.id),1)
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));
    rand_amt = zeros(top_replays_num);
%     rand_amt = rand(size(temporal_compression_all_clusters{clus_id}))/100;
%     rand_amt = zeros(size(temporal_compression_all_clusters{clus_id}));
    flight_dur_plot = flight_dur_clus ;
    hold on; scatter (mean(abs(temporal_compression_all_clusters{clus_id})), flight_dur_clus,'d','filled','DisplayName',['cluster ' num2str(clus_id)]);xlabel('temporal compression of replay (mean)'); ylabel('flight duration');
end
legend();
legend('Location','southeast');
% xlim([0 50]);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== plotting all replay speeds
figure('visible',options.figure_show);
for clus_id = setdiff(unique(f_clus.id),1)
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));
    flight_len_clus = mean(f_clus.length(f_clus.id==clus_id));
    rand_amt = zeros(top_replays_num);
%     rand_amt = rand(size(temporal_compression_all_clusters{clus_id}))/100;
%     rand_amt = zeros(size(temporal_compression_all_clusters{clus_id}));
    hold on; scatter (mean(abs(replay_velocity_all_clusters{clus_id})), flight_dur_clus,'d','filled','DisplayName',['cluster ' num2str(clus_id)]);xlabel('mean replay speed of all replays (m/s)'); ylabel('flight duration');
end
legend('Location','northwest');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== plotting temporal compression of top replays 
figure('visible',options.figure_show);
for clus_id = setdiff(unique(f_clus.id),1)
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));% getting the duration of flight

    rand_amt = zeros(top_replays_num,1);
    flight_dur_plot = flight_dur_clus + rand_amt;
    hold on; scatter (abs(temporal_compression_all_clusters{clus_id}([1:top_replays_num])),flight_dur_plot,'filled','DisplayName',['cluster ' num2str(clus_id)]);xlabel('temporal compression of top replays'); ylabel('flight duration');
end
legend();
xlim([0 50]);
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

%=== plotting top replay speeds
figure('visible',options.figure_show);
for clus_id = setdiff(unique(f_clus.id),1)
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));
    flight_len_clus = mean(f_clus.length(f_clus.id==clus_id));
    rand_amt = zeros(top_replays_num);
%     rand_amt = rand(size(temporal_compression_all_clusters{clus_id}))/100;
%     rand_amt = zeros(size(temporal_compression_all_clusters{clus_id}));
    hold on; scatter (mean(abs(replay_velocity_all_clusters{clus_id}([1:top_replays_num]))), flight_dur_clus,'d','filled','DisplayName',['cluster ' num2str(clus_id)]);xlabel('mean replay speed of top replays(m/s)'); ylabel('flight duration');
end
legend('Location','northwest');
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);



end
%% plotting when replay is happening 
for hide = 1
flight_st = f_clus.strt_frame;
flight_en = f_clus.stop_frame;
flight_st_t = (t(flight_st))';
flight_en_t = (t(flight_en))';

% plotting start of all flights in grey 
figure('units','normalized','outerposition',[.2 0 .6 1],'visible',options.figure_show);
n_row = n_surv_clusters*3;
tiledlayout(n_row,1,'TileSpacing','none');
col_clus = hsv(n_surv_clusters);

for i = 1: n_surv_clusters
    ax1 = nexttile;
    ids = find(f_clus.id==i);
    curr_flight_st = (t(f_clus.strt_frame(ids)))';
    replay_st_times = replay_times_all_clusters{i}(:,1);
%     replay_st_times = replay_st_times(find (abs(weighted_corr_all_clusters{i})> 0.8));
    stem(flight_st_t, ones(numel(flight_st_t)),".","Color","#d3d3d3"); 
    hold on;
    stem(curr_flight_st, ones(numel(curr_flight_st)),".","Color",col_clus(i,:)); 
    xlim([min(t) max(t)]);
    set(gca, 'XTick', []);
    ylabel(["cluster ", num2str(i), ' flights']);
    ax2 = nexttile;
    stem(replay_st_times, ones(numel(replay_st_times)),".","Color",col_clus(i,:)); 
    xlim([min(t) max(t)]);
    ylabel(["cluster ", num2str(i), ' replay']);
    ax3 = nexttile;
    histogram(replay_st_times,250,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    hold on;
    ylabel(["cluster ", num2str(i), ' replay distribution']);
    linkaxes([ax1 ax2 ax3],'x')
end

xlabel("time (s)");
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% examining the distribution of flights and replays (same plots)
for hide = 1
all_replay_times = []; % all replay times excluding first cluster 
all_flight_times = t(f_clus.strt_frame)';
all_replay_wcorr = [];
all_replay_clus_id = [];
all_flight_ids = [];

for i = 2: n_surv_clusters
    ids = find(f_clus.id==i);
    replay_st_times = replay_times_all_clusters{i}(:,1);
    replay_en_times = replay_times_all_clusters{i}(:,2);
%     curr_flight_st = (t(f_clus.strt_frame(ids)))';
    all_replay_wcorr = [all_replay_wcorr weighted_corr_all_clusters{i}];
%     all_flight_times = [all_flight_times curr_flight_st];
    all_replay_times = [all_replay_times replay_st_times'];
    cluster_id_curr = cell(1,numel(replay_st_times));
    [cluster_id_curr{:}] = deal(i);
    all_replay_clus_id = [all_replay_clus_id cell2mat(cluster_id_curr)];
    
    
    cluster_id_curr = cell(1,numel(curr_flight_st));
    [cluster_id_curr{:}] = deal(i);
    all_flight_ids = [all_flight_ids cell2mat(cluster_id_curr)];
end
%=== plotting things seperately 
% figure;
% tiledlayout(4,1,'TileSpacing','none');
% ax1 = nexttile;
% histogram(all_flight_times,50,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;  
% ylabel("distribution of flights");
% ax2 = nexttile;
% histogram(all_replay_times,100,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    hold on;
% ylabel("distribution of replays");
% ax3 = nexttile;
% histogram(RPL.t,100,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    hold on;
% ylabel("distribution of ripples");
% % ax4 = nexttile;
% % plot(NP_imu.t,a_abs_NP);
% % ylable("Acc.");
% linkaxes([ax1 ax2 ax3 ax4],'x')
% sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
% fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

figure('units','normalized','outerposition',[0 0.3 1 0.3],'visible',options.figure_show);

bin_size = 100;
[flight_t_count, flight_t_edge] = histcounts(all_flight_times,500);
[replay_t_c, replay_t_e] = histcounts(all_replay_times,500);
[RPL_t_c, RPL_t_e] = histcounts(RPL.t,500);


% stem(flight_t_edge(1:end-1)+diff(flight_t_edge/2), normalize(flight_t_count,"range"),".");hold on;
% %== peri/nonperi sanity checks 
% area(predicted_position_t,flight_times_peri);hold on;
% area(predicted_position_t,flight_times_nonperi);hold on;
% 
area(replay_t_e(1:end-1)+diff(replay_t_e/2), normalize(replay_t_c,"range"));hold on;
area(RPL_t_e(1:end-1)+diff(RPL_t_e/2), -normalize(RPL_t_c,"range"));hold on,
plot(NP_imu.t,-normalize(abs(a_abs_NP-1),"range")); hold on;
% area(NP_imu.t,-normalize(abs(a_abs_NP),"range"),"FaceColor",[0 0 0]); hold on;
stem(all_flight_times, ones(numel(all_flight_times)),".","Color","b"); hold on;


ylabel("Rate (a.u.)");
xlabel("Time (s)");
yticks([]); 
% legend("Replays", "Ripples", "Flights")
legend("Replays", "Ripples", "Acc.","Flights")
% legend("Peri","Nonperi","Replays", "Ripples", "Acc.","Flights")
xlim([0,min([max(NP_imu.t), max(RPL_t_e),max(replay_t_e)])])
% xlim([4600,6000]);
% % histogram(normalize(RPL.t,"range"),100,'facealpha',.5,'edgecolor','none','FaceColor','k');  
end
%% replay rate of clusters near flight 
for hide = 1
% replay_times_all_clusters % all replay times 
% % all flight times 
t_total = 20; % time before/after flight 
t_bin = 4; % time bin 
t_absolute_pre = [-t_total+t_bin/2:t_bin:-t_bin/2];
t_absolute_post = [t_bin/2:t_bin:t_total-t_bin/2];

replay_rate_cluster_all_flights = {}; % {flight}(flight cluster of replay,time)
for i=1:numel(f_clus.strt_frame)
    flight_st = t(f_clus.strt_frame(i));
    if flight_st <0
        continue;
    end
    t_pre_flight = [flight_st-t_total:t_bin:flight_st];
    replay_rate_cluster = zeros(numel(imp_bat_clusters), t_total/t_bin);
    for c = 1:numel(imp_bat_clusters)
        for t_i = 1:numel(t_pre_flight)-1
            % count the amount of replay in each time bin
            replay_rate_cluster (c,t_i) = numel(find(replay_times_all_clusters{c}(:,1)>=t_pre_flight(t_i) &  replay_times_all_clusters{c}(:,1)<t_pre_flight(t_i+1)));
        end
    end
%     replay_rate_cluster=replay_rate_cluster./t_bin;
    replay_rate_cluster_all_flights{i} = replay_rate_cluster;
end

replay_rate_cluster_all_flights_post = {}; % {flight}(flight cluster of replay,time)
for i=1:numel(f_clus.stop_frame)
    flight_en = t(f_clus.stop_frame(i));
    if flight_st <0
        continue;
    end
    t_post_flight = [flight_en:t_bin:flight_en + t_total];
    replay_rate_cluster = zeros(numel(imp_bat_clusters), t_total/t_bin);
    for c = 1:numel(imp_bat_clusters)
        for t_i = 1:numel(t_post_flight)-1
            % count the amount of replay in each time bin
            replay_rate_cluster (c,t_i) = numel(find(replay_times_all_clusters{c}(:,1)>=t_post_flight(t_i) &  replay_times_all_clusters{c}(:,1)<t_post_flight(t_i+1)));
        end
    end
%     replay_rate_cluster=replay_rate_cluster./t_bin;
    replay_rate_cluster_all_flights_post{i} = replay_rate_cluster;
end

%=== merge flights of same cluster
replay_rate_merge_flight_type = {}; % {flight cluster} (replay cluster, time bins)
for f = imp_bat_clusters
   replay_rate_temp =  zeros(numel(imp_bat_clusters), t_total/t_bin);
   for i=find(f_clus.id==f)
        if numel(replay_rate_cluster_all_flights{i})==0
            continue;
        end
        replay_rate_temp = replay_rate_temp +  replay_rate_cluster_all_flights{i};
   end
   replay_rate_merge_flight_type{f} = replay_rate_temp;
end

replay_rate_merge_flight_type_post = {}; % {flight cluster} (replay cluster, time bins)
for f = imp_bat_clusters
   replay_rate_temp =  zeros(numel(imp_bat_clusters), t_total/t_bin);
   for i=find(f_clus.id==f)
        if numel(replay_rate_cluster_all_flights_post{i})==0
            continue;
        end
        replay_rate_temp = replay_rate_temp +  replay_rate_cluster_all_flights_post{i};
   end
   replay_rate_merge_flight_type_post{f} = replay_rate_temp;
end

%=== visualization
for i=2:numel(imp_bat_clusters)% flight clusters
    figure('visible',options.figure_show);
    tiledlayout(1,2);
    for c = 2:numel(imp_bat_clusters)%replay clusters 
        subplot(1,2,1);
        hold on; plot(t_absolute_pre, replay_rate_merge_flight_type{i}(c,:),"Color",col_clus(c,:));
        hold on; stem(0,0, "Color", col_clus(i,:))
        subplot(1,2,2);
        hold on; plot(t_absolute_post, replay_rate_merge_flight_type_post{i}(c,:),"Color",col_clus(c,:));
        hold on; stem(0,0, "Color", col_clus(i,:))
        sgtitle(['Replay rate before/after cluster ' num2str(i)])
    end
end
end
%% calculate distribution of replay peri- vs far away from flight 
for hide = 1
peri_time = 30; % in units of seconds 
flight_times_peri = zeros(1,numel(predicted_position_t));
flight_times_nonperi = zeros(1,numel(predicted_position_t));

% setting peri times to 1 in a vector that has every time point
for i=f_clus.strt_frame
    flight_t = t(i);
    flight_t_pre = t(i)-peri_time;
    current_times = find(predicted_position_t >= flight_t_pre & predicted_position_t <= flight_t);  
    flight_times_peri(current_times) = 1;
end
for i=f_clus.stop_frame
    flight_t = t(i);
    flight_t_post = t(i)+peri_time;
    current_times = find(predicted_position_t >= flight_t & predicted_position_t <= flight_t_post);  
    flight_times_peri(current_times) = 1;
end
% excluding flight times 
for ii=1:numel(f_clus.strt_frame)
    flight_st = t(f_clus.strt_frame(ii));
    flight_en = t(f_clus.stop_frame(ii));
    current_times = find(predicted_position_t >= flight_st & predicted_position_t <= flight_en);  
    flight_times_peri(current_times) = 0;
end
peri_st = strfind(flight_times_peri, [0 1]);
peri_en = strfind(flight_times_peri, [1 0]);
if flight_times_peri(1)==1
    peri_st = [1 peri_st];
end
if flight_times_peri(end)==1
    peri_en = [peri_en numel(flight_times_peri)];
end
% acocunting for start and end 
keep_idx = ones(1,numel(peri_st));
session_st = predicted_position_t(find(predicted_position_t>=0,1,"first"));
session_st_idx = find(predicted_position_t>=0,1,"first");
session_en = min([max(NP_imu.t), max(RPL_t_e)]);
session_en_idx = find(predicted_position_t<session_en,1,"last");
% session_en = predicted_position_t(find(predicted_position_t<session_end_t,1,"last"));
for i=1:numel(peri_st)
    if predicted_position_t(peri_st (i)) <session_st & predicted_position_t(peri_en(i))<session_st
        keep_idx(i) = 0;
    elseif predicted_position_t(peri_st(i)) <session_st
        peri_st(i) = session_st_idx;
    end
    if predicted_position_t(peri_st (i)) >session_en & predicted_position_t(peri_en(i))>session_en
        keep_idx(i) = 0;
    elseif predicted_position_t(peri_en(i)) >session_en
        peri_en(i) = session_en_idx;
    end
end
peri_st = peri_st(find(keep_idx));
peri_en = peri_en(find(keep_idx));

% generating all replay times 
replay_times_merged = [];
for i = 1:numel(imp_bat_clusters)
    replay_times_merged = [replay_times_merged replay_times_all_clusters{i}(:,1)']; % getting all start times of replays of a cluster
%     replay_times_all_clusters{i}(:,2) % alternatively, getting all end times 
end
accounted_replays = zeros(1,numel(replay_times_merged));
% counting peri replays 
peri_replay_count = 0;
peri_dur = 0;
pnp_label_all_RP = string(num2str(nan(1,numel(replay_times_merged))));
for i=1:numel(peri_st)
    peri_st_t = predicted_position_t(peri_st(i));
    peri_en_t = predicted_position_t(peri_en(i));
    peri_replay_count = peri_replay_count + sum(replay_times_merged >= peri_st_t & replay_times_merged <= peri_en_t);  
    pnp_label_all_RP (find(replay_times_merged >= peri_st_t & replay_times_merged <= peri_en_t)) = "P";
    peri_dur = peri_dur + peri_en_t-peri_st_t;
    accounted_replays(find(replay_times_merged >= peri_st_t & replay_times_merged <= peri_en_t))=1;
end
% counting non-peri replays 
non_peri_st =[session_st_idx peri_en];
non_peri_en =[peri_st session_en_idx];

for i=1:numel(non_peri_st)
    flight_times_nonperi(non_peri_st(i):non_peri_en(i)) = 1;
end

nonperi_replay_count = 0;
nonperi_dur = 0;
for i=1:numel(non_peri_st)
    non_peri_st_t = predicted_position_t(non_peri_st(i));
    non_peri_en_t = predicted_position_t(non_peri_en(i));
    nonperi_replay_count = nonperi_replay_count + sum(replay_times_merged >= non_peri_st_t & replay_times_merged <= non_peri_en_t);  
    pnp_label_all_RP (find(replay_times_merged >= non_peri_st_t & replay_times_merged <= non_peri_en_t)) = "NP";
    nonperi_dur = nonperi_dur + non_peri_en_t-non_peri_st_t;
    accounted_replays(find(replay_times_merged >= non_peri_st_t & replay_times_merged <= non_peri_en_t))=1;
end
% sanity checks 
% peri_replay_count + nonperi_replay_count+sum(replay_times_merged>=session_en)+sum(replay_times_merged<=session_st)
% size(replay_times_merged)
% replay_times_merged(find(accounted_replays==0))

peri_replay_rate = peri_replay_count/peri_dur;
nonperi_replay_rate = nonperi_replay_count/nonperi_dur;

end
%% peri vs non-peri for different clusters 
for hide = 1
pnp_replay_count_by_cluster = zeros(numel(imp_bat_clusters),2); % (cluster, peri vs nonperi)
pnp_replay_rate_cluster =  zeros(numel(imp_bat_clusters),2);
% counting peri replays 
peri_replay_count = 0;
for c = 1:numel(imp_bat_clusters)
    for i=1:numel(peri_st)
        peri_st_t = predicted_position_t(peri_st(i));
        peri_en_t = predicted_position_t(peri_en(i));
        pnp_replay_count_by_cluster(c,1) = pnp_replay_count_by_cluster(c,1) + sum(replay_times_all_clusters{c}(:,1)>=peri_st_t & replay_times_all_clusters{c}(:,2)<peri_en_t & replay_score_all_clusters{c}'>0.4);
    end
    for i=1:numel(non_peri_st)
        non_peri_st_t = predicted_position_t(non_peri_st(i));
        non_peri_en_t = predicted_position_t(non_peri_en(i));
        pnp_replay_count_by_cluster(c,2) = pnp_replay_count_by_cluster(c,2) + sum(replay_times_all_clusters{c}(:,1)>=non_peri_st_t & replay_times_all_clusters{c}(:,2)<non_peri_en_t & replay_score_all_clusters{c}'>0.4);
    end
end

pnp_replay_rate_cluster(:,1) = pnp_replay_count_by_cluster(:,1)/peri_dur;
pnp_replay_rate_cluster(:,2) = pnp_replay_count_by_cluster(:,2)/nonperi_dur;

PlotDistr_AF_v1(pnp_replay_count_by_cluster',turbo(size(pnp_replay_count_by_cluster,2)),"Replay count (#replays)",["Peri" "Nonperi"]);
figure('visible',options.figure_show);
PlotDistr_AF_v1(pnp_replay_rate_cluster',turbo(size(pnp_replay_rate_cluster,2)),"Replay rate (#replays/second)",["Peri" "Nonperi"]);
end
%% === visualizing flights by total replay rate
for hide = 1
total_replay_rate = sum(pnp_replay_rate_cluster,2);
% total_replay_rate = total_replay_rate*60; % converts from replay/sec to replay/min
[~,sort_by_rate] = sort(total_replay_rate,"descend");
sort_by_rate = sort_by_rate';
count_subplot = 1;

figure('units','normalized','outerposition',[0 0.3 1 0.3],'visible',options.figure_show);
for clus_id = sort_by_rate
    if clus_id ==1
        continue
    end
    id = find(f_clus.id==clus_id); % find the index of all the trajectories with a certain cluster id
    % plotting trajectory of all flights
    subplot(1,numel(clusters_of_interest-1),count_subplot);
    count_subplot = count_subplot+1;
    plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8 0.2],'MarkerSize',0.001);
    xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
    xlabel('x');    ylabel('y');    hold on;
    avg_take_off = [];

    for ii=1:size(id,2)
%         title(['Cluster' num2str(clus_id) ' (' num2str(size(id,2)) ' flights),'])
        title([num2str(total_replay_rate(clus_id)*60,2) ' replay/min']);
        plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(clus_id,:));
        avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
    end
    take_off = mean(avg_take_off,2);    
    textscatter(take_off(1),take_off(2),"Take-off");     
    hold off;   axis equal;
end
sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% plotting sequences of replays
for hide =[]
%=== using just replays 
for hide = []
% sorting top_idx by selected variable
sort_by_variable = all_replay_times;
% sort_by_variable = slope_replay;
[~, sorted_index]=sort(sort_by_variable, 'ascend');
all_replay_times_sorted = all_replay_times(sorted_index);
all_replay_wcorr_sorted = all_replay_wcorr(sorted_index);
all_replay_clus_id_sorted = all_replay_clus_id(sorted_index);
x_sequence = 1: numel(sorted_index);

% basic plot 
figure('units','normalized','outerposition',[0 0.3 1 0.3],'visible',options.figure_show);

for i = 2: n_surv_clusters
    ids = find (all_replay_clus_id_sorted==i);
    for curr_id = ids
        hold on; rectangle ("Position",[curr_id, 0,1,1],"FaceColor",col_clus(i,:),'EdgeColor',col_clus(i,:),'LineWidth',0.01);
    end
end
xlim([0,numel(all_replay_clus_id_sorted)])
% forward and rev replay 
figure('units','normalized','outerposition',[0 0.3 1 0.3],'visible',options.figure_show);

for i = 2: n_surv_clusters
    ids_f = find (all_replay_clus_id_sorted==i & all_replay_wcorr_sorted > 0); 
    for curr_id = ids_f
        hold on; rectangle ("Position",[curr_id, 0,1,1],"FaceColor",[col_clus(i,:)],'LineStyle', "none");
    end
    ids_r = find (all_replay_clus_id_sorted==i & all_replay_wcorr_sorted < 0); 
    for curr_id = ids_r
        hold on; rectangle ("Position",[curr_id, 0,1,1],"FaceColor",[col_clus(i,:) 0.4],'LineStyle', "none");
    end
end
xlim([0,numel(all_replay_clus_id_sorted)]);
end
%=== using replays AND flights 
% sorting top_idx by selected variable
all_replay_flight_times = [all_replay_times all_flight_times];
all_replay_flight_wcorr = [all_replay_wcorr NaN(size(all_flight_times))];
all_replay_flight_id = [all_replay_clus_id,all_flight_ids];

sort_by_variable = all_replay_flight_times;
% sort_by_variable = slope_replay;
[~, sorted_index]=sort(sort_by_variable, 'ascend');
all_replay_flight_times_sorted = all_replay_flight_times(sorted_index);
all_replay_flight_wcorr_sorted = all_replay_flight_wcorr(sorted_index);
all_replay_flight_id_sorted = all_replay_flight_id(sorted_index);

% forward vs rev replay and flight vs replay
figure('units','normalized','outerposition',[0 0.3 1 0.3],'visible',options.figure_show);

for i = 2: n_surv_clusters
    ids_f = find (all_replay_flight_id_sorted==i & all_replay_flight_wcorr_sorted > 0); 
    for curr_id = ids_f
        hold on; rectangle ("Position",[curr_id, 0,1,1],"FaceColor",[col_clus(i,:)],'LineStyle', "none");
    end
    ids_r = find (all_replay_flight_id_sorted==i & all_replay_flight_wcorr_sorted < 0); 
    for curr_id = ids_r
        hold on; rectangle ("Position",[curr_id, 0,1,1],"FaceColor",[col_clus(i,:) 0.4],'LineStyle', "none");
    end
    ids_flight = find (all_replay_flight_id_sorted==i & isnan(all_replay_flight_wcorr_sorted)); 
    for curr_id = ids_flight
        hold on; rectangle ("Position",[curr_id, 0,1,1.2],"FaceColor",[col_clus(i,:)],'LineStyle', "none");
    end
end
xlim([0,numel(all_replay_flight_id_sorted)]);
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% cross correlation/auto-correlation of replay times & SWR
for hide = []
% [a,b]=cross_correlogram_AF_v1(replay_times_all_clusters, replay_st_times,1000,10);
max_lag = 0.5;
bin_size = 0.025;
% max_lag = 20;
% bin_size = 1;
replay_dur = replay_en_times-replay_st_times;
%=== Plot cross correlation of replay and SWR
figure('units','normalized','outerposition',[0.25 .5 .5 .4],'visible',options.figure_show);
tiledlayout(1,3,'TileSpacing','tight');
[RR_cross_corr,RR_bin_centers] = cross_correlogram_AF_v1(replay_st_times,replay_st_times,max_lag,bin_size);
[RS_cross_corr,RS_bin_centers] = cross_correlogram_AF_v1(replay_st_times+replay_dur/2,RPL.t,max_lag,bin_size);
[SS_cross_corr,SS_bin_centers] = cross_correlogram_AF_v1(RPL.t,RPL.t,max_lag,bin_size);
nexttile; bar(RR_bin_centers,RR_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); xlabel('Time(s)');  ylabel('Fraction'); title('Replay-Replay CC');
nexttile; bar(RS_bin_centers,RS_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); xlabel('Time(s)');  ylabel('Fraction'); title('Replay-SWR CC');
nexttile; bar(SS_bin_centers,SS_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); xlabel('Time(s)');  ylabel('Fraction'); title('SWR-SWR CC');

sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

end
%% cross correlation of replay clusters
for hide = []
% max_lag = 3;
% bin_size = max_lag/20;

max_lag = 10;
bin_size = max_lag/30;
% max_lag = 200;
% bin_size = 1;
% selected_clusters = 1:numel(clusters_of_interest);
selected_clusters = 2:numel(clusters_of_interest);
figure('units','normalized','outerposition',[0.2 0.05 0.7 0.95],'visible',options.figure_show);
tiledlayout(numel(selected_clusters),numel(selected_clusters),'TileSpacing','tight');
%=== all replays 
for x = selected_clusters
    for y=selected_clusters
        nexttile;
        [RR_cross_corr,RR_bin_centers] = cross_correlogram_AF_v1(replay_times_all_clusters{x}(:,1),replay_times_all_clusters{y}(:,1),max_lag,bin_size);
        if x==y
                bar(RR_bin_centers,RR_cross_corr,'FaceColor',col_clus(x,:),'EdgeColor','none','FaceAlpha',1); 
        else
                bar(RR_bin_centers,RR_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); 
        end
        title(['cluster' num2str(x) ' - cluster' num2str(y)]);
    end
end
sgtitle(['All clusters CC ',unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

% === all forward replays 
figure('units','normalized','outerposition',[0.2 0.05 0.7 0.95],'visible',options.figure_show);
tiledlayout(numel(selected_clusters),numel(selected_clusters),'TileSpacing','tight');

for x = selected_clusters
    for y=selected_clusters
        nexttile;
        x_fwd_idx = find(weighted_corr_all_clusters{x}>0);
        y_fwd_idx = find(weighted_corr_all_clusters{y}>0);
        [RR_cross_corr,RR_bin_centers] = cross_correlogram_AF_v1(replay_times_all_clusters{x}(x_fwd_idx,1),replay_times_all_clusters{y}(y_fwd_idx,1),max_lag,bin_size);
        if x==y
                bar(RR_bin_centers,RR_cross_corr,'FaceColor',col_clus(x,:),'EdgeColor','none','FaceAlpha',1); 
        else
                bar(RR_bin_centers,RR_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); 
        end
        title(['cluster' num2str(x) ' - cluster' num2str(y)]);
    end
end
sgtitle(['Fwd Replay CC ',unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
% === all reverse replays 
figure('units','normalized','outerposition',[0.2 0.05 0.7 0.95],'visible',options.figure_show);
tiledlayout(numel(selected_clusters),numel(selected_clusters),'TileSpacing','tight');

for x = selected_clusters
    for y=selected_clusters
        nexttile;
        x_fwd_idx = find(weighted_corr_all_clusters{x}<0);
        y_fwd_idx = find(weighted_corr_all_clusters{y}<0);
        [RR_cross_corr,RR_bin_centers] = cross_correlogram_AF_v1(replay_times_all_clusters{x}(x_fwd_idx,1),replay_times_all_clusters{y}(y_fwd_idx,1),max_lag,bin_size);
        if x==y
                bar(RR_bin_centers,RR_cross_corr,'FaceColor',col_clus(x,:),'EdgeColor','none','FaceAlpha',1); 
        else
                bar(RR_bin_centers,RR_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); 
        end
        title(['cluster' num2str(x) ' - cluster' num2str(y)]);
    end
end
sgtitle(['Rev Replay CC ',unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
% === all reverse replays 
figure('units','normalized','outerposition',[0.2 0.05 0.7 0.95],'visible',options.figure_show);
tiledlayout(numel(selected_clusters),numel(selected_clusters),'TileSpacing','tight');
for x = selected_clusters
    for y=selected_clusters
        nexttile;
        x_fwd_idx = find(weighted_corr_all_clusters{x}>0);
        y_fwd_idx = find(weighted_corr_all_clusters{y}<0);
        [RR_cross_corr,RR_bin_centers] = cross_correlogram_AF_v1(replay_times_all_clusters{x}(x_fwd_idx,1),replay_times_all_clusters{y}(y_fwd_idx,1),max_lag,bin_size);
        if x==y
                bar(RR_bin_centers,RR_cross_corr,'FaceColor',col_clus(x,:),'EdgeColor','none','FaceAlpha',1); 
        else
                bar(RR_bin_centers,RR_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); 
        end
        title(['cluster' num2str(x) ' - cluster' num2str(y)]);
    end
end
sgtitle(['Fwd-Rev Replay CC ',unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);
end
%% plotting replay's time distance to prev/next flight of the SAME replay cluster 
for hide = 1 
%finding the distance from the replay st to the closest flight start 
for omit =[]
replay_dist_flight_start= {};
for i = 1: n_surv_clusters
    ids = find(f_clus.id==i);
    replay_st_times = replay_times_all_clusters{i}(:,1);
    dist_curr_clus = [];
    curr_flight_st = (t(f_clus.strt_frame(ids)))';
    for ii = 1:numel(replay_st_times)
        [~, min_idx] = min(abs(replay_st_times(ii)-curr_flight_st));
        dist_curr_clus = [dist_curr_clus replay_st_times(ii)-curr_flight_st(min_idx)]; % negative: replay before flight start.
    end
    replay_dist_flight_start{i}=dist_curr_clus;
end

%finding the distance from the replay st to the closest flight end 
replay_dist_flight_end= {};
for i = 1: n_surv_clusters
    ids = find(f_clus.id==i);
    replay_st_times = replay_times_all_clusters{i}(:,1);
    dist_curr_clus = [];
    curr_flight_en = (t(f_clus.stop_frame(ids)))';
    for ii = 1:numel(replay_st_times)
        [~, min_idx] = min(abs(replay_st_times(ii)-curr_flight_en));
        dist_curr_clus = [dist_curr_clus replay_st_times(ii)-curr_flight_en(min_idx)]; % negative: replay before flight end.
    end
    replay_dist_flight_end{i}=dist_curr_clus;
end

% plotting all times
% for i = 1: n_surv_clusters
%     figure;
%     current_dist = replay_dist_flight_start{i};
%     hold on; histogram(current_dist(find(abs(weighted_corr_all_clusters{i})> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
%     xlabel ("replay time to closest flight"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
% end

%=== plotting reverse vs forward

for i = 1: n_surv_clusters
    figure('visible',options.figure_show);
    current_dist = replay_dist_flight_start{i};
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to closest flight start"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
end

for i = 1: n_surv_clusters
    figure('visible',options.figure_show);
    current_dist = replay_dist_flight_end{i};
    hold on; histogram(current_dist(weighted_corr_all_clusters{i}> 0.8),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to closest flight end"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
end


%=== visualize replay time to closest flight before/after flight (within a few minutes)
% time_range = 1000;
% for i = 1: n_surv_clusters
%     figure;
%     current_dist = replay_dist_flight_start{i};
%     hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8 & abs(current_dist)<time_range)),100,'facealpha',.5,'edgecolor','none','FaceColor','k');
%     hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8& abs(current_dist)<time_range)),100,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
%     xlabel ("replay time to closest flight start"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
%     legend ("forward", "reverse")
%     xlim([-time_range time_range]);
% end
% 
% for i = 1: n_surv_clusters
%     figure;
%     current_dist = replay_dist_flight_end{i};
%     hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8 & abs(current_dist)<time_range)),100,'facealpha',.5,'edgecolor','none','FaceColor','k');
%     hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8& abs(current_dist)<time_range)),100,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
%     xlabel ("replay time to closest flight end"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
%     legend ("forward", "reverse")
%     xlim([-time_range time_range]);
% end
end

%=== finding the distance of replay to the closest previous/next flight
replay_dist_flight_prev= {};
replay_dist_flight_next= {};
for i = 1: n_surv_clusters
    ids = find(f_clus.id==i);
    replay_st_times = replay_times_all_clusters{i}(:,1);
    curr_flight_st = (t(f_clus.strt_frame(ids)))';
    curr_flight_en = (t(f_clus.stop_frame(ids)))';
    time_intervals = [min(min(t),min(replay_st_times))-1 curr_flight_st max(max(t),max(replay_st_times))+1];
    
    dist_prev_curr = [];
    dist_next_curr = [];
    for ii = 1:(numel(time_intervals)-1)
        current_replay_times = replay_st_times(find(replay_st_times > time_intervals (ii) & replay_st_times < time_intervals (ii+1)))';
        if isempty(current_replay_times)
           continue;
        end
        if ii == 1 % for the very first interval, which is before the first flight
            dist_prev_curr = [dist_prev_curr NaN(1,numel(current_replay_times))];
            dist_next_curr = [dist_next_curr current_replay_times-time_intervals(ii+1)];
        elseif ii == (numel(time_intervals)-1) % last interval, after last flight 
            dist_prev_curr = [dist_prev_curr current_replay_times-curr_flight_en(ii-1)];
            dist_next_curr = [dist_next_curr NaN(1,numel(current_replay_times))];
        else
            dist_prev_curr = [dist_prev_curr current_replay_times-curr_flight_en(ii-1)];
            dist_next_curr = [dist_next_curr current_replay_times-time_intervals(ii+1)];
        end
    
    end
    replay_dist_flight_prev{i} = dist_prev_curr;
    replay_dist_flight_next{i} = dist_next_curr; 
end


for i = 1: n_surv_clusters
    figure('visible',options.figure_show);
    subplot(1,2,2);
    current_dist = replay_dist_flight_prev{i};
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to previous flight (same cluster)"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
    subplot(1,2,1);
    current_dist = replay_dist_flight_next{i};
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to next flight (same cluster)"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
end


sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

end
%% plotting replay's time distance to prev/next flight, ANY flight
for hide = 1 
%=== finding the distance of replay to the closest previous/next flight
replay_dist_flight_prev_any= {};
replay_dist_flight_next_any= {};
for i = 1: n_surv_clusters
    replay_st_times = replay_times_all_clusters{i}(:,1);
    curr_flight_st = (t(f_clus.strt_frame))';
    curr_flight_en = (t(f_clus.stop_frame))';
    time_intervals = [min(min(t),min(replay_st_times))-1 curr_flight_st max(max(t),max(replay_st_times))+1];
    
    dist_prev_curr = [];
    dist_next_curr = [];
    for ii = 1:(numel(time_intervals)-1)
        current_replay_times = replay_st_times(find(replay_st_times > time_intervals (ii) & replay_st_times < time_intervals (ii+1)))';
        if isempty(current_replay_times)
           continue;
        end
        if ii == 1 % for the very first interval, which is before the first flight
            dist_prev_curr = [dist_prev_curr NaN(1,numel(current_replay_times))];
            dist_next_curr = [dist_next_curr current_replay_times-time_intervals(ii+1)];
        elseif ii == (numel(time_intervals)-1) % last interval, after last flight 
            dist_prev_curr = [dist_prev_curr current_replay_times-curr_flight_en(ii-1)];
            dist_next_curr = [dist_next_curr NaN(1,numel(current_replay_times))];
        else
            dist_prev_curr = [dist_prev_curr current_replay_times-curr_flight_en(ii-1)];
            dist_next_curr = [dist_next_curr current_replay_times-time_intervals(ii+1)];
        end
    
    end
    replay_dist_flight_prev_any{i} = dist_prev_curr;
    replay_dist_flight_next_any{i} = dist_next_curr; 
end


for i = 1: n_surv_clusters
    figure('visible',options.figure_show);
    subplot(1,2,2);
    current_dist = replay_dist_flight_prev_any{i};
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to any previous flight"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
    subplot(1,2,1);
    current_dist = replay_dist_flight_next_any{i};
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}> 0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','k');
    hold on; histogram(current_dist(find(weighted_corr_all_clusters{i}< -0.8)),30,'facealpha',.5,'edgecolor','none','FaceColor','r'); 
    xlabel ("replay time to any next flight"); ylabel("count");title(['cluster ' num2str(i) ' replays']);
    legend ("forward", "reverse")
end



sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}, ', cluster ', num2str(clus_id)],'Interpreter','none');
fig_count = saveFig(analysis_directory,[cell2mat(unique_ID) '_flight_' num2str(clus_id)],fig_count,options.savefigures);

end
%% saving detected replays
for hide = 1
num_replays = 0;
cluster_id = [];
posts = [];
comps = [];
vel = [];
corrs = [];

flight_durations = [];
flight_lens = [];
num_replays = 0;

avg_j = [];
max_j = [];
scores = [];
slopes = [];


flight_type=[];
flight_active_cell = [];
replay_type = [];
replay_pos_occ = [];
replay_rate_by_type = [];
replay_st = [];
replay_en = [];
replay_flight_prev_dist = [];
replay_flight_next_dist = [];
replay_flight_prev_dist_any = [];
replay_flight_next_dist_any = [];
flight_errors = [];
peri_rate = [];
nonperi_rate = [];
peri_duration = [];
nonperi_duration=[];
replay_position = [];
%=== new variables 
peri_rate_cluster = [];
nonperi_rate_cluster = [];
session_dur = [];
total_rate_cluster = [];

%=== new variables 
mean_path3D =[];
mean_path1D=[];
mean_time1D=[];

wc_p_1_all = [];
wc_p_2_all = [];
segment_frac_all = [];

for clus_id = clusters_of_interest
    curr_replay_num = numel(avg_jump_all_clusters{clus_id});
    num_replays = num_replays + curr_replay_num;
    
    wc_p_1_all = [wc_p_1_all wc_p_1_all_clusters{clus_id}];
    wc_p_2_all = [wc_p_2_all wc_p_2_all_clusters{clus_id}];
    segment_frac_all = [segment_frac_all replay_segment_len_all_clusters{clus_id}];
    posts = [posts replay_post_all_clusters{clus_id}];
    comps = [comps temporal_compression_all_clusters{clus_id}];
    vel = [vel replay_velocity_all_clusters{clus_id}];
    corrs = [corrs weighted_corr_all_clusters{clus_id}];
    replay_type = [replay_type replay_type_all_clusters{clus_id}];
    replay_pos_occ = [replay_pos_occ replay_position_occupancy_all_clusters{clus_id}];
    replay_rate_by_type = [replay_rate_by_type replay_rate_all_clusters{clus_id}];
%     replay_flight_st_dist = [replay_flight_st_dist replay_dist_flight_start{clus_id}];
%     replay_flight_en_dist = [replay_flight_en_dist replay_dist_flight_end{clus_id}];
    replay_flight_prev_dist = [replay_flight_prev_dist replay_dist_flight_prev{clus_id}];
    replay_flight_next_dist = [replay_flight_next_dist replay_dist_flight_next{clus_id}];
    replay_flight_prev_dist_any = [replay_flight_prev_dist_any replay_dist_flight_prev_any{clus_id}];
    replay_flight_next_dist_any = [replay_flight_next_dist_any replay_dist_flight_next_any{clus_id}];
    replay_position = [replay_position; replay_position_all_clusters{clus_id}];
   
    avg_j = [avg_j avg_jump_all_clusters{clus_id}];
    max_j = [max_j max_jump_all_clusters{clus_id}];
    scores = [scores replay_score_all_clusters{clus_id}];
    slopes = [slopes slope_all_clusters{clus_id}];
    
    % split variable
    replay_time_curr = replay_times_all_clusters{clus_id};
    replay_st = [replay_st replay_time_curr(:,1)'];
    replay_en = [replay_en replay_time_curr(:,2)'];
    
    %=== saving information about cluster
    flight_dur_clus = mean(f_clus.dur(f_clus.id==clus_id));
    flight_dur_curr = cell(1,curr_replay_num);
    [flight_dur_curr{:}] = deal(flight_dur_clus);
    flight_durations = [flight_durations flight_dur_curr];
    
    flight_len_clus = mean(f_clus.length(f_clus.id==clus_id));
    flight_len_curr = cell(1,curr_replay_num);
    [flight_len_curr{:}] = deal(flight_len_clus);
    flight_lens = [flight_lens flight_len_curr];
    
    flight_type_curr = cell(1,curr_replay_num);
    [flight_type_curr{:}] = deal(flight_type_all_clusters{clus_id});
    flight_type = [flight_type flight_type_curr];

    flight_act_curr = cell(1,curr_replay_num);
    [flight_act_curr{:}] = deal(flight_active_cells_all_clusters{clus_id});
    flight_active_cell = [flight_active_cell flight_act_curr];

    flight_err_curr = cell(1,curr_replay_num);
    [flight_err_curr{:}] = deal(flight_error_all_clusters{clus_id});
    flight_errors = [flight_errors flight_err_curr];
    
    cluster_id_curr = cell(1,curr_replay_num);
    [cluster_id_curr{:}] = deal(clus_id);
    cluster_id = [cluster_id cluster_id_curr];
    
    peri_rate_cluster_curr = cell(1,curr_replay_num);
    [peri_rate_cluster_curr{:}] = deal(pnp_replay_rate_cluster(clus_id,1));
    peri_rate_cluster = [peri_rate_cluster peri_rate_cluster_curr];
    
    nonperi_rate_cluster_curr = cell(1,curr_replay_num);
    [nonperi_rate_cluster_curr{:}] = deal(pnp_replay_rate_cluster(clus_id,2));
    nonperi_rate_cluster = [nonperi_rate_cluster nonperi_rate_cluster_curr];
    
    total_rate_cluster_curr = cell(1,curr_replay_num);
    [total_rate_cluster_curr{:}] = deal(total_replay_rate(clus_id));
    total_rate_cluster = [total_rate_cluster total_rate_cluster_curr];
    
    placeholder = cell(1,curr_replay_num);
    [placeholder{:}] = deal(mean_path_3d_all_clusters{clus_id});
    mean_path3D = [mean_path3D placeholder];
    
    placeholder = cell(1,curr_replay_num);
    [placeholder{:}] = deal(mean_len_1d_all_clusters{clus_id});
    mean_path1D = [mean_path1D placeholder];
    
    placeholder = cell(1,curr_replay_num);
    [placeholder{:}] = deal(mean_time_1d{clus_id});
    mean_time1D = [mean_time1D placeholder];
    
end
%=== saving variables on the flight level
bat_ids = cell(1,num_replays);
[bat_ids{:}] = deal(unique_ID{2});

Fs_session = cell(1,num_replays);
[Fs_session{:}] = deal(Fs);

r_session = cell(1,num_replays);
[r_session{:}] = deal(r);

session_ids = cell(1,num_replays);
[session_ids{:}] = deal([unique_ID{1} '_' unique_ID{2} '_' unique_ID{3}]);

peri_rate = cell(1,num_replays);
[peri_rate{:}] = deal(peri_replay_rate);

nonperi_rate = cell(1,num_replays);
[nonperi_rate{:}] = deal(nonperi_replay_rate);

session_dur = cell(1,num_replays);
[session_dur{:}] = deal (peri_dur+nonperi_dur);

peri_duration = cell(1,num_replays);
[peri_duration{:}] = deal (peri_dur);

nonperi_duration = cell(1,num_replays);
[nonperi_duration{:}] = deal (nonperi_dur);


%% checking for duplicates 
session_replay_t = zeros(1,numel(t)); % varibale (0 or 1) spanning whole session to see whether replay is happening at that time point 
idx_keep = [];
% construct vector of when replay is happening 
for i = 1: numel(replay_st)
    if  cell2mat(cluster_id (i))==1
        continue;
    end
    st = replay_st (i);
    en = replay_en (i);
    session_replay_t (find(t>=st & t<=en)) = 1;
end
% check which replays fall into these categories 
session_replay_t = conv2(session_replay_t, [1 1 1], "same");
event_st = (strfind(session_replay_t, [0 1]));
event_en = (strfind(session_replay_t, [1 0]));
if session_replay_t(1)==1
    event_st = [1 event_st];
end
if session_replay_t(end)==1
    event_en = [event_en numel(session_replay_t)];
end

% find duplicate replays 
replays_per_time = {};                      % the index of all replays of a specific chunk 
duplicates = zeros(1,numel(event_st));       % whether a chunk of time contain duplicates 
for i = 1: numel (event_st)
    event_st_t = t(event_st(i));
    event_en_t = t(event_en(i));
    current_replays = find (replay_st> event_st_t & replay_en <event_en_t & cell2mat(cluster_id) ~=1);
    if numel(current_replays)>1
       duplicates (i) = 1; 
    end
    replays_per_time {i} = current_replays;
end
for i = 1: numel (event_st)
    current_replays = replays_per_time {i};
    if numel(current_replays)==1
       idx_keep = [idx_keep current_replays];
    end
    if numel(current_replays)>1
       all_seg = segment_frac_all(current_replays);
       [~, max_i] = max(all_seg);
       idx_keep = [idx_keep current_replays(max_i)];
    end
end

duplicates_replay = zeros(1,numel(replay_st));
duplicates_replay (idx_keep) =1;
%% assigning variables 
detected_replays.session_id = session_ids;
detected_replays.session_dur = session_dur;
detected_replays.flight_duration = flight_durations;
detected_replays.flight_length = flight_lens;    
detected_replays.cluster_id = cluster_id;
detected_replays.post = posts;
detected_replays.bat_id = bat_ids;
detected_replays.comp = comps;
detected_replays.velocity = vel;
detected_replays.corr = corrs;

detected_replays.avg_j = avg_j;
detected_replays.max_j = max_j;
detected_replays.scores = scores;
detected_replays.slopes = slopes;

detected_replays.flight_type=flight_type;
detected_replays.flight_active_cell = flight_active_cell;
detected_replays.replay_type = replay_type;
detected_replays.replay_pos_occ = replay_pos_occ;
detected_replays.replay_rate_by_type = replay_rate_by_type;
detected_replays.replay_st = replay_st;
detected_replays.replay_en = replay_en;
detected_replays.replay_flight_prev_dist = replay_flight_prev_dist;
detected_replays.replay_flight_next_dist = replay_flight_next_dist;
detected_replays.replay_flight_prev_dist_any = replay_flight_prev_dist_any;
detected_replays.replay_flight_next_dist_any = replay_flight_next_dist_any;
detected_replays.flight_errors = flight_errors;
detected_replays.peri_rate=peri_rate;
detected_replays.nonperi_rate=nonperi_rate;

detected_replays.peri_rate_cluster=peri_rate_cluster;
detected_replays.nonperi_rate_cluster=nonperi_rate_cluster;
detected_replays.total_rate_cluster = total_rate_cluster;
detected_replays.pnp_label_all_RP = pnp_label_all_RP;
detected_replays.peri_duration = peri_duration;
detected_replays.nonperi_duration = nonperi_duration;

detected_replays.mean_path3D=mean_path3D;
detected_replays.mean_path1D=mean_path1D;
detected_replays.mean_time1D = mean_time1D;

detected_replays.wc_p_1_all = wc_p_1_all;
detected_replays.wc_p_2_all = wc_p_2_all;
detected_replays.duplicates_replay =  duplicates_replay;
detected_replays.replay_position = replay_position';
detected_replays.Fs_session = Fs_session;
% detected_replays.r_session = r_session;
end
%% saving position and time vectors 
save([analysis_directory,'/Analyzed_Replays_',current_version, unique_ID{1}, '_', unique_ID{2}, '_' unique_ID{3},'.mat'],'detected_replays');
session_info.curr_session_id = [unique_ID{1} '_' unique_ID{2} '_' unique_ID{3}];
session_info.t = {t};
session_info.r = {r};
session_info.f_clus = f_clus;
session_info.same_clus_error = {same_clus_error};
session_info.diff_clus_error = {diff_clus_error};
session_info.peri_dur = peri_dur;
session_info.nonperi_dur = nonperi_dur;
session_info.session_st = session_st;
session_info.session_en = session_en;
save([analysis_directory,'/session_info_',current_version, unique_ID{1}, '_', unique_ID{2}, '_' unique_ID{3},'.mat'],'session_info','-v7.3');

