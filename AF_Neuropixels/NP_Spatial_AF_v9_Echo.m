function NP_Spatial_AF_v9_Echo(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

fig_visibility = 'on';

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    mrgFile = dir('merges_*');              load([mrgFile.folder '/' mrgFile.name]);    % This contains the ids of potential merges              
    nrnFile = dir('SU_kilosort*');                                                      % Single Units
    unique_ID = options.unique_ID;                                                      % Session identifier
    load('LFP_probe2.mat');                                                             % This file contains LFP data from probe 2
    load('RPL_probe2.mat');                                                             % This file contains Ripple data from probe 2
    NP_unit = table();                  MUA_unit = table();                             % Single Units and MUA structures
    n_probes = size(nrnFile,1);                                                         % Number of probes
    merge_cells = 1;
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        
        %=== If merging potential duplicates
        if merge_cells
            if ~isempty(merges{i,1})
                for jj=1:size(merges{i,1},1)
                    toBeMerged = find(any(out.good_units.cluster_id == merges{i,1}{jj,1},2));
                    spikeTimes_usec = vertcat(out.good_units.spikeTimes_usec{toBeMerged,1});
                    temp_unit = out.good_units(end,:);
                    %temp_unit.cluster_id = temp_unit.cluster_id+1;
                    temp_unit.cluster_id = out.good_units.cluster_id(toBeMerged(1));
                    temp_unit.spikeTimes_usec = {spikeTimes_usec};
                    out.good_units = [out.good_units;temp_unit];
                end
                toBeDeleted = find(any(out.good_units.cluster_id == horzcat(merges{i,1}{:}),2));
                out.good_units(toBeDeleted,:) = [];
            end
        end
        
        NP_unit = [NP_unit; out.good_units];                                            % Single units and MUA are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;                                                                          % Clear unnecessary variables
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);             % Add average firing rate
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);           % Add average firing rate
    
    %=== Echolocation Clicks
    if ~isempty(dir(fullfile('Detected_Clicks_*')))
        Clck_file = dir(fullfile('Detected_Clicks_*'));                     % File with times of the detected echolocation clicks
        load([Clck_file.folder,'\',Clck_file.name]);
    else
        Detected_Clicks =[];
        Detected_Clicks.times = [];
    end
    
    %=== Create analysis folder for storing the results
    options.savedata = 1;                                   % If creating folder for saving the data or not
    options.savefigures = 1;                                % Save Figures
    fig_count = 1;                                          % Id of the first figure
    if ~exist('folder_name') && options.savedata
        analysis_directory = [replace(pwd,'NP_datasets','NP_datasets_Analysis_Figures'),'\NP&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
        if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
    elseif exist('folder_name') && options.savedata
        analysis_directory = [replace(pwd,'NP_datasets',folder_name),'\NP&BHv_Analysis_',datestr(now, 'yymmdd_HHMM')];
        if ~exist(analysis_directory,'dir');mkdir(analysis_directory);end
    elseif ~options.savedata
        analysis_directory = 'none';
        options.savefigures = 0;
    end
    
end

%% MAKE SURE BEHAVIOR AND EPHYS HAVE THE SAME DURATION
for hide=1
    
    % Make sure t is a column vector
    if isrow(t),t=t';end
    
    %=== Show difference between recording times
    dt_str = NP_imu.t(1)-t(1);
    dt_end = NP_imu.t(end)-t(T);
    disp(['NP started ',num2str(-dt_str),' s before behavior']);
    disp(['NP ended '  ,num2str(dt_end),' s after behavior']);
    
    %=== Cut the end of recordings accordingly
    if dt_end<0
        %=== Find end of recording
        t_max = NP_imu.t(end);
        T_max = find(diff(t<t_max)==-1);
        
        %=== Cut behavioral variables
        T = T_max;
        bflying = bflying(1:T_max,1);
        r = r(1:T_max,:);
        v = v(1:T_max,:);
        t = t(1:T_max,1);
        v_abs = v_abs(1:T_max,1);
        f_smp = f_smp(all(f_smp<T_max,2),:);
        f_num = size(f_smp,1);
    end
    
end

%% BASIC PARAMS, ASSIGNMENTS
for hide=1
    
    %=== Recluster Flights (default is .8)
    alpha_clus = .8;
    if strcmp(unique_ID{1},'Dataset_2'),alpha_clus = 1.2;end    % More liberal clustering for Field Station Dataset
    clear f_clus;
    f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
    while numel(unique(f_clus.id))==1 && alpha_clus<1.5         % Try increasing the clustering param if only 1 cluster
        alpha_clus = alpha_clus+0.1;
        clear f_clus;   close all;
        f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Params
    if isrow(t),t=t';end                                    % Make sure t is a column vector
    NP_unit_FS = NP_unit(NP_unit.fr>2,:);                   % Store high firing cells (putative interneurons) in a different variable (def 1 Hz)
    NP_unit = NP_unit(NP_unit.fr<2,:);                      % Exclude neurons with high-firing rate
    n_cells = size(NP_unit,1);                              % Number of units
    n_cells_FS = size(NP_unit_FS,1);                        % Number of units (putative interneurons)
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 10;                                             % Number of repetitions for shuffling (Place Cells)
    n_rep_sweeps = 3;                                       % Number of repetitions for shuffling (Theta Sweeps)
    times_to_shift = t(randi([3*Fs T-3*Fs],1,n_rep));       % Time intervals for circshift
    smpls_to_shift = randi([3*Fs T-3*Fs],1,n_rep);          % Sample intervals for circshift
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
    
    spk_dsty_smooth_t = 0.10;                               % Smoothing factor for the spike density
    max_Replay_dur = 1;                                     % Maximal duration for Replay events
    Rate_treshold_SD = 2;                                   % Can try to decrease this one
    Replay_db_time = 0.2;                                   % Debouncing time for Replay detection
    Replay_min_dur = 0.05;                                  % Min Replay duration
    seq_min_spikes = 3;                                     % Minimum number of spikes in the sequence
    seq_min_fr_active = 0.1;                                % Minimum fraction of active cells in the sequence
    
    t_bin_dur = 0.020;                                      % Time bin duration for decoding
    t_bin_ovl = 0.015;                                      % Time bin overlap for decoding
    t_Hi = NP_imu.t;                                        % Time for fast sampling (500 Hz, IMU and LFP)
    bflying_Hi = logical(interp1(t,double(bflying),t_Hi,'nearest',0));
    Fs_Hi = NP_imu.Fs;                                      % Sampling frequency for fast sampling (500 Hz, IMU and LFP)
    Rpl_table = table();                                    % Initialize the Replay table
    SWP_table = table();                                    % Initialize the Sweep table
    
    %=== Populate the s cell with the single units and calculate firing rate
    s = cell(n_cells,n_rep);            % Single units' spikes
    Rate = zeros(length(t),n_cells);    % Smoothed Firing Rate
    for nc = 1:n_cells
        s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
        %=== Clean up duplicated spikes
        too_close = find(diff(s{nc,1})<0.001);
        if ~isempty(too_close)
            %disp([num2str(numel(too_close)), ' cleaned spikes']);
            s{nc,1}(too_close+1) = [];
        end
        Rate(2:end,nc) = histcounts(s{nc,1},t)*Fs;                      % Calculate firing rate at behavioral samples
        Rate(:,nc) = smoothdata(Rate(:,nc),1,'gaussian',round(Fs*0.1)); % Smooth firing rate
        for n = 1:n_rep;    s{nc,n+1} = mod(s{nc,1}+times_to_shift(n),t(T)); end        % Shuffling
    end
    Rate_norm = normalize(Rate,1,'range');
    NP_unit = table2struct(NP_unit);                % Convert NP unit to structure
    [~,rate_sorted_ids] = sort([NP_unit.fr]);       % Get sorting based on firing frequency
    
    %=== Populate the mua cell with the single units and calculate firing rate
    n_mua = size(MUA_unit,1);              % Number of mua units
    s_mua = cell(n_mua,1);                 % Single units' spikes
    for nc = 1:n_mua
        s_mua{nc,1} = MUA_unit.spikeTimes_usec{nc,1}/1e6;                   % Each row contains the time (in s) of the spikes from a single unit
        too_close = find(diff(s_mua{nc,1})<0.001);                          % Clean up duplicated spikes
        if ~isempty(too_close)
            %disp([num2str(numel(too_close)), ' cleaned spikes']);
            s_mua{nc,1}(too_close+1) = [];
        end
    end
    MUA_unit = table2struct(MUA_unit);                    % Convert MUA unit to structure
    
    %=== Populate the s cell for putative interneurons with the single units and calculate firing rate
    s_FS = cell(n_cells_FS,n_rep);               % Single units' spikes
    Rate_FS = zeros(length(t),n_cells_FS);    % Smoothed Firing Rate
    for nc = 1:n_cells_FS
        s_FS{nc,1} = NP_unit_FS.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
        %=== Clean up duplicated spikes
        too_close = find(diff(s_FS{nc,1})<0.001);
        if ~isempty(too_close)
            %disp([num2str(numel(too_close)), ' cleaned spikes']);
            s_FS{nc,1}(too_close+1) = [];
        end
        Rate_FS(2:end,nc) = histcounts(s_FS{nc,1},t)*Fs;                      % Calculate firing rate at behavioral samples
        Rate_FS(:,nc) = smoothdata(Rate_FS(:,nc),1,'gaussian',round(Fs*0.1)); % Smooth firing rate
        for n = 1:n_rep;    s_FS{nc,n+1} = mod(s_FS{nc,1}+times_to_shift(n),t(T)); end        % Shuffling
    end
    Rate_FS_norm = normalize(Rate_FS,1,'range');
    NP_unit_FS = table2struct(NP_unit_FS);                % Convert NP unit to structure
    [~,rate_sorted_ids_FS] = sort([NP_unit_FS.fr]);       % Get sorting based on firing frequency
    
    %=== Transform the LFP into double and interpolate at 500 Hz sampling
    LFP = interp1(red_out.t_ds,double(red_out.lfp).*red_out.voltage_scaling,t_Hi,'linear','extrap');     % Interpolate LFP at the accelerometer time samples
    LFP_mn = normalize(mean(LFP,2),'range');                                                             % Calculate average LFP
    LFP_th = bandpass(LFP_mn,[4 12],Fs_Hi);                                                              % Filtered in the theta band
    LFP_lp = bandpass(LFP_mn,[0.5 20],Fs_Hi);                                                            % Low pass
    LFP_rp = bandpass(RPL_out.LFP_trace,[100 200],RPL_out.Fs);                                           % Ripple trace
    LFP_th_power = smoothdata(zscore(abs(hilbert(LFP_th)),[],1),1,'gaussian',0.05*Fs_Hi);                % Theta power
    LFP_lw = interp1(t_Hi,LFP_mn,t,'linear',median(LFP_mn));                                             % LFP at t sampling
    [PS_LFP,freq_SG] = cwt(LFP_lw,Fs);                                                                   % Spectrogram
    PSD_low = sum(abs(PS_LFP(knnsearch(freq_SG,11):knnsearch(freq_SG,1),:)),1);                          % Spectral power at low frequences
    
    %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                     % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_Hi);              % Filtered at the wing-beat frequency
    [up,lo] = envelope(a_flt_NP,round(0.06*Fs_Hi),'peak');  % Upper and lower envelopes
    env = normalize(up - lo,'range');                       % Amplitude of the envelope
    env_th = otsuthresh(histcounts(env));                   % Threshold (based on Otsu method). Can be set at 0.35
    wBeats = movsum(env>env_th,2*Fs_Hi)>Fs_Hi/5;            % Euristic criterion for flight detection
    a_mdg = interp1(t_Hi,movmean(abs(a_abs_NP-1),Fs_Hi*0.5),t,'linear',0);  % mean deviation from g in a 0.5s integration window
    immobile = a_mdg<0.04;                                  % Immobility vector
    wb_phase = interp1(t_Hi,angle(hilbert(a_flt_NP)),t,'linear',0);  % Wingbeat phase at t sampling
    %wb_phase = angle(hilbert(a_flt_NP));
    
    %=== Calculate characteristic frequency of wingbeat signal
    n = 2^nextpow2(numel(a_flt_NP(bflying)));           % Optimal number of points for FFT
    Y = fft(a_flt_NP(bflying),n);                       % Calculate fft
    P = abs(Y/n).^2;                                    % Power density at all frequences (positive and negative)
    f = Fs_Hi*(0:n/2)/n;                                % Fs is the sampling frequency
    PSD = P(1:n/2+1);                                   % Power spectral density
    sm_PSD = smoothdata(PSD,'movmedian',n/Fs_Hi*3);     % Smoothed
    [~,max_loc] = max(sm_PSD); f_wBeats = f(max_loc);   % Characteristic frequency of wingbeat signal
    %     plot(f,sm_PSD,'LineWidth',2,'Color','r'); set(gca, 'YScale', 'log');
    
    %=== Calculate Population Spike Density
    all_s = sort(vertcat(s{:,1}));
    all_rate = kernel_rate_AF_v1(all_s,0.1,t)/n_cells;
    
    %=== Calculate Population Spike Density for putative interneurons
    all_s_FS = sort(vertcat(s_FS{:,1}));
    all_rate_FS = kernel_rate_AF_v1(all_s_FS,0.1,t)/n_cells_FS;
    
    %=== Cluster Positional Data and find feeders
    r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64];                                       % Feeder coordinates
    if strcmp(unique_ID{1,1},'Dataset_2')
        r_fd = [-3,0.82,2; 100,100,100];                                            % Correct coordinates for Field Station
    end
    X = r(~bflying & v_abs<0.5,:);                                                  % Consider only stationary epochs
    [~,centroid,~,~,~] = Cluster3D_AF_v1(X,30*Fs,0.3,100);                          % Perform clustering
    for i=1:2                                                                       % Correct feeder position(s)
        [feeder_distance,fd_idx] = min(vecnorm(r_fd(i,:)-centroid,2,2));
        if feeder_distance<0.2,r_fd(i,:) =  centroid(fd_idx,:);end                  % Do not correct if further than 20 cm
    end
    bat_on_feeder = vecnorm(r_fd(1,:)-r,2,2)<0.3 | vecnorm(r_fd(2,:)-r,2,2)<0.3;    % Bat on feeder
    
    %=== Transform start/stop samples of flight to the IMU sampling
    f_smp_Hi = f_smp;
    for i=1:f_num
        if t(f_smp(i,1))>=NP_imu.t(1) && t(f_smp(i,2))<=NP_imu.t(end)
            f_smp_Hi(i,1) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,1)),0.1);
            f_smp_Hi(i,2) = knnsearch_fast_AF_v0(NP_imu.t',t(f_smp(i,2)),0.1);
        else
            f_smp_Hi(i,1) = NaN;
            f_smp_Hi(i,2) = NaN;
        end
    end
    f_smp_Hi = f_smp_Hi(~isnan(f_smp_Hi(:,1)),:);    f_num_Hi = size(f_smp_Hi,1);
    
    %=== Transform the s cell into a cell with the t-samples of each spike
    valid_s = find(all_s>0 & all_s<t(T));
    st1 = all_s(valid_s(1)  );
    st2 = all_s(valid_s(end));
    ss1 = knnsearch_fast_AF_v0(t,st1,1);
    ss2 = knnsearch_fast_AF_v0(t,st2,1);
    m_slope = (ss2-ss1)/(st2-st1);
    q_intcp = ss1-m_slope*st1;
    s_smp = s;           % Single units' spikes samples
    for nc = 1:n_cells
        s_smp{nc,1} = round(s{nc,1}*m_slope+q_intcp);
        for n = 1:n_rep;    s_smp{nc,n+1} = mod(s_smp{nc,1}+smpls_to_shift(n),T); end        % Shuffling
    end
    
    %=== Create structure for storing miscellanea variables
    MSC = struct();
end

%% CLUSTER ANALYSIS
for hide=1
    
    max_curv = NaN(numel(imp_bat_clusters),1);
    gaussR2 = NaN(numel(imp_bat_clusters),1);
    figure('units','normalized','outerposition',[0 .3 1 .4],'Visible',fig_visibility);
    tiledlayout(1,numel(imp_bat_clusters),'TileSpacing','tight');
    for j = setdiff(imp_bat_clusters,1)
        
        %=== Average trajectory shape 3D for this cluster
        flight_pos = f_clus.pos(:,:,f_clus.id==j);
        flight_cell = mat2cell(permute(flight_pos,[2 1 3]),size(flight_pos,2),3,ones(1,size(flight_pos,3)));
        flight_rsmp = cell2mat(cellfun(@(x) interp1(linspace(0,1,numel(x(~isnan(x(:,1)),1))),x(~isnan(x(:,1)),:),linspace(0,1,100)),flight_cell,'UniformOutput',false));
        mean_path3D = squeeze(mean(flight_rsmp,3));
        
        %=== Average trajectory shape 1D for this cluster
        fc_lin = {f_clus.lin_tr{1,f_clus.id == j}}';                                                    % 1D normalized trajectories
        fc_len = f_clus.length(1,f_clus.id == j)';                                                      % Trajectory lenght
        fc_all = cellfun(@(x,y) y.*x,fc_lin,num2cell(fc_len),'UniformOutput',false);                    % Multiply 1D trajectories by their lenght
        fc_min = cellfun(@(x) x(1:min(cellfun(@(x) size(x,1), fc_lin)))',fc_all,'UniformOutput',false); % Keep the first n points that are in common between all the trajectories in the cluster
        mean_path1D = mean(cell2mat(fc_min),1);
        mean_time1D = [1:numel(mean_path1D)]./Fs;
        
        %=== Calculate curvature
        norm_grad3D = diff(mean_path3D,1,1)./vecnorm(diff(mean_path3D,1,1),2,2);
        curv3D = vecnorm(diff(norm_grad3D,1,1),2,2);
        max_curv(j) = max(curv3D);
        
        %=== Calculate fit R2 with 1 or 2 gaussians
        mean_v1D = diff(mean_path1D)./diff(mean_time1D);
        [fg1,gof1] = fit(mean_time1D(2:end)',mean_v1D','gauss1','StartPoint',[6 1 1]);
        [fg2,gof2] = fit(mean_time1D(2:end)',mean_v1D','gauss2','StartPoint',[6 1 1 6 1 1]);
        %plot(fg1,mean_time1D(2:end),mean_v1D);  hold on;    plot(fg2,mean_time1D(2:end),mean_v1D);
        gaussR2(j) = gof1.rsquare;
        
        %=== Visualize cluster
        nexttile;   plot3(mean_path3D(:,1),mean_path3D(:,2),mean_path3D(:,3),'-','LineWidth',1,'Color', col_clus(j,:));
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   zlim(r_lim(3,:));  view(0,90);
        xlabel(['x, Cluster', num2str(j)]);
        axis equal; ylabel(['Max Curvature: ',num2str(max_curv(j),2),'  R2 Gauss: ',num2str(gaussR2(j),2)]);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        
    end
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
end

%% FIND LONG PAUSES BETWEEN FLIGHTS
for hide=1
    
    %=== Find the median interflight interval (or the max of the ifi distribution)
    ifi = diff(t(f_smp(:,1)));
    %[N_ifi,ifi_bins] = histcounts(ifi,[0:3:600]);
    %ifi_ctrs = ifi_bins(1:end-1)+1.5;
    %[~,max_ind] = max(N_ifi);
    %opt_ifi = ifi_ctrs(max_ind);
    opt_ifi = median(ifi);
    
    %=== Pauses between bouts are separated by 2*opt_ifi
    sl_smp_1 = f_smp(find(ifi>2*opt_ifi),2);
    sl_smp_2 = f_smp(find(ifi>2*opt_ifi)+1,1);
    
    %=== Define the sleep vector
    bsleeping = zeros(T,1);
    for ff = 1:numel(sl_smp_1)
        bsleeping(sl_smp_1(ff):sl_smp_2(ff)) = 1;
    end
    
end

%% ECHOLOCATION ANALYSIS
for hide=1
    
    %=== Raster plots of echolocation clicks
    figure('units','normalized','outerposition',[0.2 0.3 0.3 0.5],'Visible',fig_visibility);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   Raster_AF_v3(Detected_Clicks.times,t(f_smp(:,1)),[],[],[-1 1],'Clicks-A-Tko',1,1);
    nexttile;   Raster_AF_v3(Detected_Clicks.times,t(f_smp(:,2)),[],[],[-1 1],'Clicks-A-Lnd',1,1);
    nexttile;   Raster_TimeWrap_AF_v1(Detected_Clicks.times,t(f_smp(:,1)),t(f_smp(:,2)),[],[],1,'Clicks-T-Wrp',1);
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Phase of the wingbeat
    [click_phase,~] = spike_phase_AF_v0(Detected_Clicks.times,a_flt_NP,t_Hi);
    polarhistogram(click_phase,90,'Normalization','probability','DisplayStyle','stairs','EdgeColor','k');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Define an echolocation rate vector, on the same time base of the accelerometer
    click_rate = kernel_rate_AF_v1(Detected_Clicks.times,0.05,t_Hi);
    %plot(t_Hi,normalize(click_rate,'range')); hold on; plot(t_Hi,a_flt_NP);
    
    %=== Quantify acceleration and click rate during rest vs during flight and their profiles
    aa1 = [];   aa2 = [];   cr1 = [];   cr2 = [];
    interval = [round(-0.5*Fs_Hi):round(+0.5*Fs_Hi)-1];
    for i=1:f_num_Hi
        aa1 = [aa1, abs(a_abs_NP(f_smp_Hi(i,1)+interval)-1)];
        aa2 = [aa2, abs(a_abs_NP(f_smp_Hi(i,2)+interval)-1)];
        cr1 = [cr1, click_rate(f_smp_Hi(i,1)+interval)'];
        cr2 = [cr2, click_rate(f_smp_Hi(i,2)+interval)'];
    end
    MSC.aa1 = mean(aa1,2);  MSC.aa2 = mean(aa2,2);                                          % Accelerometer profile around takeoff/landing
    MSC.cr1 = mean(cr1,2);  MSC.cr2 = mean(cr2,2);                                          % Echolocation profile around takeoff/landing
    MSC.tme = interval'/Fs_Hi;                                                              % Time interval around takeoff/landing
    MSC.acc_distribution = abs(a_abs_NP(~bflying_Hi)-1);                                    % Accelerometer data during rest (useful for distribution)
    MSC.abs_acc = [mean(abs(a_abs_NP(bflying_Hi)-1)),mean(abs(a_abs_NP(~bflying_Hi)-1))];   % Mean accelerometer during flight vs rest
    MSC.echo_rate = [mean(click_rate(bflying_Hi)),mean(click_rate(~bflying_Hi))];           % Mean click rate during flight vs rest
    
end

%% PLACE CELLS ON FLIGHT CLUSTERS
for hide=1
    disp('Extracting Place Cells...');
    for nc=1:n_cells
        %=== Loop across fligth clusters (excluding cluster 1)
        for j = setdiff(imp_bat_clusters,1)
            
            smp1_clus = f_clus.strt_frame(f_clus.id == j)';   % Takeoff sample
            smp2_clus = f_clus.stop_frame(f_clus.id == j)';   % Landing sample
            cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
            smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
            n_fclus = numel(smp1_clus);                                         % Number of surviving flights in the cluster
            d_fclus = mean(t(smp2_clus)-t(smp1_clus));                          % Mean duration of the flight
            
            %=== Calculate min samples for optimal shuffling
            min_sh = round(mean(smp2_clus-smp1_clus));
            
            %=== Linearize flights from start to stop, stack them
            fc_lgt = f_clus.length(f_clus.id == j);       % Length in m of the paths
            fc_lin = {f_clus.lin_tr{1,f_clus.id == j}}';  % Flights linearized between 0 and 1
            flight_pos_1D = []; batinclus = zeros(T,1);
            for m = 1:n_fclus
                flight_pos_1D = [flight_pos_1D; fc_lin{m, 1}];          %=== Stack together flights, parametrized between 0 and 1
                batinclus(smp1_clus(m,:):smp2_clus(m,:)) = 1;           %=== Define logical vector when bat is in cluster [0,T]
            end
            
            %=== Calculate linearized trajectories and average 3d path (this is only for calculating correlation and plotting)
            flight_pos = f_clus.pos(:,:,f_clus.id==j);
            flight_pos_rsh = reshape(flight_pos,3,[]);
            flight_pos_rsh = flight_pos_rsh(:,all(~isnan(flight_pos_rsh),1));
            oned_edgs = [0:bin_size_1D:mean(fc_lgt)];
            oned_ctrs = oned_edgs(1:end-1)+diff(oned_edgs)/2;
            flight_cell = mat2cell(permute(flight_pos,[2 1 3]),size(flight_pos,2),3,ones(1,size(flight_pos,3)));
            flight_rsmp = cell2mat(cellfun(@(x) interp1(linspace(0,1,numel(x(~isnan(x(:,1)),1))),x(~isnan(x(:,1)),:),linspace(0,1,100)),flight_cell,'UniformOutput',false));
            length_rsmp = cell2mat(cellfun(@(x) interp1(linspace(0,1,numel(x(~isnan(x(:,1)),1))),x(~isnan(x(:,1)),:),linspace(0,1,100)),fc_lin','UniformOutput',false)');
            
            %=== Get the spikes happening during cluster j
            [num_spikes_per_flight,s_flight] = count_spikes_AF_v1(s{nc,1},t,[smp1_clus smp2_clus+1]);
            s_temp=histcounts(s_flight,t);                                         % Number of above spikes in each sample [0,T]
            s_cut = s_temp(batinclus==1)';                                         % Number of above spikes in each valid flight j sample
            s_pos = [];
            for nr=1:max(s_cut)                                                             % Get the linearized positions associated with those spikes
                s_pos = [s_pos; flight_pos_1D(s_cut>0,:)];
                s_cut = s_cut-1;
            end
            spikes_pos_1D = s_pos;                                                          % Store the data
            
            %=== Initialize and Calculate place properties
            map_1 = NaN(n_rep+1,numel(bin_size_1D:bin_size_1D:mean(fc_lgt)));SI_1 = NaN(1,n_rep+1);SP_1 = NaN(1,n_rep+1);
            prob_x = [];    % Probability of being in spatial bin x
            [map_1(1,:),SI_1(1),SP_1(1),prob_x] = RateMap_AF_v4(spikes_pos_1D,flight_pos_1D,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
            %=== Surrogate generation, only if a minimum of 4 flights with spikes
            if sum(num_spikes_per_flight>0)>(min_flights_with_spikes-1)
                spikes_pos_sh1 = repmat(spikes_pos_1D,1,1,n_rep+1);
                s_cut = s_temp(batinclus==1)';
                for n = 2:n_rep+1
                    s_cut_sh = circshift(s_cut,randi([min_sh size(s_cut,1)-min_sh],1),1); % Circshift
                    s_pos = [];
                    for nr=1:max(s_cut)                                 % Get the positions
                        s_pos = [s_pos; flight_pos_1D(s_cut_sh>0,:)];
                        s_cut_sh = s_cut_sh-1;
                    end
                    spikes_pos_sh1(1:size(s_pos,1),:,n) = s_pos;                      % Store the surrogate
                    
                    [map_1(n,:),SI_1(n),SP_1(n),~] = RateMap_AF_v4(spikes_pos_sh1(:,:,n),flight_pos_1D,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
                end
            end
            
            %=== Split trials into even and odd (or just random half)
            idx_lo = [1:2:n_fclus]';                 % sort(datasample([1:size(F_temp,1)],ceil(size(F_temp,1)/2),'Replace',false))';
            idx_hi = setdiff([1:n_fclus]',idx_lo);
            
            %=== Calculate correlation between paths falling into those 2 categories
            path_corrOdEv = trace(corr(mean(flight_rsmp(:,:,idx_lo),3),mean(flight_rsmp(:,:,idx_hi),3)))/3;
            
            %=== Map Calculation for odd/even trials
            [spikes_pos_1DOd,flight_pos_1DOd] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_lo',t);
            [map_1Od,~,~,~] = RateMap_AF_v4(spikes_pos_1DOd,flight_pos_1DOd,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            [spikes_pos_1DEv,flight_pos_1DEv] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_hi',t);
            [map_1Ev,~,~,~] = RateMap_AF_v4(spikes_pos_1DEv,flight_pos_1DEv,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
            %=== Calculate correlation between maps
            map_corrOdEv = corr(map_1Ev',map_1Od','Type','Spearman');
            
            %=== Calculate distance between maps
            map_distOdEv = pdist2(map_1Ev,map_1Od)/size(oned_ctrs,2);
            
            %=== Split trials into first half and second half
            idx_lo = [1:floor(n_fclus/2)]';
            idx_hi = setdiff([1:n_fclus]',idx_lo);
            
            %=== Calculate correlation between paths falling into those 2 categories
            path_corr1h2h = trace(corr(mean(flight_rsmp(:,:,idx_lo),3),mean(flight_rsmp(:,:,idx_hi),3)))/3;
            
            %=== Map Calculation for first half and second half
            [spikes_pos_1D1h,flight_pos_1D1h] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_lo',t);
            [map_11h,~,~,~] = RateMap_AF_v4(spikes_pos_1D1h,flight_pos_1D1h,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            [spikes_pos_1D2h,flight_pos_1D2h] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_hi',t);
            [map_12h,~,~,~] = RateMap_AF_v4(spikes_pos_1D2h,flight_pos_1D2h,Fs,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
            %=== Calculate correlation between maps
            map_corr1h2h = corr(map_12h',map_11h','Type','Spearman');
            
            %=== Calculate distance between maps
            map_dist1h2h = pdist2(map_12h,map_11h)/size(oned_ctrs,2);
            
            %=== Save 1D field information
            NP_unit(nc).f_clus(j).cell_id = nc;
            NP_unit(nc).f_clus(j).id = j;
            NP_unit(nc).f_clus(j).n = size(flight_pos,3);
            NP_unit(nc).f_clus(j).f_lenght = mean(fc_lgt);
            NP_unit(nc).f_clus(j).f_duration = d_fclus;
            NP_unit(nc).f_clus(j).SI = SI_1;
            NP_unit(nc).f_clus(j).SP = SP_1;
            NP_unit(nc).f_clus(j).map = map_1;
            NP_unit(nc).f_clus(j).binC = oned_ctrs;
            NP_unit(nc).f_clus(j).prob_x = prob_x;
            NP_unit(nc).f_clus(j).prc_lg = mean(length_rsmp,1)*mean(fc_lgt);
            NP_unit(nc).f_clus(j).prc_3d = mean(flight_rsmp,3);
            NP_unit(nc).f_clus(j).Q_corr = median(f_clus.corr(1,f_clus.id==j));
            NP_unit(nc).f_clus(j).Q_dist = median(f_clus.dist(1,f_clus.id==j));
            NP_unit(nc).f_clus(j).corr_map_OdEv = map_corrOdEv;
            NP_unit(nc).f_clus(j).dist_map_OdEv = map_distOdEv;
            NP_unit(nc).f_clus(j).corr_pth_OdEv = path_corrOdEv;
            NP_unit(nc).f_clus(j).corr_map_1h2h = map_corr1h2h;
            NP_unit(nc).f_clus(j).dist_map_1h2h = map_dist1h2h;
            NP_unit(nc).f_clus(j).corr_pth_1h2h = path_corr1h2h;
            NP_unit(nc).f_clus(j).sum_spk = sum(num_spikes_per_flight);
            
            %=== Calculate p_value SI
            non_nan = nnz(~isnan(SI_1(2:end)));
            if isnan(SI_1(1)), p_val_l = NaN;   p_val_r = NaN;
            else, p_val_r =  nnz(SI_1(2:end)>SI_1(1))/non_nan;p_val_l =  nnz(SI_1(2:end)<SI_1(1))/non_nan;end
            NP_unit(nc).f_clus(j).SI_value = SI_1(1,1);
            NP_unit(nc).f_clus(j).SI_p_val = p_val_r;
            
            %=== Add the firing rate at the peak and the location of the peak
            NP_unit(nc).f_clus(j).field = map_1(1,:);
            [NP_unit(nc).f_clus(j).peakHz,NP_unit(nc).f_clus(j).field_loc] = max(map_1(1,:));
            
            %=== Get the number of fields and width of the place field
            warning('off');
            [heights,locs,width] = findpeaks(NP_unit(nc).f_clus(j).field,linspace(0,mean(fc_lgt),numel(NP_unit(nc).f_clus(j).field)),'MinPeakDistance',0.3*mean(fc_lgt),'MinPeakHeight',0.5*NP_unit(nc).f_clus(j).peakHz,'SortStr','descend');
            NP_unit(nc).f_clus(j).n_fields = numel(locs);
            if numel(locs)>1
                sh = sort(heights,'descend');
                NP_unit(nc).f_clus(j).sff = sh(2)/sh(1);    % Single field fraction
            else
                NP_unit(nc).f_clus(j).sff = 0;
            end
            
            if numel(locs),NP_unit(nc).f_clus(j).f_width = width(1);
            else,NP_unit(nc).f_clus(j).f_width = 0;end
            warning('on');
            
            %=== Get the field shape in percentages of flight (Interpolate between takeoff and landing, 100 pts)
            NP_unit(nc).f_clus(j).map_interp = interp1(NP_unit(nc).f_clus(j).binC',NP_unit(nc).f_clus(j).map(1,:)',linspace(NP_unit(nc).f_clus(j).binC(1),NP_unit(nc).f_clus(j).binC(end),100)','linear')';
            [~,NP_unit(nc).f_clus(j).phase_max] = max(NP_unit(nc).f_clus(j).map_interp);  %=== Get the location of the maximum, in terms of flight phase
            
        end
    end
    
    %=== Rearrange Place Cell Data
    NP_unitOnClus = cell(1,n_surv_clusters);
    for j = setdiff(imp_bat_clusters,1)
        NP_unitOnClus{1, j} = struct(); % Initialize as struct array
        
        %=== Store features of the unit within a cluster of interest
        for nc = 1:n_cells
            NP_unitOnClus{1, j}(nc).cell_id = nc;
            NP_unitOnClus{1, j}(nc).clus_id = j;
            NP_unitOnClus{1, j}(nc).f_lenght = NP_unit(nc).f_clus(j).f_lenght;
            NP_unitOnClus{1, j}(nc).f_duration = NP_unit(nc).f_clus(j).f_duration;
            NP_unitOnClus{1, j}(nc).fr = NP_unit(nc).fr;
            NP_unitOnClus{1, j}(nc).SI = NP_unit(nc).f_clus(j).SI_value;
            NP_unitOnClus{1, j}(nc).p_val = NP_unit(nc).f_clus(j).SI_p_val;
            NP_unitOnClus{1, j}(nc).spkPerflight = NP_unit(nc).f_clus(j).sum_spk/NP_unit(nc).f_clus(j).n;
            NP_unitOnClus{1, j}(nc).plc_map = NP_unit(nc).f_clus(j).map(1,:)';
            NP_unitOnClus{1, j}(nc).plc_ctr = NP_unit(nc).f_clus(j).binC';
            NP_unitOnClus{1, j}(nc).prob_x = NP_unit(nc).f_clus(j).prob_x';
            NP_unitOnClus{1, j}(nc).corr_m = NP_unit(nc).f_clus(j).corr_map_OdEv;
            NP_unitOnClus{1, j}(nc).dist_m = NP_unit(nc).f_clus(j).dist_map_OdEv;
            NP_unitOnClus{1, j}(nc).stab_m = NP_unit(nc).f_clus(j).corr_map_1h2h;
            NP_unitOnClus{1, j}(nc).stbd_m = NP_unit(nc).f_clus(j).dist_map_1h2h;
            NP_unitOnClus{1, j}(nc).peakHz = NP_unit(nc).f_clus(j).peakHz;
            NP_unitOnClus{1, j}(nc).field_loc = NP_unit(nc).f_clus(j).field_loc;
            NP_unitOnClus{1, j}(nc).field_loc_m = NP_unit(nc).f_clus(j).field_loc*bin_size_1D;
            NP_unitOnClus{1, j}(nc).n_fields = NP_unit(nc).f_clus(j).n_fields;
            NP_unitOnClus{1, j}(nc).f_width = NP_unit(nc).f_clus(j).f_width;
            NP_unitOnClus{1, j}(nc).phase_max = NP_unit(nc).f_clus(j).phase_max;
            NP_unitOnClus{1, j}(nc).sff = NP_unit(nc).f_clus(j).sff;
            NP_unitOnClus{1, j}(nc).asymm_frct = sum(NP_unitOnClus{1, j}(nc).plc_map(1:round(NP_unitOnClus{1, j}(nc).field_loc)))./sum(NP_unitOnClus{1, j}(nc).plc_map);
            NP_unitOnClus{1, j}(nc).map_interp = NP_unit(nc).f_clus(j).map_interp;
            
        end
    end
end

%% AUTOCORRELATION AND THETA MODULATION INDEX
for hide=1
    
    %=== Params
    bin_TMI = 0.01;
    theta_range = [5 11];
    max_lag_TMI = 0.5;
    
    %=== Autocorrelation and TMI for principal cells
    for nc=1:n_cells
        
        %=== Initialize Autocorrelation, TMI and p value
        NP_unit(nc).AC = [];
        NP_unit(nc).AC_bins = [];
        NP_unit(nc).TMI = NaN(1,n_rep_sweeps);
        NP_unit(nc).p_val_TMI = NaN;
        
        
        %=== Loop across shuffles
        for jj=1:n_rep_sweeps+1
            [~,s_flight] = count_spikes_AF_v1(s{nc,jj},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
            if ~isempty(s_flight)
                if jj==1
                    [NP_unit(nc).AC,NP_unit(nc).AC_bins] = cross_correlogram_AF_v0(s_flight,s_flight,max_lag_TMI,bin_TMI);
                end
                NP_unit(nc).TMI(jj) = Theta_mod_idx_AF_v0(s_flight,bin_TMI,theta_range,max_lag_TMI);     % Theta modulation index
            end
        end
        
        %=== Calculate p value
        if ~isnan(NP_unit(nc).TMI(1)), NP_unit(nc).p_val_TMI = sum(NP_unit(nc).TMI>NP_unit(nc).TMI(1))/sum(~isnan(NP_unit(nc).TMI));end
        
    end
    
    %=== Plot the 20 best cells
    figure('units','normalized','outerposition',[0.1 0.3 0.8 0.4],'Visible',fig_visibility);
    tiledlayout(2,10,'TileSpacing','tight');
    [~,tmp_idx] = sort([NP_unit.p_val_TMI]);
    for i=1:20
        nexttile;
        area(NP_unit(tmp_idx(i)).AC_bins,NP_unit(tmp_idx(i)).AC,'FaceColor','k');   yticks([]); xlabel('Time lag (s)'); title(['TMI: ',num2str(NP_unit(tmp_idx(i)).TMI(1),2), ', p = ',num2str(NP_unit(tmp_idx(i)).p_val_TMI,1)]);
        hold on; plot(repmat(1/f_wBeats*[-3:3],2,1),repmat(ylim,7,1)','r--');
        ylim([0 1.3*max(NP_unit(tmp_idx(i)).AC(NP_unit(tmp_idx(i)).AC_bins>bin_TMI))]);
        
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% LOOK AT POPULATION RASTER PLOT
for hide=1
    
    %=== Plot raster with sorted units (by firing frequency), LFP and Spike Density
    figure('units','normalized','outerposition',[0 .3 1 .6],'Visible',fig_visibility);
    tiledlayout(13,1,'TileSpacing','none');
    cx(1) = nexttile(1,[4 1]);
    for nc= 1:n_cells
        plot(s{rate_sorted_ids(nc)}, nc*ones(size(s{rate_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(5,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(6,[2 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_lp);  xticks([]); legend('All','Theta');  ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
    cx(4) = nexttile(8,[2 1]);
    imagesc([t(1),t(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[.5 100])); shading interp;  colormap(hot);
    set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([0.5 50]);  yticks([1 5 10 20 50]);    ylabel('Freq (Hz)');
    %cx(5) = nexttile(10,[2 1]); plot(t,smoothdata(all_rate,'movmean',Fs*1),'k');   ylabel('Spike Density');    ylim([0 prctile(all_rate,99.99)]);   set(gca,'TickLength',[0 0]);
    %cx(5) = nexttile(10,[2 1]); plot(t,mean(Rate./max(Rate),2));   ylabel('Spike Density');    set(gca,'TickLength',[0 0]);
    cx(5) = nexttile(10,[2 1]); plot(t,smoothdata(mean(zscore(Rate,1),2),'movmedian',5*Fs));   ylabel('Spike Density');    set(gca,'TickLength',[0 0]);
    cx(6) = nexttile(12,[2 1]); plot(t,smoothdata(PSD_low, 'movmedian',Fs*5),'k');    ylabel('LP Power');   hold on;
    area(t,bsleeping*max(PSD_low),0,'FaceColor',[1 0 0],'FaceAlpha',0.1,'LineStyle','none');
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %     figure
    %     [~,all_s_imm] = count_spikes_AF_v1(all_s,t,[knnsearch(t,4943) knnsearch(t,5591)]);
    %     [ACS,ACS_bins] = cross_correlogram_AF_v0(all_s_imm,all_s_imm,10,0.05);
    %     bar(ACS_bins,ACS,'FaceColor','k');  xlim([0 10]);
    
    %=== Plot raster with sorted units (II VERSION)
    %=== Make bflying_c vector
    bflying_c = zeros(T,n_surv_clusters);
    for j=1:n_surv_clusters
        smp1_clus = f_clus.strt_frame(f_clus.id == j)';   % Takeoff sample
        smp2_clus = f_clus.stop_frame(f_clus.id == j)';   % Landing sample
        for ff = 1:numel(smp1_clus)
            bflying_c(smp1_clus(ff):smp2_clus(ff),j) = 1;
        end
    end
    figure('units','normalized','outerposition',[0 .3 1 .6],'Visible',fig_visibility);
    tiledlayout(9,1,'TileSpacing','none');
    dx(1) = nexttile(1,[4 1]);
    for nc= 1:n_cells
        plot(s{rate_sorted_ids(nc)}, nc*ones(size(s{rate_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    dx(2) = nexttile(5,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
    dx(3) = nexttile(6,[2 1]);   plot(t_Hi,LFP_mn); hold on;    xticks([]);  ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
    dx(4) = nexttile(8,[2 1]);
    for j=1:n_surv_clusters
        area(t,bflying_c(:,j),0,'FaceColor',col_clus(j,:),'FaceAlpha',1,'LineStyle','none');  hold on;
    end
    linkaxes(dx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% SHOW EXAMPLE PLACE CELLS AND LOOK AT REPLAY
for hide=1
    
    for j = setdiff(imp_bat_clusters,1)
        
        %=== Average trajectory shape 3D for this cluster
        flight_pos = f_clus.pos(:,:,f_clus.id==j);
        flight_cell = mat2cell(permute(flight_pos,[2 1 3]),size(flight_pos,2),3,ones(1,size(flight_pos,3)));
        flight_rsmp = cell2mat(cellfun(@(x) interp1(linspace(0,1,numel(x(~isnan(x(:,1)),1))),x(~isnan(x(:,1)),:),linspace(0,1,100)),flight_cell,'UniformOutput',false));
        mean_path3D = squeeze(mean(flight_rsmp,3));
        dF = [vecnorm(mean_path3D(1,:)-r_fd,2,2),vecnorm(mean_path3D(end,:)-r_fd,2,2)];                 % Takeoff/Landing distances from F1,F2
        
        %=== Average trajectory shape 1D for this cluster
        fc_lin = {f_clus.lin_tr{1,f_clus.id == j}}';                                                    % 1D normalized trajectories
        fc_len = f_clus.length(1,f_clus.id == j)';                                                      % Trajectory lenght
        fc_all = cellfun(@(x,y) y.*x,fc_lin,num2cell(fc_len),'UniformOutput',false);                    % Multiply 1D trajectories by their lenght
        fc_min = cellfun(@(x) x(1:min(cellfun(@(x) size(x,1), fc_lin)))',fc_all,'UniformOutput',false); % Keep the first n points that are in common between all the trajectories in the cluster
        mean_path1D = mean(cell2mat(fc_min),1);
        mean_time1D = [1:numel(mean_path1D)]./Fs;
        
        %=== Get the subtable of cells for template matching and cells for decoding
        NP_table = struct2table(NP_unitOnClus{1, j});
        
        %=== For template matching
        NP_table.place_cond = NP_table.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
            NP_table.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
            NP_table.stab_m>0.4 &...      % Min Stability (DEF: .4)
            NP_table.sff<0.5;             % Min Peakyness (DEF: .7)
        NP_subtable = NP_table(NP_table.place_cond,:);
        n_place_cells = size(NP_subtable,1);
        meters_per_cell = range(NP_subtable.field_loc_m)/n_place_cells;  % Average space covered by one place cell (useful later)
        flight_bin_centers = NP_subtable.plc_ctr{1,1};                   % Positions at bin centers (useful later)
        
        %=== Extract some variables for plotting
        smp1_clus = f_clus.strt_frame(f_clus.id == j)';                     % Takeoff sample
        smp2_clus = f_clus.stop_frame(f_clus.id == j)';                     % Landing sample
        cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
        smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
        n_fclus = numel(smp1_clus);                                         % Number of surviving flights in the cluster
        id = find(f_clus.id==j);                                            % ids of the flight clusters
        avg_takeoff = mean(squeeze(f_clus.ds_pos(:,1,id)),2);               % Average takeoff position
        avg_landing = mean(squeeze(f_clus.ds_pos(:,end,id)),2);             % Average landing position
        d_fclus = mean(t(smp2_clus)-t(smp1_clus));                          % Mean duration of the flight
        all_s_plc = sort(vertcat(s{NP_subtable.cell_id,1}));                % Define rate by pooling spikes from place cells
        all_rate_plc = kernel_rate_AF_v1(all_s_plc,spk_dsty_smooth_t,t)/n_place_cells;
        t_ic = [-1 1];                                                      % Time for plotting around the flight
        smooth_f = [.01 1];                                                 % [space bin,time bin]
        
        %=== Target sequence determination
        %sorted_tko_plcCells = intersect(sorted_tko,NP_subtable.cell_id,'stable');   % Restrict to place cells only
        %NP_subtable_sorted = sortrows(NP_subtable,'field_loc');
        NP_subtable_sorted = sortrows(NP_subtable,'phase_max');
        sorted_tko_plcCells = NP_subtable_sorted.cell_id;
        
        %=== Get the spike times happening around flight (plus t_ic)
        s_flight = cell(n_place_cells,1);       % Spikes during the flight cluster (sorted)
        for nc = 1:n_place_cells,[~,s_flight{nc,1}] = count_spikes_AF_v1(s{sorted_tko_plcCells(nc),1},t,[smp1_clus smp2_clus]);end
        s_flight1 = cell(n_cells,1);            % Spikes during the flight cluster (unsorted, extended)
        for nc = 1:n_cells,[~,s_flight1{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);end
        s_smp_flight = cell(n_place_cells,n_rep);% Spike samples during the flight cluster (sorted)
        for nr= 1:n_rep+1
            for nc = 1:n_place_cells
                [~,~,valid_samples] = count_spikes_AF_v2(s{sorted_tko_plcCells(nc),nr},t,[smp1_clus smp2_clus]);
                s_smp_flight{nc,nr} = s_smp{sorted_tko_plcCells(nc),nr}(valid_samples);
            end
        end
        
        %=== Visualize cluster
        figure('units','normalized','outerposition',[0.35 .3 0.2 .4],'Visible',fig_visibility);
        plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   zlim(r_lim(3,:));  view(0,90);
        xlabel(['x, Cluster' num2str(j) ' (' num2str(n_fclus) ' flights),']);    ylabel('y');    hold on;
        for ii=1:n_fclus
            plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(j,:));
        end
        textscatter(avg_takeoff(1),avg_takeoff(2),"Take-off");     hold off;   axis equal;
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(6,4*ceil(n2show/6),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable_sorted.cell_id(ncc);
            
            %=== Calculate linearized trajectories and average 3d path
            flight_pos = f_clus.pos(:,:,f_clus.id==j);
            flight_pos_rsh = reshape(flight_pos,3,[]);
            %flight_pos_rsh = flight_pos_rsh(:,all(~isnan(flight_pos_rsh),1));
            
            %=== Find the closest positional sample to each spike (plotting purposes)
            if ~isempty(s_flight{ncc,1})
                k = round((s_flight{ncc,1}-t(1))*Fs);    k = k(k<T);
                spikes_pos = r(k,:);
            else, spikes_pos = [nan nan nan];end
            
            %=== Plot
            nexttile([1 2]);   plot(flight_pos_rsh(1,:),flight_pos_rsh(2,:),'Color',[.9 .9 .9]);
            hold on;    textscatter(flight_pos_rsh(1,1)-0.05,flight_pos_rsh(2,1),"**");
            scatter(spikes_pos(:,1),spikes_pos(:,2),6,'MarkerFaceColor', 'r','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
            hold off;    axis equal;
            xlim(r_lim(1,:)); ylim(r_lim(2,:)); 
            xticks([]); yticks([]); title(['Unit ' num2str(nc)]);
            nexttile([1 2]);   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],.5,[],1);
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells sorted by peak location
        figure('units','normalized','outerposition',[.3 .1 .1 .35],'Visible',fig_visibility);
        tiledlayout(n_place_cells,1,'TileSpacing','none');
        for ncc=flip(1:n_place_cells)
            nexttile;   Raster_TimeWrap_AF_v1(s{sorted_tko_plcCells(ncc),1},t(smp1_clus),t(smp2_clus),[],[],0.1,[],0);    box('Off');
            if ncc<n_place_cells
                set(gca, 'YTick', [],'YColor', 'none','XLabel',[],'Box','off','XTick', []);
            else
                set(gca, 'YTick', [],'YColor', 'none','XLabel',[],'Box','off');
            end
        end
        xticks([0, 1]); xticklabels({'Takeoff','Landing'});
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot sorted activity around each flight
        figure('units','normalized','outerposition',[0 .2 .8 .3],'Visible',fig_visibility);
        tiledlayout(1,n_fclus,'TileSpacing','none');
        ax = [];
        for i=1:n_fclus
            ax(i) = nexttile;
            for nc = 1:n_place_cells
                plot(s_flight1{sorted_tko_plcCells(nc)}, nc*ones(size(s_flight1{sorted_tko_plcCells(nc)})), 'k|','MarkerSize', max(round(n_place_cells*0.01),1));    hold on;
                xlim([t(smp1_clus(i))+t_ic(1) t(smp1_clus(i))+d_fclus+t_ic(2)]);
            end
            area(t,bflying*n_place_cells,0,'FaceColor',col_clus(j,:),'FaceAlpha',.1,'LineStyle','none');
            yticks([]); xticks([]); ylim([0 n_place_cells]);   box off;
            if i==1, ylabel('Unit #'); yticks([20:20:n_place_cells]); end
        end
        linkaxes(ax,'y');   title(['Activity during cluster ',num2str(j), ', Time to takeoff: [', num2str(t_ic(1)), ',' num2str(d_fclus+t_ic(2),2), '] s']);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Make snake plot
        figure('units','normalized','outerposition',[.5 .3 .1 .4],'Visible',fig_visibility);
        Snake = NP_subtable.map_interp;         Snake = Snake./max(Snake,[],2);
        [~,order] = sort(NP_subtable.phase_max,'descend');
        imagesc(Snake(order,:));    colormap(viridis);  ylabel('1D Field');   xlabel('Flight Phase (%)');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        
        %==================================================== REPLAY ANALYSIS STARTS HERE ========================================================
        %=== Detect candidate events from spike density from place cells only
        trace_for_candidate = all_rate_plc;
        trace_for_candidate(bflying) = 0;
        [candidate_event_vector,~,candidate_event_str,candidate_event_stp,candidate_event_width] = BiLevel_Segm_AF_v0(zscore(trace_for_candidate)',Rate_treshold_SD,Fs,Replay_db_time,Replay_min_dur);
        t_temp= t(candidate_event_str);
        candidate_event_str =     candidate_event_str(candidate_event_width<max_Replay_dur & t_temp>0);
        candidate_event_stp =     candidate_event_stp(candidate_event_width<max_Replay_dur & t_temp>0);
        candidate_event_width = candidate_event_width(candidate_event_width<max_Replay_dur & t_temp>0);
        
        %=== Find Candidate Replays
        cRt = [t(candidate_event_str),t(candidate_event_stp)];
        n_candidate_events = size(cRt,1);
        candidate_perc_act = zeros(n_candidate_events,1);
        candidate_rkcorr_c = zeros(n_candidate_events,1);
        candidate_rkcorr_p = NaN(n_candidate_events,1);
        candidate_rkcorr_v = NaN(n_candidate_events,1);
        candidate_rkcorr_shuffles = NaN(n_candidate_events,100);
        candidate_s_all = cell(n_candidate_events,1);
        candidate_pairs_cell_ids = cell(n_candidate_events,1);
        candidate_pairs_timing = cell(n_candidate_events,1);
        
        %=== Assess event quality
        warning('off');
        for i=1:n_candidate_events
            tL = cRt(i,1);
            tR = cRt(i,2);
            candidate_s = cellfun(@(x) x(x>tL & x<tR), s(sorted_tko_plcCells,1),'UniformOutput', false);   % Extract the spikes happening during the event
            candidate_perc_act(i) =1-sum(cellfun(@isempty,candidate_s))/n_place_cells;
            if  numel(vertcat(candidate_s{:}))>seq_min_spikes && candidate_perc_act(i)>seq_min_fr_active
                [candidate_rkcorr_c(i),candidate_rkcorr_p(i),candidate_rkcorr_shuffles(i,:),candidate_rkcorr_v(i),candidate_pairs_cell_ids(i),candidate_pairs_timing(i)] = spikeSeq_analysis_AF_v1(candidate_s,100);
                candidate_s_all(i) = {candidate_s};
            end
        end
        warning('on');
        
        %=== Calculate pairwise distances between field peaks
        cell_pairs = nchoosek(sorted_tko_plcCells,2);
        pw_distance = cell_pairs(:,1)*0;
        for i=1:numel(pw_distance)
            cell1 = find(NP_subtable.cell_id==cell_pairs(i,1));
            cell2 = find(NP_subtable.cell_id==cell_pairs(i,2));
            pw_distance(i) = NP_subtable.plc_ctr{1,1}(NP_subtable.field_loc(cell2))-NP_subtable.plc_ctr{1,1}(NP_subtable.field_loc(cell1));
        end
        
        %=== Extract echolocation rate and accelerometer around candidate replays plus/minus 0.5s(and a random point as control) 
        interval = [round(-0.5*Fs_Hi):round(+0.5*Fs_Hi)-1];
        smpL_Hi = knnsearch(t_Hi',cRt(:,1));
        smpR_Hi = knnsearch(t_Hi',cRt(:,2));
        smpC_Hi = smpL_Hi+round((smpR_Hi-smpL_Hi)/2);
        smpC_Sf = datasample(intersect(find(~bflying_Hi),[round(0.5*Fs_Hi)+1:(numel(bflying_Hi)-round(0.5*Fs_Hi)-1)]),n_candidate_events,'Replace',false);
        candidate_acc = cell(n_candidate_events,2);
        candidate_ech = cell(n_candidate_events,2);
        for i=1:n_candidate_events
            if smpC_Hi(i)>round(0.5*Fs_Hi) && smpC_Hi(i)<numel(bflying_Hi)-round(0.5*Fs_Hi)
                candidate_acc{i,1} = abs(a_abs_NP(smpC_Hi(i)+interval)'-1);
                candidate_ech{i,1} = click_rate(smpC_Hi(i)+interval);
                candidate_acc{i,2} = abs(a_abs_NP(smpC_Sf(i)+interval)'-1);
                candidate_ech{i,2} = click_rate(smpC_Sf(i)+interval);
            end
        end
%         plot(mean(vertcat(candidate_acc{:,1})));    hold on;    plot(mean(vertcat(candidate_acc{:,2})))
%         plot(mean(vertcat(candidate_ech{:,1})));    hold on;    plot(mean(vertcat(candidate_ech{:,2})))
        
        %=== Define Replay Table for this cluster and add to the general table
        Rpl_table_tmp = table();
        Rpl_table_tmp.unique_ID = repmat(unique_ID,n_candidate_events,1);
        Rpl_table_tmp.clus_id = j*ones(n_candidate_events,1);                                                       % Cluster Id.
        Rpl_table_tmp.flight_dur = NP_table.f_duration(1)*ones(n_candidate_events,1);                               % Flight Duration
        Rpl_table_tmp.flight_len = NP_table.f_lenght(1)*ones(n_candidate_events,1);                                 % Flight Lenght
        Rpl_table_tmp.n_flights = n_fclus*ones(n_candidate_events,1);                                               % Number of flights
        Rpl_table_tmp.bin1Ds = bin_size_1D*ones(n_candidate_events,1);                                              % Bin size for posterior (space)
        Rpl_table_tmp.bin1Dt = (t_bin_dur-t_bin_ovl)*ones(n_candidate_events,1);                                                         % Bin size for posterior (time)
        Rpl_table_tmp.mean_path3D = repmat({mean_path3D},n_candidate_events,1);                                     % Average Path 3D (resampled in 100 pts)
        Rpl_table_tmp.mean_path1D = repmat({mean_path1D},n_candidate_events,1);                                     % Average Path 1D (original sampling)
        Rpl_table_tmp.mean_time1D = repmat({mean_time1D},n_candidate_events,1);                                     % Time for the 1D path
        Rpl_table_tmp.posterior = repmat({mean_time1D},n_candidate_events,1);                                       % Allocate space for the posterior
        Rpl_table_tmp.fdist = repmat({dF},n_candidate_events,1);                                                    % Distance from feeders at takeoff/landing
        Rpl_table_tmp.tL = cRt(:,1);    Rpl_table_tmp.tR = cRt(:,2);                                                % Start Stop times
        Rpl_table_tmp.dt = diff(cRt,1,2);                                                                           % Replay duration
        Rpl_table_tmp.smpL = knnsearch(t,Rpl_table_tmp.tL); Rpl_table_tmp.smpR = knnsearch(t,Rpl_table_tmp.tR);     % Start Stop samples
        Rpl_table_tmp.bat_on_feeder = bat_on_feeder(Rpl_table_tmp.smpL);                                            % If the bat was on the feeder or not
        Rpl_table_tmp.bat_mdg = a_mdg(Rpl_table_tmp.smpL);                                                          % Mean deviation from g (accelerometer)
        Rpl_table_tmp.spikes = candidate_s_all;                                                                     % Spikes within the event
        Rpl_table_tmp.n_place_cells = n_place_cells*ones(n_candidate_events,1);                                     % Number of place cells
        Rpl_table_tmp.f = candidate_perc_act;                                                                       % Fraction of active cells
        Rpl_table_tmp.n = candidate_perc_act*n_place_cells;                                                         % Number of active cells
        Rpl_table_tmp.C = candidate_rkcorr_c;                                                                       % Rank correlation
        Rpl_table_tmp.p = candidate_rkcorr_p;                                                                       % p value
        Rpl_table_tmp.v = meters_per_cell*candidate_rkcorr_v;                                                       % Replay speed
        Rpl_table_tmp.c_ratio =  Rpl_table_tmp.v/(NP_subtable.f_lenght(1)/NP_subtable.f_duration(1));               % Compression Ratio
        Rpl_table_tmp.shuffles = mat2cell(candidate_rkcorr_shuffles,ones(1,numel(candidate_rkcorr_p)),100);         % Shuffled Rank correlation
        Rpl_table_tmp.d1 = vecnorm(r(Rpl_table_tmp.smpL,:)-avg_takeoff',2,2);                                       % Distance from takeoff spot
        Rpl_table_tmp.d2 = vecnorm(r(Rpl_table_tmp.smpL,:)-avg_landing',2,2);                                       % Distance from landing spot
        [~,Rpl_table_tmp.t2tko] = knnsearch(t(smp1_clus),Rpl_table_tmp.tL);                                         % Time from nearest takeoff from that cluster
        [~,Rpl_table_tmp.t2lnd] = knnsearch(t(smp2_clus),Rpl_table_tmp.tL);                                         % Time from nearest landing from that cluster
        Rpl_table_tmp.cell_pairs = candidate_pairs_cell_ids;
        Rpl_table_tmp.pw_timing = candidate_pairs_timing;
        Rpl_table_tmp.pw_distance = repmat({pw_distance},n_candidate_events,1);
        for rr=1:size(Rpl_table_tmp,1)
            Rpl_table_tmp.issleep(rr) = bsleeping(knnsearch_fast_AF_v0(t,Rpl_table_tmp.tL(rr)+Rpl_table_tmp.dt(rr)/2,0.1)); % If happening during sleep or not
        end
        Rpl_table_tmp.sleep_fraction = (sum(bsleeping)/numel(bsleeping))*ones(n_candidate_events,1);                % Sleep fraction
        Rpl_table_tmp.place_phases = repmat({NP_subtable_sorted.phase_max},n_candidate_events,1);                   % Store the place field location (flight phase) for each cell in the sequence
        Rpl_table_tmp.abs_acc = candidate_acc;
        Rpl_table_tmp.ech_rat = candidate_ech;
        
        %=== Define subset of interest
        rng(1);
        cRT_good = find(Rpl_table_tmp.n>4 & Rpl_table_tmp.p<0.05 & Rpl_table_tmp.f>0.1);
        n_event2show = min(100,size(cRT_good,1));
        cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
        Rpl_subtable = Rpl_table_tmp(cRT_good,:);
        
        %=== Plot sorted activity (entire session)
        figure('units','normalized','outerposition',[0 .2 1 .6],'Visible',fig_visibility);
        tiledlayout(7,1,'TileSpacing','none');
        bx(1) = nexttile([2 1]);
        for nc= 1:n_place_cells
            plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'k|','MarkerSize', 5);  hold on;             % Raster for each cell
        end
        area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.2,'LineStyle','none');       % Plot flights
        area(t,bat_on_feeder*n_cells,0,'FaceColor',[0 1 0],'FaceAlpha',0.2,'LineStyle','none'); % Plot when at feeder
        area(t,immobile*n_cells,0,'FaceColor',[1 0 0],'FaceAlpha',0.1,'LineStyle','none');      % Plot when immobile
        plot((t(f_clus.strt_frame(f_clus.id == j)')*[1 1])',(ones(size(f_clus.strt_frame(f_clus.id == j)'))*[0 n_cells])','r');
        ylim([0 n_place_cells]);  ylabel('Unit #');    xticks([]);
        bx(2) = nexttile([1 1]);  plot(t,all_rate_plc,'k');     ylim([0 prctile(all_rate_plc,99.99)]);   set(gca,'TickLength',[0 0]);   ylabel('Spike Density');
        bx(3) = nexttile([1 1]);  area(t,candidate_event_vector);                                        set(gca,'TickLength',[0 0]);   ylabel('Candidate events');
        bx(4) = nexttile([1 1]);  plot(t_Hi,a_abs_NP,'k');                                               set(gca,'TickLength',[0 0]);   ylabel('Accelerometer');
        bx(5) = nexttile([1 1]);  plot(RPL_out.time_vector,RPL_out.LFP_trace,'b');                       set(gca,'TickLength',[0 0]);   ylabel('LFP');
        bx(6) = nexttile([1 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_th);  legend('All','Theta');  ylabel('LFP (norm)');
        %bx(6) = nexttile([1 1]); plot(linspace(0,numel(data)/fs,numel(data)),data);   ylabel('Mic');
        xlabel('Time (s)'); linkaxes(bx,'x');   xlim('tight');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot pairwise intervals between spikes during replay vs pairwise distances between field peaks.
        figure('units','normalized','outerposition',[.2 .3 .25 .33],'Visible',fig_visibility);
        tiledlayout(1,2,'TileSpacing','tight');
        FRpl_subtable = Rpl_subtable(Rpl_subtable.C>0,:);
        RRpl_subtable = Rpl_subtable(Rpl_subtable.C<0,:);
        nexttile;
        for i=1:size(FRpl_subtable,1)
            scatter(FRpl_subtable.pw_distance{i,1},FRpl_subtable.pw_timing{i,1},15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   hold on;
        end
        ylabel('dt during F replay (s)'); xlabel('Place field distance (m)');
        nexttile;
        for i=1:size(RRpl_subtable,1)
            scatter(RRpl_subtable.pw_distance{i,1},RRpl_subtable.pw_timing{i,1},15,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   hold on;
        end
        ylabel('dt during R replay (s)'); xlabel('Place field distance (m)');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot a few features for all events, good and shuffled (when available)
        figure('units','normalized','outerposition',[.2 .3 .5 .35],'Visible',fig_visibility);
        tiledlayout(1,5,'TileSpacing','tight');
        cRT_else = setdiff(1:size(Rpl_table_tmp,1),cRT_good);
        maxD = max([Rpl_table_tmp.d1;Rpl_table_tmp.d2]);
        C_edges = [-1:0.2:1];
        D_edges =[0,10.^[-.5:0.1:log10(maxD)]];
        T_edges = [0,10.^[0:0.1:log10(3600)]];
        nexttile;
        histogram(Rpl_table_tmp.C(cRT_good),C_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table_tmp.C(cRT_else),C_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');
        histogram(horzcat(Rpl_table_tmp.shuffles{:}),C_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','k');   hold off;
        ylabel('Fraction'); xlabel('Rank Correlation'); legend('Good','Excluded','Shuffled','Location','northoutside');
        nexttile;
        histogram(Rpl_table_tmp.d1(cRT_good),D_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table_tmp.d1(cRT_else),D_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Distance from Takeoff (m)');
        nexttile;
        histogram(Rpl_table_tmp.d2(cRT_good),D_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table_tmp.d2(cRT_else),D_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Distance from Landing (m)');
        nexttile;
        histogram(Rpl_table_tmp.t2tko(cRT_good),T_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table_tmp.t2tko(cRT_else),T_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Time from Nrst-Takeoff (s)');
        nexttile;
        histogram(Rpl_table_tmp.t2lnd(cRT_good),T_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table_tmp.t2lnd(cRT_else),T_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Time from Nrst-Landing (s)');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot replays
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(5,20,'TileSpacing','tight');
        for zz=cRT_subset'
            tL = cRt(zz,1);
            tR = cRt(zz,2);
            nexttile;  
            for nc= 1:n_place_cells
                %plotSpikes_AF_v0(s{sorted_tko_plcCells(nc)}, nc, 'r', 1);
                plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
            end
            ylim([0 n_place_cells+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[num2str(candidate_rkcorr_c(zz),1)],[num2str((tR-tL)*1000,3), ' ms']});
            title(['p ', num2str(candidate_rkcorr_p(zz),1)]);
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot replays with LFP
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'Visible',fig_visibility);
        tiledlayout(15,20,'TileSpacing','tight');
        i=1;
        for zz=cRT_subset'
            tL = cRt(zz,1);
            tR = cRt(zz,2);
            nexttile(40*floor((i-1)/20)+i,[2 1]);
            for nc= 1:n_place_cells
                plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
            end
            ylim([0 n_place_cells+1]); xlim([tL, tR]); yticks([]);
            xticks([tL, tR]);   xticklabels({[num2str(candidate_rkcorr_c(zz),1)],[num2str(candidate_rkcorr_p(zz),1)]});
            nexttile(40*(1+floor((i-1)/20))+i);   plot(RPL_out.time_vector,LFP_rp,'b');    xlim([tL, tR]);   yticks([]);
            xticks([tL, tR]);   xticklabels({'',[num2str((tR-tL)*1000,3), ' ms']});
            i=i+1;
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot cross correlation of replay and SWR
        figure('units','normalized','outerposition',[.05 .1 .2 .3],'Visible',fig_visibility);
        tiledlayout(1,2,'TileSpacing','tight');
        [RS_cross_corr,RS_bin_centers] = cross_correlogram_AF_v0(Rpl_subtable.tL+Rpl_subtable.dt/2,RPL_out.table.t,0.5,0.025);
        [SS_cross_corr,RS_bin_centers] = cross_correlogram_AF_v0(RPL_out.table.t,RPL_out.table.t,0.5,0.025);
        nexttile; bar(RS_bin_centers,RS_cross_corr,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5); xlabel('Time(s)');  ylabel('Fraction'); title('Replay-SWR AC');
        nexttile; bar(RS_bin_centers,SS_cross_corr,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5); xlabel('Time(s)');  ylabel('Fraction'); title('SWR-SWR AC');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot when Replays are happening
        figure('units','normalized','outerposition',[0 .4 1 .4],'Visible',fig_visibility);
        tiledlayout(2,1,'TileSpacing','none');
        ax = [];    ax(1) = nexttile;   hold on;
        area(t,bflying,0,'FaceColor',[0 0 0],'FaceAlpha',0.2,'LineStyle','none');
        plot((t(f_clus.strt_frame(f_clus.id == j)')*[1 1])',(ones(size(f_clus.strt_frame(f_clus.id == j)'))*[0 1])','Color',col_clus(j,:));
        plot((cRt(cRT_subset,1)*[1 1])',(ones(size(cRt(cRT_subset,1)))*[0 -.5])','Color','r');    set(gca,'TickLength',[0 0]);  box('off'); ylabel('Replays/Flights');  yticks([]);
        ax(2) = nexttile; plot(t_Hi,LFP_mn);    ylabel('LFP');
        xlabel('Time (s)'); linkaxes(ax,'x');   xlim('tight');  box('off'); yticks([]);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot replays with Position of the bat
        figure('units','normalized','outerposition',[.05 .1 .9 .8],'Visible',fig_visibility);
        tiledlayout(10,20,'TileSpacing','tight');
        i=1;
        for zz=cRT_subset'
            tL = cRt(zz,1);
            tR = cRt(zz,2);
            sL = Rpl_table_tmp.smpL(zz);
            sR = Rpl_table_tmp.smpR(zz);
            
            nexttile(20*floor((i-1)/20)+i);
            for nc= 1:n_place_cells
                plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
            end
            ylim([0 n_place_cells+1]); xlim([tL, tR]); yticks([]);
            xticks([tL, tR]);   xticklabels({'',[num2str((tR-tL)*1000,3), ' ms']});
            nexttile(20*(1+floor((i-1)/20))+i);
            plot3(mean_path3D(:,1),mean_path3D(:,2),mean_path3D(:,3),'-','LineWidth',1,'Color', col_clus(j,:));
            xlim(r_lim(1,:)); ylim(r_lim(2,:));   zlim(r_lim(3,:));  view(0,90);    hold on;
            textscatter(avg_takeoff(1),avg_takeoff(2),">>");
            scatter3(r(sL:sR,1),r(sL:sR,2),r(sL:sR,3),50,'MarkerFaceColor', 'r','MarkerEdgeColor','none');
            hold off;   axis equal;  axis off;
            i=i+1;
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %==================================================== DECODING =====================================================================
        %=== For decoding
        %NP_table.use4decoding = NP_table.peakHz>0;
        NP_table.use4decoding = NP_table.peakHz>0 & NP_table.fr<5 & NP_table.corr_m>0.5 & NP_table.spkPerflight>0;
        NP_subtable_dec = NP_table(NP_table.use4decoding,:);
        n_dec_cells = size(NP_subtable_dec,1);
        
        %=== DECODING Initialization
        n_bins = floor((t(T)-t_bin_dur)/(t_bin_dur-t_bin_ovl))+1;   % Number of time bins
        st_times = (0:n_bins-1)*(t_bin_dur - t_bin_ovl);            % Start times
        ed_times = st_times + t_bin_dur;                            % End times
        t_d = zeros(1, 2 * n_bins);                                 % All times (initialize)
        t_d(1:2:end) = st_times;    t_d(2:2:end) = ed_times;        % All times (define)
        n_vector = zeros(1,n_dec_cells);                            % Initialize vector containing the number of spikes fired by each place cell
        p_x = NP_unitOnClus{1,j}(1).prob_x;                         % Probability of being on a given bin of the flight cluster
        f_x = zeros(numel(p_x),n_dec_cells);                        % Vector with the firing rates of place cells in a given bin of the flight cluster
        for nc=1:n_dec_cells, f_x(:,nc) = NP_subtable_dec.plc_map{nc,1};end
        
        %=== Generate mirrored spikes for smoothing (2 on the left and 2 on the right of each spike)
        s_smooth_1 = cellfun(@(x) x+0.005,s(:,1),'Uniformoutput',false);
        s_smooth_2 = cellfun(@(x) x+0.010,s(:,1),'Uniformoutput',false);
        s_smooth_3 = cellfun(@(x) x-0.005,s(:,1),'Uniformoutput',false);
        s_smooth_4 = cellfun(@(x) x-0.010,s(:,1),'Uniformoutput',false);
        s_smooth = cellfun(@(x,y,z,t,u) sort([x;y;z;t;u]),s(:,1),s_smooth_1,s_smooth_2,s_smooth_3,s_smooth_4,'Uniformoutput',false);
        
        %=== DECODING Plotting
        figure('units','normalized','outerposition',[.2 .3 .6 .5],'Visible',fig_visibility);
        tiledlayout(5,20,'TileSpacing','compact');
        for zz=cRT_subset'
            
            %=== Find closest bins for decoding start/stop (NEW VERSION, FAST)
            st_times_red_vec = st_times>cRt(zz,1)-.5 & st_times<cRt(zz,2)+.5;   % Restrict search to 0.5s around the times of interest
            st_times_red = st_times(st_times_red_vec);
            st_times_red_start = find(st_times_red_vec);
            st_times_red_start = st_times_red_start(1);
            i_bin_strt = st_times_red_start+knnsearch(st_times_red',cRt(zz,1))-1;
            i_bin_stop = st_times_red_start+knnsearch(st_times_red',cRt(zz,2))-1;
            
            p_dec_global = [];                              % Initialize
            %=== Loop across bins
            for i_dec = i_bin_strt:i_bin_stop
                %=== Count spikes of each unit on the time bin
                for nc=1:n_dec_cells,n_vector(nc) = histcounts(s_smooth{NP_subtable_dec.cell_id(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
                %=== Decoding
                p_x = double(p_x>-1);   % Uniform Prior
                p_dec_global = [p_dec_global, decode_1Dpos_AF_v1(n_vector,t_bin_dur,p_x,f_x)];
            end
            
            PPlim = prctile(p_dec_global,[1 98],'all')';    % Deal with imagesc limits
            if ~any(PPlim),PPlim = [0 1];end
            nexttile;   imagesc(imgaussfilt(p_dec_global,smooth_f,'FilterDomain','spatial'),PPlim);
            colormap('hot');  ylabel('Spatial bin'); axis off;  set(gca,'YDir','normal');
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Decoding all the events and storing the posterior (no smoothing)
        for zz=1:n_candidate_events
            i_bin_strt = knnsearch_fast_AF_v0(st_times',cRt(zz,1),0.5);
            i_bin_stop = knnsearch_fast_AF_v0(st_times',cRt(zz,2),0.5);
            p_dec_global = zeros(size(p_x,1),i_bin_stop-i_bin_strt+1);
            counter=1;
            for i_dec = i_bin_strt:i_bin_stop
                for nc=1:n_dec_cells,n_vector(nc) = histcounts(s_smooth{NP_subtable_dec.cell_id(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
                p_x = double(p_x>-1);   % Uniform Prior
                p_dec_global(:,counter) = decode_1Dpos_AF_v1(n_vector,t_bin_dur,p_x,f_x);
                counter=counter+1;
            end
%             tmp_sum = sum(p_dec_global);
%             if any(tmp_sum(tmp_sum>0)<0.99)
%             zz
%             end
            Rpl_table_tmp.posterior{zz,1} = p_dec_global;
        end
        
        %=== Add to the Rpl table of all clusters
        Rpl_table = [Rpl_table;Rpl_table_tmp];
    end
end

%% LOOK AT FLIGHTS, REPLAYS, SWR, ACCELEROMETER and LFP power
for hide=1
    
    %=== Define time bins
    t_bin = [t(1):10:t(end)];
    t_ctr = t_bin(1:end-1)+mean(diff(t_bin))/2;
    
    %=== Count Replays and Ripples
    %[Nr_pb] = histcounts(Rpl_table.tL,t_bin);
    if ~isempty(Rpl_table)
        Nr_pb = histcounts(Rpl_table.tL(Rpl_table.n>4 & Rpl_table.p<0.05),t_bin);
    else
        Nr_pb = zeros(1,numel(t_ctr));
    end
    Ns_pb = histcounts(RPL_out.table.t,t_bin);
    
    figure('units','normalized','outerposition',[0 .3 1 .4],'Visible',fig_visibility);
    %tiledlayout(13,1,'TileSpacing','none');
    area(t,bflying,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');   hold on;
    area(t_ctr,normalize(Nr_pb,'range'),'FaceAlpha',0.5);
    area(t_ctr,-normalize(Ns_pb,'range'),'FaceAlpha',0.5);
    plot(t,-normalize(a_mdg,'range'));
    xlim('tight');    legend('Flights','Replays','Ripples','Acc.');
    box('off'); set(gca,'TickLength',[0 0]);
    xlabel('Time'); yticks([]); ylabel('Rate (a.u.)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% SAVE THE DATA
for hide=1
    if options.savedata
        save([analysis_directory,'/Analyzed_NPs_', unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},'.mat'],'Rpl_table','NP_unit','MSC');
    end
end

