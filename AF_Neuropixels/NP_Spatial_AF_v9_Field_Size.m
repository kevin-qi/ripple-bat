function NP_Spatial_AF_v9_Field_Size(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

fig_visibility = 'off';

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    mrgFile = dir('merges_*');              load([mrgFile.folder '/' mrgFile.name]);    % This contains the ids of potential merges
    nrnFile = dir('SU_kilosort*');                                                      % Single Units
    unique_ID = options.unique_ID;                                                      % Session identifier
    load('LFP_probe2.mat');                                                             % This file contains LFP data from probe 1
    load('RPL_probe2.mat');                                                             % This file contains Ripple data from probe 1
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

%% SHOW EXAMPLE PLACE CELLS 
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
        
        %=== Create the cluster table
        NP_clus = struct();
        NP_clus.unique_ID = unique_ID;
        NP_clus.clus_id = j;
        NP_clus.n_cells = n_cells;
        NP_clus.length = fc_len(end);
        NP_clus.table = NP_table;
        
        if j==2
            NP_cluster = NP_clus;
        else
            NP_cluster = [NP_cluster; NP_clus];
        end
        
    end
end

%% SAVE THE DATA
for hide=1
    if options.savedata
        save([analysis_directory,'/Analyzed_NPs_', unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},'.mat'],'NP_cluster');
    end
end

