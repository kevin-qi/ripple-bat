function NP_Spatial_Theta_AF_v1_Ray(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

rng(1); % For repeatibility
fig_visibility = 'off';

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    mrgFile = dir('merges_*');              load([mrgFile.folder '/' mrgFile.name]);    % This contains the ids of potential merges
    nrnFile = dir('SU_kilosort*');                                                      % Single Units
    bpcFile = dir('BPC_Dataset*');          load([bpcFile.folder '/' bpcFile.name]);    % Best probe and channel
    unique_ID = options.unique_ID;                                                      % Session identifier
    load(['LFP_probe',num2str(best_prb_ch(1)),'.mat']);                                 % This file contains LFP data from best probe (as indicated in the dataset table)
    load(['RPL_probe',num2str(best_prb_ch(1)),'.mat']);                                 % This file contains Ripple data from best probe (as indicated in the dataset table)
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
    NP_unit.unique_ID = repmat(unique_ID,size(NP_unit,1),1);
    
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
    
    %=== Recluster Flights (default 0.8)
    alpha_clus = .6;
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
    NP_unit_FS = NP_unit(NP_unit.fr>2,:);                   % Store high firing cells (putative interneurons) in a different variable (def 5 Hz)
    NP_unit = NP_unit(NP_unit.fr<2,:);                      % Exclude neurons with high-firing rate
    n_cells = size(NP_unit,1);                              % Number of units
    n_cells_FS = size(NP_unit_FS,1);                        % Number of units (putative interneurons)
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 10;                                             % Number of repetitions for shuffling (Place Cells)
    n_rep_mod_idx = 5;                                      % Number of repetitions for shuffling (Theta Sweeps)
    times_to_shift = t(randi([3*Fs T-3*Fs],1,n_rep));       % Time intervals for circshift
    smpls_to_shift = randi([3*Fs T-3*Fs],1,n_rep);          % Sample intervals for circshift
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
    t_Hi = NP_imu.t;                                        % Time for fast sampling (500 Hz, IMU and LFP)
    bflying_Hi = logical(interp1(t',double(bflying),t_Hi,'nearest',0))'; % Flight vector at higher sampling
    Fs_Hi = NP_imu.Fs;                                      % Sampling frequency for fast sampling (500 Hz, IMU and LFP)
    opt_int = best_prb_ch(2)+[-2:2];                        % Optimal channel interval for getting the LFP
    min_prctile_power = 25;                                 % Min Percentile for power of a good oscillation (tmr)
    min_prctile_power_wbt = 0;                             % Min Percentile for power of a good oscillation (wbt)
    n_bins_phase = 10;
    bin_TMI = 0.01;
    max_lag_TMI = 0.5;
    Tht_rng = [4 11];                                       % Theta range
    Swr_rng = [100 200];                                    % Sharp-wave-ripple range
    Dlt_rng = [1 4];                                        % Delta range
    Lwp_rng = [1 10];                                       % Low Pass Range (from Eliav et al., 2018)
    Sgm_rng = [30 60];                                      % Slow gamma
    SWP_table = [];
    n_rep_phase = 100;                                      % Shuffles for circular tests of phase locking
    
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
    
    %=== Calculate Population Spike Density
    all_s = sort(vertcat(s{:,1}));
    all_rate = kernel_rate_AF_v1(all_s,0.1,t)/n_cells;
    
    %=== Calculate Population Spike Density for putative interneurons
    all_s_FS = sort(vertcat(s_FS{:,1}));
    all_rate_FS = kernel_rate_AF_v1(all_s_FS,0.1,t)/n_cells_FS;
    
    %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                                     % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_Hi);                              % Filtered at the wing-beat frequency
    [up,lo] = envelope(a_flt_NP,round(0.06*Fs_Hi),'peak');                  % Upper and lower envelopes
    env = normalize(up - lo,'range');                                       % Amplitude of the envelope
    env_th = otsuthresh(histcounts(env));                                   % Threshold (based on Otsu method). Can be set at 0.35
    wBeats = movsum(env>env_th,2*Fs_Hi)>Fs_Hi/5;                            % Euristic criterion for flight detection
    a_mdg = interp1(t_Hi,movmean(abs(a_abs_NP-1),Fs_Hi*0.5),t,'linear',0);  % mean deviation from g in a 0.5s integration window
    immobile = a_mdg<0.04;                                                  % Immobility vector
    
    %=== Calculate characteristic frequency of wingbeat signal
    n = 2^nextpow2(numel(a_flt_NP(bflying)));           % Optimal number of points for FFT
    Y = fft(a_flt_NP(bflying),n);                       % Calculate fft
    P = abs(Y/n).^2;                                    % Power density at all frequences (positive and negative)
    f = Fs_Hi*(0:n/2)/n;                                % Fs is the sampling frequency
    PSD = P(1:n/2+1);                                   % Power spectral density
    sm_PSD = smoothdata(PSD,'movmedian',n/Fs_Hi*3);     % Smoothed
    [~,max_loc] = max(sm_PSD); f_wBeats = f(max_loc);   % Characteristic frequency of wingbeat signal
    
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
    
    %=== Create click trace, useful later
    click_trace = zeros(numel(a_abs_NP),1);
    for i=1:numel(Detected_Clicks.times)
        click_trace(knnsearch_fast_AF_v0(NP_imu.t',Detected_Clicks.times(i),0.1)) = 1;
    end
    
    %=== Transform the s cell into a cell with the t_Hi-samples of each spike
    valid_s = find(all_s>0 & all_s<t_Hi(end));
    st1 = all_s(valid_s(1)  );
    st2 = all_s(valid_s(end));
    ss1 = knnsearch_fast_AF_v0(t_Hi,st1,1);
    ss2 = knnsearch_fast_AF_v0(t_Hi,st2,1);
    m_slope = (ss2-ss1)/(st2-st1);
    q_intcp = ss1-m_slope*st1;
    s_smp_hi = s(:,1);           % Single units' spikes samples
    s_smp_flight_hi = s(:,1);    % Single units' spikes samples (flight)
    for nc = 1:n_cells
        s_smp_hi{nc,1} = round(s{nc,1}*m_slope+q_intcp);
        [~,s_flight,~] = count_spikes_AF_v3(s{nc,1},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
        s_flight = s_flight(s_flight>0 & s_flight<t_Hi(end));
        s_smp_flight_hi{nc,1} = round(s_flight*m_slope+q_intcp);
    end
    
    %=== Transform the s cell into a cell with the t-samples of each spike
    valid_s = find(all_s>0 & all_s<t(T));
    st1 = all_s(valid_s(1)  );
    st2 = all_s(valid_s(end));
    ss1 = knnsearch_fast_AF_v0(t,st1,1);
    ss2 = knnsearch_fast_AF_v0(t,st2,1);
    m_slope = (ss2-ss1)/(st2-st1);
    q_intcp = ss1-m_slope*st1;
    s_smp_lo = s(:,1);           % Single units' spikes samples
    s_smp_flight_lo = s(:,1);   % Single units' spikes samples (flight)
    for nc = 1:n_cells
        s_smp_lo{nc,1} = round(s{nc,1}*m_slope+q_intcp);
        [~,s_flight,~] = count_spikes_AF_v3(s{nc,1},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
        s_flight = s_flight(s_flight>0 & s_flight<t(end));
        s_smp_flight_lo{nc,1} = round(s_flight*m_slope+q_intcp);
    end
    
    %=== Get the transformation from s_smp_lo to s_smp_hi (useful later)
    s_smp_lo_all = vertcat(s_smp_lo{:});
    s_smp_hi_all = vertcat(s_smp_hi{:});
    p_lohi = polyfit(s_smp_hi_all, s_smp_lo_all, 1);
    m_slope = p_lohi(1);
    q_intcp = p_lohi(2);
    
end

%% LFP PROCESSING
for hide=1
    
    %=== Process the LFP (transform into double and interpolate at 500 Hz sampling)
    LFP = mean(interp1(red_out.t_ds,double(red_out.lfp(:,opt_int)).*red_out.voltage_scaling,t_Hi,'linear',0),2);     % Interpolate LFP at the accelerometer time samples
    
    %=== Extract different frequency bands
    LFP_tht = bandpass(LFP,Tht_rng,Fs_Hi);       % Filtered in the theta band
    LFP_dlt = bandpass(LFP,Dlt_rng,Fs_Hi);       % Filtered in the delta band
    LFP_lwp = bandpass(LFP,Lwp_rng,Fs_Hi);       % Filtered in the low pass band (as Eliav et al., 2018)
    LFP_sgm = bandpass(LFP,Sgm_rng,Fs_Hi);       % Filtered in the slow gamma band
    
    %=== Extract relevant phases and amplitudes
    wbt_phase_tmp = wrapTo2Pi(angle(hilbert(a_flt_NP)))-pi;   wbt_power_tmp = abs(hilbert(a_flt_NP)).^2;    wbt_freq_tmp = instfreq(a_flt_NP,Fs_Hi,'Method','hilbert');
    tht_phase = wrapTo2Pi(angle(hilbert(LFP_tht)))-pi;        tht_power = abs(hilbert(LFP_tht)).^2;         tht_freq = instfreq(LFP_tht,Fs_Hi,'Method','hilbert');
    dlt_phase = wrapTo2Pi(angle(hilbert(LFP_dlt)))-pi;        dlt_power = abs(hilbert(LFP_dlt)).^2;         dlt_freq = instfreq(LFP_dlt,Fs_Hi,'Method','hilbert');
    lwp_phase = wrapTo2Pi(angle(hilbert(LFP_lwp)))-pi;        lwp_power = abs(hilbert(LFP_lwp)).^2;         lwp_freq = instfreq(LFP_lwp,Fs_Hi,'Method','hilbert');
    sgm_phase = wrapTo2Pi(angle(hilbert(LFP_sgm)))-pi;        sgm_power = abs(hilbert(LFP_sgm)).^2;         sgm_freq = instfreq(LFP_sgm,Fs_Hi,'Method','hilbert');
    
    %=== Invalidate wingbeat outside of flight
    wbt_phase = zeros(size(wbt_phase_tmp));
    wbt_power = zeros(size(wbt_phase_tmp));
    wbt_freq = NaN(size(wbt_phase_tmp));
    for i=1:f_num_Hi
        wbt_phase(f_smp_Hi(i,1):f_smp_Hi(i,2))= wbt_phase_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
        wbt_power(f_smp_Hi(i,1):f_smp_Hi(i,2))= wbt_power_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
        wbt_freq(f_smp_Hi(i,1):f_smp_Hi(i,2))=  wbt_freq_tmp(f_smp_Hi(i,1):f_smp_Hi(i,2));
    end
    
    %=== Extract phase and power of Eliav's style LFP
    [~,min_locs] = findpeaks(-LFP_lwp,'MinPeakDistance',round(Fs_Hi*0.01)); % Find minima
    pi_tmp_phase = 2*pi*ones(numel(min_locs),1);                            % Assign phase 0 to minima
    pi_tmp_phase = cumsum(pi_tmp_phase);
    tmr_phase = interp1(min_locs,pi_tmp_phase,1:numel(t_Hi),'linear',0)';   % Interpolate
    tmr_phase = wrapTo2Pi(tmr_phase)-pi;                                    % Wrap to -pi:pi
    tmr_phase = wrapTo2Pi(tmr_phase)-pi;                                    % Shift, such that minima have 0 phase
    %plot(tmr_phase);    hold on;    plot(normalize(LFP_lwp));  plot(normalize(LFP));  
    %plot(wbt_phase);    hold on;    plot(normalize(a_flt_NP));
    tmr_power = lwp_power;
    for i=1:numel(min_locs)-1
        tmr_power(min_locs(i):min_locs(i+1)-1)  = mean(lwp_power(min_locs(i):min_locs(i+1)-1),'omitnan');
    end
    
    %=== Find theta bouts (fast way)
    avg_tht_power = smoothdata(tht_power,'gaussian',round(Fs_Hi*2));
    avg_dlt_power = smoothdata(dlt_power,'gaussian',round(Fs_Hi*2));
    Theta2Delta = avg_tht_power./avg_dlt_power;
    [btheta,thtBt_num,tht_bt_str,tht_bt_stp,~] = BiLevel_Segm_AF_v0(Theta2Delta,3,Fs_Hi,0.1,1);
    thtBt_smp = [tht_bt_str,tht_bt_stp];
    
    %=== Accumulate delta power time course around theta bouts
    interval = [round(-5*Fs_Hi):round(5*Fs_Hi)-1];
    dlt_around_bout = [];    dlt_around_ctrl = [];
    tht_around_bout = [];    tht_around_ctrl = [];
    for i=1:numel(tht_bt_str)
        if tht_bt_str(i)+round(-5*Fs_Hi)>0 && tht_bt_str(i)+round(5*Fs_Hi)<numel(t_Hi)
            dlt_around_bout = [dlt_around_bout, avg_dlt_power(tht_bt_str(i)+interval)];
            dlt_around_ctrl = [dlt_around_ctrl, avg_dlt_power(-interval(1) + randi(numel(t_Hi)-2*interval(end))+interval)];
            tht_around_bout = [tht_around_bout, avg_tht_power(tht_bt_str(i)+interval)];
            tht_around_ctrl = [tht_around_ctrl, avg_tht_power(-interval(1) + randi(numel(t_Hi)-2*interval(end))+interval)];
        end
    end
    tme_around_bout = interval/Fs_Hi;
    
    %=== LFP and PSD
    LFP_r = LFP(~bflying_Hi);      [PSD_r,f_PSD_r] = pwelch(LFP_r,[],[],[],Fs_Hi);  PSD_r = PSD_r./sum(PSD_r,1);  f_PSD_smpl_r = mean(diff(f_PSD_r)); % Rest
    LFP_f = LFP( bflying_Hi);      [PSD_f,f_PSD_f] = pwelch(LFP_f,[],[],[],Fs_Hi);  PSD_f = PSD_f./sum(PSD_f,1);  f_PSD_smpl_f = mean(diff(f_PSD_f)); % Flight
    LFP_b = LFP( btheta);      [PSD_b,f_PSD_b] = pwelch(LFP_b,[],[],[],Fs_Hi);  PSD_b = PSD_b./sum(PSD_b,1);  f_PSD_smpl_b = mean(diff(f_PSD_b)); % Flight
    
    %=== PSD of the accelerometer signal during flight
    [PSD_a,f_PSD_a] = pwelch(a_abs_NP(bflying_Hi),[],[],[],Fs_Hi);  PSD_a = PSD_a./sum(PSD_a,1);  f_PSD_smpl_a = mean(diff(f_PSD_a));
    
    %=== Average, smooth and normalize PSDs
    scaled_PSD_r = normalize(smoothdata(mean(PSD_r,2),'movmedian',1/f_PSD_smpl_r),'range');
    scaled_PSD_f = normalize(smoothdata(mean(PSD_f,2),'movmedian',1/f_PSD_smpl_f),'range');
    scaled_PSD_b = normalize(smoothdata(mean(PSD_b,2),'movmedian',1/f_PSD_smpl_b),'range');
    scaled_PSD_a = normalize(smoothdata(mean(PSD_a,2),'movmedian',1/f_PSD_smpl_a),'range');
    
    %=== Define pre-flight, flight and postflight time windows (using same duration)
    pre_flgt = false(size(wBeats)); for i=1:f_num_Hi,pre_flgt(f_smp_Hi(i,1)+[round(-2*Fs_Hi) :        0      ]) = 1;    end
    dur_flgt = false(size(wBeats)); for i=1:f_num_Hi,dur_flgt(f_smp_Hi(i,1)+[      0         : round(2*Fs_Hi)]) = 1;    end
    pst_flgt = false(size(wBeats)); for i=1:f_num_Hi,pst_flgt(f_smp_Hi(i,2)+[      0         : round(2*Fs_Hi)]) = 1;    end
    LFP_pre = LFP(pre_flgt);      [PSD_pre,f_PSD_pre] = pwelch(LFP_pre,[],[],[],Fs_Hi);  PSD_pre = PSD_pre./sum(PSD_pre,1);  f_PSD_smpl_pre = mean(diff(f_PSD_pre)); % Rest
    LFP_dur = LFP(dur_flgt);      [PSD_dur,f_PSD_dur] = pwelch(LFP_dur,[],[],[],Fs_Hi);  PSD_dur = PSD_dur./sum(PSD_dur,1);  f_PSD_smpl_dur = mean(diff(f_PSD_dur)); % Flight
    LFP_pst = LFP(pst_flgt);      [PSD_pst,f_PSD_pst] = pwelch(LFP_pst,[],[],[],Fs_Hi);  PSD_pst = PSD_pst./sum(PSD_pst,1);  f_PSD_smpl_pst = mean(diff(f_PSD_pst)); % Flight
    scaled_PSD_pre = normalize(smoothdata(mean(PSD_pre,2),'movmedian',1/f_PSD_smpl_pre),'range');
    scaled_PSD_dur = normalize(smoothdata(mean(PSD_dur,2),'movmedian',1/f_PSD_smpl_dur),'range');
    scaled_PSD_pst = normalize(smoothdata(mean(PSD_pst,2),'movmedian',1/f_PSD_smpl_pst),'range');
    
    %=== Quantify difference in theta power distribution during - before
    p_val_tht_DurPre = zeros(f_num_Hi,1);
    delta_tht_DurPre = zeros(f_num_Hi,1);
    theta2delta_flgt = zeros(f_num_Hi,1);
    for i=1:f_num_Hi
        tht_pwr_pre = tht_power(f_smp_Hi(i,1)+[round(-3*Fs_Hi):0]);
        tht_pwr_dur = tht_power(f_smp_Hi(i,1)+[0: round(3*Fs_Hi)]);
        p_val_tht_DurPre(i) = ranksum(tht_pwr_pre, tht_pwr_dur);
        delta_tht_DurPre(i) = (mean(tht_pwr_dur)-mean(tht_pwr_pre))/mean(tht_pwr_pre);
        theta2delta_flgt(i) = median(tht_power(f_smp_Hi(i,1):f_smp_Hi(i,2)))/median(dlt_power(f_smp_Hi(i,1):f_smp_Hi(i,2)));
    end
    
    %=== Store the (not normalized) PSDs in a variable for saving
    LFP_PSD = struct();
    LFP_PSD.unique_ID = unique_ID;
    LFP_PSD.PSD_pre = PSD_pre.*sum(PSD_pre,1);  LFP_PSD.f_PSD_pre = f_PSD_pre;
    LFP_PSD.PSD_dur = PSD_dur.*sum(PSD_dur,1);  LFP_PSD.f_PSD_dur = f_PSD_dur;
    LFP_PSD.PSD_pst = PSD_pst.*sum(PSD_pst,1);  LFP_PSD.f_PSD_pst = f_PSD_pst;
    LFP_PSD.PSD_r   = PSD_r.*sum(PSD_r,1);      LFP_PSD.f_PSD_r   = f_PSD_r;
    LFP_PSD.PSD_f   = PSD_f.*sum(PSD_f,1);      LFP_PSD.f_PSD_f   = f_PSD_f;
    LFP_PSD.PSD_b   = PSD_b.*sum(PSD_b,1);      LFP_PSD.f_PSD_b   = f_PSD_b;
    LFP_PSD.PSD_a   = PSD_a.*sum(PSD_a,1);      LFP_PSD.f_PSD_a   = f_PSD_a;
    
    %=== Store additional features of the LFP
    LFP_PSD.perc_thtpower_b = sum(scaled_PSD_b(f_PSD_b>Tht_rng(1) & f_PSD_b<Tht_rng(2)))/sum(scaled_PSD_b); % Percent power in the theta band during bouts
    LFP_PSD.perc_thtpower_f = sum(scaled_PSD_f(f_PSD_f>Tht_rng(1) & f_PSD_f<Tht_rng(2)))/sum(scaled_PSD_f); % Percent power in the theta band during flight
    LFP_PSD.f_num = f_num_Hi;
    LFP_PSD.delta_tht_DurPre = delta_tht_DurPre;
    LFP_PSD.p_val_tht_DurPre = p_val_tht_DurPre;
    LFP_PSD.theta2delta_flgt = theta2delta_flgt;
    LFP_PSD.fract_bouts_flight = sum(bflying_Hi(thtBt_smp(:,1)+diff(thtBt_smp,1,2)))./numel(bflying_Hi(thtBt_smp(:,1)+diff(thtBt_smp,1,2)));
    LFP_PSD.tht_bouts_dur = diff(thtBt_smp')'./Fs_Hi;
    LFP_PSD.dlt_around_bout = dlt_around_bout;
    LFP_PSD.dlt_around_ctrl = dlt_around_ctrl;
    LFP_PSD.tht_around_bout = tht_around_bout;
    LFP_PSD.tht_around_ctrl = tht_around_ctrl;
    LFP_PSD.tme_around_bout = tme_around_bout;
    LFP_PSD.tmr_power = [median(tmr_power(bflying_Hi)),median(tmr_power(~bflying_Hi)),median(tmr_power(datasample(find(~bflying_Hi),sum(bflying_Hi),'Replace',false)))];
    LFP_PSD.tmr_power1 = [mean(tmr_power(bflying_Hi)),mean(tmr_power(~bflying_Hi)),mean(tmr_power(datasample(find(~bflying_Hi),sum(bflying_Hi),'Replace',false)))];
    LFP_PSD.avg_tht_power = avg_tht_power;
    LFP_PSD.avg_dlt_power = avg_dlt_power;
    LFP_PSD.bflying_Hi = bflying_Hi;
    LFP_PSD.Fs_Hi = Fs_Hi;
    
end

%% LFP PLOTTING
%while 0
for hide=1
    
    %=== Plot the PSDs
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.4],'Visible',fig_visibility);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   plot(f_PSD_r,scaled_PSD_r,'LineWidth',2);    hold on;     plot(f_PSD_f,scaled_PSD_f,'LineWidth',2); plot(f_PSD_b,scaled_PSD_b,'LineWidth',2);  plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);
    rectangle('Position',[Tht_rng(1) 0 diff(Tht_rng) 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Rest)','LFP (Flight)','LFP (Bouts)','Accelerometer (Flight)');
    nexttile;   plot(f_PSD_r,scaled_PSD_r,'LineWidth',2);    hold on;     plot(f_PSD_f,scaled_PSD_f,'LineWidth',2); plot(f_PSD_b,scaled_PSD_b,'LineWidth',2);  plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);
    [~,tmp_idx] = max(scaled_PSD_a);    plot(f_PSD_a(tmp_idx)*[1 1],ylim,'k--');
    xlim([3 12]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  set(gca, 'YScale', 'log');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot the PSDs (pre, dur, post flight)
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.4],'Visible',fig_visibility);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   plot(f_PSD_pre,scaled_PSD_pre,'LineWidth',2);    hold on;     plot(f_PSD_dur,scaled_PSD_dur,'LineWidth',2); plot(f_PSD_pst,scaled_PSD_pst,'LineWidth',2);
    rectangle('Position',[Tht_rng(1) 0 diff(Tht_rng) 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Pre)','LFP (Dur)','LFP (Post)');
    nexttile;   plot(f_PSD_pre,scaled_PSD_pre,'LineWidth',2);    hold on;     plot(f_PSD_dur,scaled_PSD_dur,'LineWidth',2); plot(f_PSD_pst,scaled_PSD_pst,'LineWidth',2);   plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);
    rectangle('Position',[Tht_rng(1) 0 diff(Tht_rng) 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Pre)','LFP (Dur)','LFP (Post)','Accelerometer');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot cross correlation between LFP and accelerometer during flight
    [cor_LFP2acc,lags] = xcorr(LFP_f,a_abs_NP(bflying_Hi),round(2*Fs_Hi));
    [coh, freq_coh] = mscohere(LFP_f,a_abs_NP(bflying_Hi), [], [], [], Fs_Hi);
    wb_maxima = abs(wbt_phase)>0.99*pi;
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.3],'Visible',fig_visibility);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot(lags./Fs_Hi,cor_LFP2acc);  hold on;    plot([0 0],ylim);   xlim('tight');  xlabel('Time lag'); ylabel('Correlation Acc-LFP');
    nexttile;   plot(freq_coh,smoothdata(coh,'movmedian',Fs_Hi*0.1));  xlabel('Frequency (Hz)');   ylabel('Coherence Acc-LFP');    xlim([0 20]);    hold on;    plot([f_wBeats f_wBeats],ylim,'k--');
    nexttile;   scatter(a_abs_NP(wb_maxima),-LFP(wb_maxima),6,'MarkerFaceColor', 'r','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
    xlim(prctile(a_abs_NP(wb_maxima),[2 98]));  ylim(prctile(-LFP(wb_maxima),[2 98]));  axis square;  lsline;
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Accumulate the LFP for each flight, aligned to takeoff and landing
    [~,tmp_idx] = sort(diff(f_smp_Hi,1,2));
    int_t = [-3 6]; smht = [0.5 1]; sat_c = [2 98];
    interval = [round(int_t(1)*Fs_Hi) : round(int_t(2)*Fs_Hi)];
    LFP_htmp = zeros(f_num_Hi,numel(interval),2);
    LTH_htmp = zeros(f_num_Hi,numel(interval),2);
    ACC_htmp = zeros(f_num_Hi,numel(interval),2);
    for i=1:f_num_Hi
        LFP_htmp(i,:,1) = LFP(f_smp_Hi(tmp_idx(i),1)+interval)';
        LTH_htmp(i,:,1) = LFP_tht(f_smp_Hi(tmp_idx(i),1)+interval)';
        ACC_htmp(i,:,1) = a_abs_NP(f_smp_Hi(tmp_idx(i),1)+interval)';
        LFP_htmp(i,:,2) = LFP(f_smp_Hi(tmp_idx(i),2)+interval)';
        LTH_htmp(i,:,2) = LFP_tht(f_smp_Hi(tmp_idx(i),2)+interval)';
        ACC_htmp(i,:,2) = a_abs_NP(f_smp_Hi(tmp_idx(i),2)+interval)';
    end
    
    %=== Plot the LFP for each flight, aligned to takeoff and landing
    figure('units','normalized','outerposition',[0.4 0.3 0.35 0.5],'Visible',fig_visibility);
    tiledlayout(2,3,'TileSpacing','tight');
    ax(1) = nexttile;   imagesc(int_t,[],imgaussfilt(ACC_htmp(:,:,1),smht),prctile(ACC_htmp,sat_c,'all')');  colormap(ax(1),gray);  hold on;
    plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
    ax(2) = nexttile;   imagesc(int_t,[],imgaussfilt(LFP_htmp(:,:,1),smht),prctile(LFP_htmp,sat_c,'all')');  colormap(ax(2),redblue); hold on;
    plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('LFP (All Freq.)');
    ax(3) = nexttile;   imagesc(int_t,[],imgaussfilt(LTH_htmp(:,:,1),smht),prctile(LTH_htmp,sat_c,'all')');  colormap(ax(3),redblue); hold on;      colorbar;
    plot([0 0],ylim,'k--');   xlabel('Time to takeoff (s)');  ylabel('Sorted Flight #'); title('LFP (Theta)');
    ax(4) = nexttile;   imagesc(int_t,[],imgaussfilt(ACC_htmp(:,:,2),smht),prctile(ACC_htmp,sat_c,'all')');  colormap(ax(4),gray);  hold on;
    plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #');
    ax(5) = nexttile;   imagesc(int_t,[],imgaussfilt(LFP_htmp(:,:,2),smht),prctile(LFP_htmp,sat_c,'all')');  colormap(ax(5),redblue); hold on;
    plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #');
    ax(6) = nexttile;   imagesc(int_t,[],imgaussfilt(LTH_htmp(:,:,2),smht),prctile(LTH_htmp,sat_c,'all')');  colormap(ax(6),redblue); hold on;
    plot([0 0],ylim,'k--');   xlabel('Time to landing (s)');  ylabel('Sorted Flight #');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');  colorbar;
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot each flight (70 max)
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(7,20,'TileSpacing','compact');
    int_t = [-2 4]; interval = [round(int_t(1)*Fs_Hi) : round(int_t(2)*Fs_Hi)];
    for i=1:min(f_num_Hi,70)
        zoom_int = f_smp_Hi(i,1)+interval;
        nexttile;   plot(t_Hi(zoom_int),normalize(LFP_tht(zoom_int)'),'r');
        hold on;    plot(t_Hi(zoom_int),normalize(a_abs_NP(zoom_int)')+5,'k');
        plot(t_Hi(zoom_int),2*click_trace(zoom_int)'-5,'b');
        plot(t_Hi(f_smp_Hi(i,1))*[1 1],ylim,'k--'); yticks([]);
        if i==min(f_num_Hi,70),xlabel('Time (s)');else,xticks([]); end
        if i==1,ylabel('Flight aligned LFP');end
        
        nexttile;
        [PS_LFP,freq_SG] = cwt(LFP(zoom_int),Fs_Hi);
        imagesc(int_t,[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[.5 100])); shading interp;  colormap(hot); hold on;
        plot([0 0],ylim,'k--');
        set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([0.5 20]);
        if i==1,yticks([1 5 10 20]); ylabel('Freq (Hz)');else, yticks([]); end
        
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot theta bouts (70 max)
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(7,20,'TileSpacing','compact');
    int_t = [-2 4]; interval = [round(int_t(1)*Fs_Hi) : round(int_t(2)*Fs_Hi)];
    for i=1:min(thtBt_num,70)
        zoom_int = thtBt_smp(i,1)+interval;
        nexttile;   plot(t_Hi(zoom_int),normalize(LFP_tht(zoom_int)'),'r');
        hold on;    plot(t_Hi(zoom_int),normalize(a_abs_NP(zoom_int)')+5,'k');
        plot(t_Hi(zoom_int),2*click_trace(zoom_int)'-5,'b');
        plot(t_Hi(thtBt_smp(i,1))*[1 1],ylim,'k--');
        plot(t_Hi(thtBt_smp(i,2))*[1 1],ylim,'k--');
        if i==min(thtBt_num,70),xlabel('Time (s)');else,xticks([]); end
        if i==1,ylabel('Theta Bouts');end
        
        nexttile;
        [PS_LFP,freq_SG] = cwt(LFP(zoom_int),Fs_Hi);
        imagesc(int_t,[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[.5 100])); shading interp;  colormap(hot); hold on;
        plot([0 0],ylim,'k--');
        set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([0.5 20]);    title(num2str(i));
        yticks([1 5 10 20]); ylabel('Freq (Hz)');
        %if i==1,yticks([1 5 10 20]); ylabel('Freq (Hz)');else, yticks([]); end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot delta and theta around theta bouts
    figure('units','normalized','outerposition',[.3 0.3 .25 .35],'Visible',fig_visibility);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    data_tmp = dlt_around_ctrl';
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');  hold on;
    data_tmp = dlt_around_bout';
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    xlabel('Time from bout onset (s)'); ylabel('Delta power');  plot([0 0],ylim,'k--');
    
    nexttile;
    data_tmp = tht_around_ctrl';
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');  hold on;
    data_tmp = tht_around_bout';
    plotWinterval_AF_v0(tme_around_bout,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    xlabel('Time from bout onset (s)'); ylabel('Theta power');  plot([0 0],ylim,'k--');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot the accelerometer and the LFP
    figure('units','normalized','outerposition',[0 0.3 1 0.5],'Visible',fig_visibility);
    tiledlayout(5,1,'TileSpacing','none');
    bx(1) = nexttile;       plot(t_Hi,a_flt_NP,'k');    xticks([]); set(gca,'TickLength',[0 0]);    ylabel('Accelerometer');
    bx(2) = nexttile;       plot(t_Hi,LFP,'k');         xticks([]); set(gca,'TickLength',[0 0]);    ylabel('LFP All');
    bx(3) = nexttile;       plot(t_Hi,LFP_tht,'r');     xticks([]); set(gca,'TickLength',[0 0]);    ylabel('LFP Theta');
    bx(4) = nexttile;       plot(t_Hi,LFP_dlt,'b');     xticks([]); set(gca,'TickLength',[0 0]);    ylabel('LFP Delta');
    bx(5) = nexttile;       plot(t_Hi,Theta2Delta,'k'); hold on;
    area(t_Hi,btheta*max(Theta2Delta),0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    box('off');                                 ylabel('Theta/Delta ratio');
    linkaxes(bx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot distribution of Tamir's phase during flight vs rest (equal number of samples)
    figure('units','normalized','outerposition',[.3 .2 .2 .4],'Visible',fig_visibility);
    pwr_bins = linspace(min(tmr_power),max(tmr_power),100);
    histogram(tmr_power(bflying_Hi),                                                   pwr_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r'); hold on;
    histogram(tmr_power(datasample(find(~bflying_Hi),sum(bflying_Hi),'Replace',false)),pwr_bins,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    legend('Flight','Rest');    xlabel('Power (Hilbert), Matching flight-rest duration');    ylabel('Fraction');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot session timeline
    figure('units','normalized','outerposition',[0 .3 1 .6],'Visible',fig_visibility);
    tiledlayout(13,1,'TileSpacing','none');
    cx(1) = nexttile(1,[4 1]);
    for nc= 1:n_cells
        plot(s{nc}, nc*ones(size(s{nc})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(5,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(6,[2 1]);   plot(t_Hi,LFP); hold on;    plot(t_Hi,LFP_lwp);  xticks([]); legend('All','Theta');  ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
   
    
end

%% PLACE CELLS ON FLIGHT CLUSTERS
%while 0
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
            NP_unit(nc).f_clus(j).spk2wbt_rho = NaN;
            NP_unit(nc).f_clus(j).spk2wbt_pvl = NaN;
            NP_unit(nc).f_clus(j).spk2wbt_dtp = NaN;
            NP_unit(nc).f_clus(j).spk2wbt_pvl_sh = NaN;
            
            %=== Calculate p_value SI
            non_nan = nnz(~isnan(SI_1(2:end)));
            if isnan(SI_1(1)), p_val_l = NaN;   p_val_r = NaN;
            else, p_val_r =  nnz(SI_1(2:end)>SI_1(1))/non_nan;p_val_l =  nnz(SI_1(2:end)<SI_1(1))/non_nan;end
            NP_unit(nc).f_clus(j).SI_value = SI_1(1,1);
            NP_unit(nc).f_clus(j).SI_p_val = p_val_r;
            
            %=== Add the firing rate at the peak and the location of the peak
            NP_unit(nc).f_clus(j).field = map_1(1,:);
            [NP_unit(nc).f_clus(j).peakHz,NP_unit(nc).f_clus(j).field_loc] = max(map_1(1,:));
            NP_unit(nc).f_clus(j).asymm_frct = sum(NP_unit(nc).f_clus(j).map(1,1:round(NP_unit(nc).f_clus(j).field_loc)))./sum(NP_unit(nc).f_clus(j).map(1,:));
            
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
%while 0
for hide=1
    
    %=== Autocorrelation and MI for principal cells
    for nc=1:n_cells
        
        %=== Initialize Autocorrelation, TMI and p value
        NP_unit(nc).AC = [];
        NP_unit(nc).AC_bins = [];
        %NP_unit(nc).AC_all = [];
        %NP_unit(nc).AC_all_bins = [];
        NP_unit(nc).AC_PSD = [];
        NP_unit(nc).AC_PSD_freq = [];
        
        %=== Get empirical values
        [~,s_flight,~] = count_spikes_AF_v3(s{nc,1},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
        [NP_unit(nc).AC,NP_unit(nc).AC_bins] = cross_correlogram_AF_v0(s_flight,s_flight,max_lag_TMI,bin_TMI);  % Autocorrelogram
        %[NP_unit(nc).AC_all,NP_unit(nc).AC_all_bins] = cross_correlogram_AF_v0(s{nc,1},s{nc,1},max_lag_TMI,bin_TMI);  % Autocorrelogram
        
        if ~isempty(s_flight)
            [NP_unit(nc).AC_PSD,NP_unit(nc).AC_PSD_freq] = Theta_mod_idx_AF_v1(s_flight,bin_TMI,max_lag_TMI);
            AC_PSD_temp = zeros(numel(NP_unit(nc).AC_PSD),n_rep_mod_idx);
            for jj=1:n_rep_mod_idx
                [~,~,s_flight_sh] = count_spikes_AF_v3(s{nc,1},t,[f_smp(:,1) f_smp(:,2)+1]);
                [AC_PSD_temp(:,jj),PSD_freq] = Theta_mod_idx_AF_v1(s_flight_sh,bin_TMI,max_lag_TMI);
            end
            
            AC_PSD_zscore = smoothdata((NP_unit(nc).AC_PSD-mean(AC_PSD_temp,2))./std(AC_PSD_temp,[],2),'gaussian',round(3/mean(diff(PSD_freq))));
            AC_PSD_pvalue = sum(NP_unit(nc).AC_PSD<AC_PSD_temp,2)./n_rep_mod_idx;
            AC_PSD_zscore(PSD_freq<Tht_rng(1) | PSD_freq>Tht_rng(2)) = -inf;                                                  % Look at the best value in the theta range
            [NP_unit(nc).mod_idx,tmp_idx] = max(AC_PSD_zscore);
            NP_unit(nc).AC_peak_freq = PSD_freq(tmp_idx);
            NP_unit(nc).AC_p_val = AC_PSD_pvalue(tmp_idx);
            NP_unit(nc).nspkpf = numel(s_flight)/f_num;
            
        else
            NP_unit(nc).mod_idx = NaN;
            NP_unit(nc).AC_peak_freq = NaN;
            NP_unit(nc).AC_p_val = NaN;
            NP_unit(nc).nspkpf = 0;
        end
        
    end
    
    %=== Plot the best cells
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.7],'Visible',fig_visibility);
    tiledlayout(4,10,'TileSpacing','tight');
    good_cells = [NP_unit.AC_p_val]<0.05 & ([NP_unit.AC_peak_freq]>Tht_rng(1) | [NP_unit.AC_peak_freq]<Tht_rng(2)) & [NP_unit.nspkpf]>1;
    good_cell_ids = find(good_cells);
    [~,tmp_idx] = sort([NP_unit(good_cells).mod_idx],'descend');
    for i=1:min(sum(good_cells),40)
        nexttile;
        jj = good_cell_ids(tmp_idx(i));
        area(NP_unit(jj).AC_bins,NP_unit(jj).AC,'FaceColor','k');   yticks([]); xlabel('Time lag (s)'); title(['MI: ',num2str(NP_unit(jj).mod_idx,2),...
            ', p = ',num2str(NP_unit(jj).AC_p_val,1),', f = ',num2str(NP_unit(jj).AC_peak_freq,2) ]);
        hold on; plot(repmat(1/f_wBeats*[-3:3],2,1),repmat(ylim,7,1)','r--');
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %     %=== Autocorrelation and TMI for interneurons
    %     for nc=1:n_cells_FS
    %
    %         %=== Initialize Autocorrelation, TMI and p value
    %         NP_unit_FS(nc).AC = [];
    %         NP_unit_FS(nc).AC_bins = [];
    %         NP_unit_FS(nc).TMI = NaN(1,n_rep_sweeps);
    %         NP_unit_FS(nc).p_val_TMI = NaN;
    %
    %
    %         %=== Loop across shuffles
    %         for jj=1:n_rep_sweeps+1
    %             [~,s_flight] = count_spikes_AF_v1(s_FS{nc,jj},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
    %             if ~isempty(s_flight)
    %                 if jj==1
    %                     [NP_unit_FS(nc).AC,NP_unit_FS(nc).AC_bins] = cross_correlogram_AF_v0(s_flight,s_flight,max_lag_TMI,bin_TMI);
    %                 end
    %                 NP_unit(nc).TMI(jj) = Theta_mod_idx_AF_v0(s_flight,bin_TMI,Tht_rng,max_lag_TMI);     % Theta modulation index
    %             end
    %         end
    %
    %         %=== Calculate p value
    %         if ~isnan(NP_unit_FS(nc).TMI(1)), NP_unit_FS(nc).p_val_TMI = sum(NP_unit_FS(nc).TMI>NP_unit_FS(nc).TMI(1))/sum(~isnan(NP_unit_FS(nc).TMI));end
    %
    %     end
    %
    %     %=== Plot the 40 best cells
    %     figure('units','normalized','outerposition',[0.1 0.2 0.8 0.7]);
    %     tiledlayout(4,10,'TileSpacing','tight');
    %     [~,tmp_idx] = sort([NP_unit_FS.p_val_TMI]);
    %     for i=1:min(n_cells_FS,40)
    %         nexttile;
    %         area(NP_unit_FS(tmp_idx(i)).AC_bins,NP_unit_FS(tmp_idx(i)).AC,'FaceColor','k');   yticks([]); xlabel('Time lag (s)'); title(['TMI: ',num2str(NP_unit_FS(tmp_idx(i)).TMI(1),2), ', p = ',num2str(NP_unit_FS(tmp_idx(i)).p_val_TMI,1)]);
    %         hold on; plot(repmat(1/f_wBeats*[-3:3],2,1),repmat(ylim,7,1)','r--');
    %         ylim([0 1.3*max(NP_unit_FS(tmp_idx(i)).AC(NP_unit_FS(tmp_idx(i)).AC_bins>bin_TMI))]);
    %
    %     end
    %     sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    %     fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% PHASE LOCKING
%while 0
for hide=1
    
    n_bins_phase = 20;
    phase_bins = linspace(-pi,pi,n_bins_phase);
    phase_ctrs = phase_bins(1:end-1)+mean(diff(phase_bins))/2;
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %=== Calculate phase locking for all neurons
    for nc=1:n_cells
        
        %=== Real data
        spk_phase2wbt_F = wbt_phase(s_smp_flight_hi{nc,1});   power_tmrAtspike_F =  tmr_power(s_smp_flight_hi{nc,1});  spk_phase2wbt_F = spk_phase2wbt_F(power_tmrAtspike_F>prctile(tmr_power,min_prctile_power_wbt));
        spk_phase2tht_F = tht_phase(s_smp_flight_hi{nc,1});   power_thtAtspike_F =  tht_power(s_smp_flight_hi{nc,1});  spk_phase2tht_F = spk_phase2tht_F(power_thtAtspike_F>prctile(tht_power,min_prctile_power));
        spk_phase2tmr_F = tmr_phase(s_smp_flight_hi{nc,1});   power_tmrAtspike_F =  tmr_power(s_smp_flight_hi{nc,1});  spk_phase2tmr_F = spk_phase2tmr_F(power_tmrAtspike_F>prctile(tmr_power,min_prctile_power));
        spk_phase2tht_A = tht_phase(s_smp_hi{nc,1});          power_thtAtspike_A =  tht_power(s_smp_hi{nc,1});         spk_phase2tht_A = spk_phase2tht_A(power_thtAtspike_A>prctile(tht_power,min_prctile_power));
        spk_phase2tmr_A = tmr_phase(s_smp_hi{nc,1});          power_tmrAtspike_A =  tmr_power(s_smp_hi{nc,1});         spk_phase2tmr_A = spk_phase2tmr_A(power_tmrAtspike_A>prctile(tmr_power,min_prctile_power));
        
        NP_unit(nc).spk_phase2wbt_F = spk_phase2wbt_F;  % Spike phase to wingbeat (Flight)
        NP_unit(nc).spk_phase2tht_F = spk_phase2tht_F;  % Spike phase to Theta (Flight)
        NP_unit(nc).spk_phase2tmr_F = spk_phase2tmr_F;  % Spike phase to Tamir's LFP (Flight)
        NP_unit(nc).spk_phase2tht_A = spk_phase2tht_A;  % Spike phase to Theta (All)
        NP_unit(nc).spk_phase2tmr_A = spk_phase2tmr_A;  % Spike phase to Tamir's LFP (All)
        
        %=== Normalized counts in each phase bin
        NP_unit(nc).spk_phase2wbt_F_counts = histcounts(spk_phase2wbt_F,phase_bins,'Normalization','probability')';
        NP_unit(nc).spk_phase2tht_F_counts = histcounts(spk_phase2tht_F,phase_bins,'Normalization','probability')';
        NP_unit(nc).spk_phase2tmr_F_counts = histcounts(spk_phase2tmr_F,phase_bins,'Normalization','probability')';
        NP_unit(nc).spk_phase2tht_A_counts = histcounts(spk_phase2tht_A,phase_bins,'Normalization','probability')';
        NP_unit(nc).spk_phase2tmr_A_counts = histcounts(spk_phase2tmr_A,phase_bins,'Normalization','probability')';
        
        NP_unit(nc).all_phase2wbt_counts = histcounts(wbt_phase,phase_bins,'Normalization','probability')';
        NP_unit(nc).all_phase2tht_counts = histcounts(tht_phase(tht_power>prctile(tht_power,min_prctile_power)),phase_bins,'Normalization','probability')';
        NP_unit(nc).all_phase2tmr_counts = histcounts(tmr_phase(tmr_power>prctile(tmr_power,min_prctile_power)),phase_bins,'Normalization','probability')';
        NP_unit(nc).phase_bins = phase_bins;
        
    end
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    %======================================================================================= Reviewer Phase Locking =======================
    spk_phase2wbt_C = wbt_phase(bflying_Hi);                    % All wingbeat phases during flight
    spk_phase2tmr_C = tmr_phase(bflying_Hi);                    % All Tamir's phases during flight
    power_tmrAtspike_C =  tmr_power(bflying_Hi);                % All Tamir's power during flight
    spk_phase2tmr_C = spk_phase2tmr_C(power_tmrAtspike_C>prctile(tmr_power,min_prctile_power));    % Keep only > 25% percentile
    spk_phase2wbt_C = spk_phase2wbt_C(~isnan(spk_phase2wbt_C) & power_tmrAtspike_C>prctile(tmr_power,min_prctile_power_wbt)); % Remove NaNs and keep same number of spikes as Tamir's
    
    %=== Calculate circular statistics for wingbeat and Tamir's phase
    for nc=1:n_cells
        
        if numel(s_smp_flight_hi{nc,1})>=30
            
            %=== Get wingbeat phase of spikes emitted during flight and calculate circular statistics
            spk_phase2wbt_F = wbt_phase(s_smp_flight_hi{nc,1});                     % Wingbeat phase during spikes
            power_tmrAtspike_F =  tmr_power(s_smp_flight_hi{nc,1});  
            spk_phase2wbt_F = spk_phase2wbt_F(power_tmrAtspike_F>prctile(tmr_power,min_prctile_power_wbt));
            r_spk_phase2wbt = circ_r(spk_phase2wbt_F);                              % mean resultant vector length
            m_spk_phase2wbt = circ_mean(spk_phase2wbt_F);                           % Circular mean
            n_spk_phase2wbt = numel(spk_phase2wbt_F);                               % Sample size
            [p_val_r_spk_phase2wbt, z_spk_phase2wbt] = circ_rtest(spk_phase2wbt_F);  % P value and Z of the Rayleigh test
            
            %=== Shuffling procedure
            r_spk_phase2wbt_Ctrl = zeros(1,n_rep_phase);
            for i=1:n_rep_phase
                r_spk_phase2wbt_Ctrl(i) = circ_r(datasample(spk_phase2wbt_C,n_spk_phase2wbt,'Replace',false));
            end
            if ~isnan(r_spk_phase2wbt)
                p_val_r_spk_phase2wbt_boots = sum(r_spk_phase2wbt<r_spk_phase2wbt_Ctrl)/n_rep_phase;
            else
                p_val_r_spk_phase2wbt_boots = NaN;
            end
                
            %=== Get Tamir's phase of spikes emitted during flight (only top 25% power) and calculate circular statistics
            spk_phase2tmr_F = tmr_phase(s_smp_flight_hi{nc,1});   power_tmrAtspike_F =  tmr_power(s_smp_flight_hi{nc,1});  spk_phase2tmr_F = spk_phase2tmr_F(power_tmrAtspike_F>prctile(tmr_power,min_prctile_power));
            r_spk_phase2tmr = circ_r(spk_phase2tmr_F);                              % mean resultant vector length
            m_spk_phase2tmr = circ_mean(spk_phase2tmr_F);                           % Circular mean
            n_spk_phase2tmr = numel(spk_phase2tmr_F);                               % Sample size
            [p_val_r_spk_phase2tmr, z_spk_phase2tmr] = circ_rtest(spk_phase2tmr_F);  % P value and Z of the Rayleigh test
            
            %=== Shuffling procedure
            r_spk_phase2tmr_Ctrl = zeros(1,n_rep_phase);
            for i=1:n_rep_phase
                r_spk_phase2tmr_Ctrl(i) = circ_r(datasample(spk_phase2tmr_C,n_spk_phase2tmr,'Replace',false));
            end
            if ~isnan(r_spk_phase2tmr)
                p_val_r_spk_phase2tmr_boots = sum(r_spk_phase2tmr<r_spk_phase2tmr_Ctrl)/n_rep_phase;
            else
                p_val_r_spk_phase2tmr_boots = NaN;
            end
            
            %=== Assign all values for saving
            NP_unit(nc).r_spk_phase2wbt = r_spk_phase2wbt;
            NP_unit(nc).m_spk_phase2wbt = m_spk_phase2wbt;
            NP_unit(nc).n_spk_phase2wbt = n_spk_phase2wbt;
            NP_unit(nc).p_val_r_spk_phase2wbt = p_val_r_spk_phase2wbt;
            NP_unit(nc).z_spk_phase2wbt = z_spk_phase2wbt;
            NP_unit(nc).p_val_r_spk_phase2wbt_boots = p_val_r_spk_phase2wbt_boots;
            NP_unit(nc).r_spk_phase2tmr = r_spk_phase2tmr;
            NP_unit(nc).m_spk_phase2tmr = m_spk_phase2tmr;
            NP_unit(nc).n_spk_phase2tmr = n_spk_phase2tmr;
            NP_unit(nc).p_val_r_spk_phase2tmr = p_val_r_spk_phase2tmr;
            NP_unit(nc).z_spk_phase2tmr = z_spk_phase2tmr;
            NP_unit(nc).p_val_r_spk_phase2tmr_boots = p_val_r_spk_phase2tmr_boots;
        else
            NP_unit(nc).r_spk_phase2wbt = NaN;
            NP_unit(nc).m_spk_phase2wbt = NaN;
            NP_unit(nc).n_spk_phase2wbt = NaN;
            NP_unit(nc).p_val_r_spk_phase2wbt = NaN;
            NP_unit(nc).z_spk_phase2wbt = NaN;
            NP_unit(nc).p_val_r_spk_phase2wbt_boots = NaN;
            NP_unit(nc).r_spk_phase2tmr = NaN;
            NP_unit(nc).m_spk_phase2tmr = NaN;
            NP_unit(nc).n_spk_phase2tmr = NaN;
            NP_unit(nc).p_val_r_spk_phase2tmr = NaN;
            NP_unit(nc).z_spk_phase2tmr = NaN;
            NP_unit(nc).p_val_r_spk_phase2tmr_boots = NaN;
        end
    end
    
    wbt_locked_cells = find([NP_unit.p_val_r_spk_phase2wbt_boots]<0.05);
    tmr_locked_cells = find([NP_unit.p_val_r_spk_phase2tmr_boots]<0.05);

    %=== Show some phase locked cells
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(10,20,'TileSpacing','none');
    for i=1:numel(wbt_locked_cells)
        nexttile;   polarhistogram( NP_unit(wbt_locked_cells(i)).spk_phase2wbt_F,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');  rticks([]);  thetaticks([]);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');   
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(10,20,'TileSpacing','none');
    for i=1:numel(tmr_locked_cells)
        nexttile;   polarhistogram( NP_unit(tmr_locked_cells(i)).spk_phase2tmr_F,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');  rticks([]);  thetaticks([]);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');   
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Show comparison between wingbeat locking and Tamir's locking
    if ~isempty(wbt_locked_cells)
        figure('units','normalized','outerposition',[.2 .3 .2 .4],'Visible',fig_visibility);
        tiledlayout(1,2,'TileSpacing','tight');
        nexttile;   plot_distr_AF_v0([NP_unit.r_spk_phase2wbt]', [NP_unit.r_spk_phase2tmr]', {'Wingbeat', 'Non-oscillatory'}, 'SEM', 'Mean resultant vector length (all)'); ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
        nexttile;   plot_distr_AF_v0([NP_unit(wbt_locked_cells).r_spk_phase2wbt], [NP_unit(tmr_locked_cells).r_spk_phase2tmr]', {'Wingbeat', 'Non-oscillatory'}, 'SEM', 'Mean resultant vector length (significant)'); ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    end
    %=======================================================================================  =======================
    
end

%% SPIKE TRIGGERED LFP AND WINGBEAT
%while 0
for hide=1
    
    interval = [round(-2*Fs_Hi):round(2*Fs_Hi)];
    
    %=== Spike triggered LFP (all spikes)
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        LFP_tmp = zeros(numel(interval),numel(s_smp_hi{nc,1}));
        for i=1:numel(s_smp_hi{nc,1})
            
            if all((s_smp_hi{nc,1}(i)+interval)>1) && all((s_smp_hi{nc,1}(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_hi{nc,1}(i)+interval);
            end
        end
        nexttile;   plot(interval/Fs,mean(LFP_tmp,2));  xticks([]); yticks([]);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Spike triggered LFP (control)
    figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        LFP_tmp = zeros(numel(interval),numel(s_smp_hi{nc,1}));
        s_smp_hi_sh = randi([min(s_smp_hi{nc,1}),max(s_smp_hi{nc,1})],numel(s_smp_hi{nc,1}),1);
        for i=1:numel(s_smp_hi{nc,1})
            
            if all((s_smp_hi_sh(i)+interval)>1) && all((s_smp_hi_sh(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_hi_sh(i)+interval);
            end
        end
        nexttile;   plot(interval/Fs,mean(LFP_tmp,2));  xticks([]); yticks([]);
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Accumulate (flight, all and rest)
    spkTgdLfP_flight = [];
    spkTgdWbt_flight = [];
    spkTgdLfP_all = [];
    spkTgdLfP_rest = [];
    for nc=1:n_cells
        
        LFP_tmp = zeros(numel(interval),numel(s_smp_flight_hi{nc,1}));
        for i=1:numel(s_smp_flight_hi{nc,1})
            if all((s_smp_flight_hi{nc,1}(i)+interval)>1) && all((s_smp_flight_hi{nc,1}(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_flight_hi{nc,1}(i)+interval);
            end
        end
        spkTgdLfP_flight = [spkTgdLfP_flight,mean(LFP_tmp,2)];
        NP_unit(nc).spkTgdLfP_flight = {mean(LFP_tmp,2)};
        
        LFP_tmp = zeros(numel(interval),numel(s_smp_hi{nc,1}));
        for i=1:numel(s_smp_hi{nc,1})
            if all((s_smp_hi{nc,1}(i)+interval)>1) && all((s_smp_hi{nc,1}(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_hi{nc,1}(i)+interval);
            end
        end
        spkTgdLfP_all = [spkTgdLfP_all,mean(LFP_tmp,2)];
        NP_unit(nc).spkTgdLfP_all = {mean(LFP_tmp,2)};
        
        s_smp_rest_hi = setdiff(s_smp_hi{nc,1},s_smp_flight_hi{nc,1});
        LFP_tmp = zeros(numel(interval),numel(s_smp_rest_hi));
        for i=1:numel(s_smp_rest_hi)
            if all((s_smp_rest_hi(i)+interval)>1) && all((s_smp_rest_hi(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_rest_hi(i)+interval);
            end
        end
        spkTgdLfP_rest = [spkTgdLfP_rest,mean(LFP_tmp,2)];
        NP_unit(nc).spkTgdLfP_rest = {mean(LFP_tmp,2)};
        
        WBT_tmp = zeros(numel(interval),numel(s_smp_flight_hi{nc,1}));
        for i=1:numel(s_smp_flight_hi{nc,1})
            if all((s_smp_flight_hi{nc,1}(i)+interval)>1) && all((s_smp_flight_hi{nc,1}(i)+interval)<numel(t_Hi))
                WBT_tmp(:,i) = a_abs_NP(s_smp_flight_hi{nc,1}(i)+interval);
            end
        end
        spkTgdWbt_flight = [spkTgdWbt_flight,mean(WBT_tmp,2)];
        NP_unit(nc).spkTgdWbt_flight = {mean(WBT_tmp,2)};
        
    end
    
    %=== Plot
    figure('units','normalized','outerposition',[.2 .1 .4 .3],'Visible',fig_visibility);
    tiledlayout(1,4,'TileSpacing','tight');
    nexttile;   data_tmp = spkTgdLfP_flight';
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('uV');   title('Flight');
    nexttile;   data_tmp = spkTgdLfP_all';
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('uV');   title('All');
    nexttile;   data_tmp = spkTgdLfP_rest';
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('uV');   title('Rest');
    nexttile;   data_tmp = spkTgdWbt_flight';
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('g');   title('Wingbeat');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    LFP_PSD.spkTgdLfP_time = {interval/Fs_Hi};
    LFP_PSD.spkTgdLfP_all = {mean(spkTgdLfP_all','omitnan')};
    LFP_PSD.spkTgdLfP_flight = {mean(spkTgdLfP_flight','omitnan')};
    LFP_PSD.spkTgdWbt_flight = {mean(spkTgdWbt_flight','omitnan')};
    LFP_PSD.spkTgdLfP_rest = {mean(spkTgdLfP_rest','omitnan')};
    
end

%% SHOW EXAMPLE PLACE CELLS AND LOOK AT PHASE PRECESSION AND AUTOCORRELATIONS
while 0
%for hide=1
    
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
        
        %=== Get the subtable of place cells
        NP_table = struct2table(NP_unitOnClus{1, j});
        
        %=== For template matching (aka stable place cells)
        NP_table.place_cond = NP_table.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
            NP_table.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
            NP_table.stab_m>0.4 &...      % Min Stability (DEF: .4)
            NP_table.sff<0.7;             % Min Peakyness (DEF: .7)
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
        smp1_clus_hi = smp1_clus;   smp2_clus_hi = smp2_clus;               % Samples at higher sampling frequency
        for nn = 1:n_fclus
            smp1_clus_hi(nn) = knnsearch_fast_AF_v0(t_Hi',t(smp1_clus(nn)),0.1);
            smp2_clus_hi(nn) = knnsearch_fast_AF_v0(t_Hi',t(smp2_clus(nn)),0.1);
        end
        id = find(f_clus.id==j);                                            % ids of the flight clusters
        avg_takeoff = mean(squeeze(f_clus.ds_pos(:,1,id)),2);               % Average takeoff position
        avg_landing = mean(squeeze(f_clus.ds_pos(:,end,id)),2);             % Average landing position
        d_fclus = mean(t(smp2_clus)-t(smp1_clus));                          % Mean duration of the flight
        all_s_plc = sort(vertcat(s{NP_subtable.cell_id,1}));                % Define rate by pooling spikes from place cells
        all_rate_plc = kernel_rate_AF_v1(all_s_plc,0.1,t)/n_place_cells;
        t_ic = [-1 1];                                                      % Time for plotting around the flight
        
        %=== Target sequence determination
        NP_subtable_sorted = sortrows(NP_subtable,'phase_max');
        sorted_tko_plcCells = NP_subtable_sorted.cell_id;
        
        %=== Get the spike times happening around flight (plus t_ic)
        s_flight = cell(n_place_cells,1);       % Spikes during the flight cluster (sorted)
        for nc = 1:n_place_cells,[~,s_flight{nc,1}] = count_spikes_AF_v1(s{sorted_tko_plcCells(nc),1},t,[smp1_clus smp2_clus]);end
        s_flight1 = cell(n_cells,1);            % Spikes during the flight cluster (unsorted, extended)
        for nc = 1:n_cells,[~,s_flight1{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);end
        
        % Get the flight phase at each t sample
        lin_flight = zeros(T,1);
        len_flight = zeros(T,1);
        for nn=1:size(fc_lin,1),lin_flight(smp1_clus(nn):smp2_clus(nn)) = fc_lin{nn,1};end
        for nn=1:size(fc_lin,1),len_flight(smp1_clus(nn):smp2_clus(nn)) = fc_all{nn,1};end
        
        %=============================================================================================================
        
        %=== Params
        find_optimal_phase = 1; % Specifies if adjusting the phase for min correlation
        n_bins_phase = 10;
        phase_bins = linspace(-pi,pi,n_bins_phase);
        phase_ctrs = phase_bins(1:end-1)+mean(diff(phase_bins))/2;
        
        %=== Get the spike samples for this cluster
        s_smp_flight_lo_clus = cell(n_place_cells,1);% Spike samples during the flight cluster (sorted), low sampling
        s_smp_flight_hi_clus = cell(n_place_cells,1);% Spike samples during the flight cluster (sorted), hi sampling
        for nc = 1:n_place_cells
            [~,~,valid_samples] = count_spikes_AF_v2(s{sorted_tko_plcCells(nc),1},t_Hi,[smp1_clus_hi smp2_clus_hi]);
            s_smp_flight_hi_clus{nc,1} = s_smp_hi{sorted_tko_plcCells(nc),1}(valid_samples);
            s_smp_flight_lo_clus{nc,1} = round(m_slope * s_smp_flight_hi_clus{nc,1}+q_intcp);
        end
        
        %=== Get the control phase distributions during the specific flight cluster
        wbt_phase_ctrl = [];
        tht_phase_ctrl = [];
        tmr_phase_ctrl = [];
        for nn = 1:n_fclus
            wbt_phase_tmp  = wbt_phase(smp1_clus_hi(nn):smp2_clus_hi(nn));
            tht_phase_tmp  = tht_phase(smp1_clus_hi(nn):smp2_clus_hi(nn));   power_thtAtspike =  tht_power(smp1_clus_hi(nn):smp2_clus_hi(nn));   tht_phase_tmp = tht_phase_tmp(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
            tmr_phase_tmp  = tmr_phase(smp1_clus_hi(nn):smp2_clus_hi(nn));   power_tmrAtspike =  tmr_power(smp1_clus_hi(nn):smp2_clus_hi(nn));   tmr_phase_tmp = tmr_phase_tmp(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
            wbt_phase_ctrl = [wbt_phase_ctrl; wbt_phase_tmp];
            tht_phase_ctrl = [tht_phase_ctrl; tht_phase_tmp];
            tmr_phase_ctrl = [tmr_phase_ctrl; tmr_phase_tmp];
        end
        
        %=== Initialize the cells containing, for each cell and spike the...
        spk2plc_dist_wbt = s_smp_flight_lo_clus;    %... distance to place field
        spk2plc_dist_tht = s_smp_flight_lo_clus;    %... distance to place field
        spk2plc_dist_tmr = s_smp_flight_lo_clus;    %... distance to place field
        spk2wbt_phase = s_smp_flight_lo_clus;   %... wingbeat phase of the spike
        spk2tht_phase = s_smp_flight_lo_clus;   %... theta phase of the spike
        spk2tmr_phase = s_smp_flight_lo_clus;   %... Tamir's phase of the spike
        
        phase_wbt_place_corr = zeros(n_place_cells,1);
        phase_tht_place_corr = zeros(n_place_cells,1);
        phase_tmr_place_corr = zeros(n_place_cells,1);
        
        %=== Extract a few variables for phase precession (NEW)
        for nc = 1:n_place_cells
            
            %d2place = lin_flight - NP_subtable_sorted.phase_max(nc)*0.01;
            d2place = len_flight;
            
            spk_phase2wbt = wbt_phase(s_smp_flight_hi_clus{nc,1});
            spk_phase2tht = tht_phase(s_smp_flight_hi_clus{nc,1});  power_thtAtspike =  tht_power(s_smp_flight_hi_clus{nc,1});  spk_phase2tht = spk_phase2tht(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
            spk_phase2tmr = tmr_phase(s_smp_flight_hi_clus{nc,1});  power_tmrAtspike =  tmr_power(s_smp_flight_hi_clus{nc,1});  spk_phase2tmr = spk_phase2tmr(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
            
            spk2plc_dist_wbt{nc,1} = d2place(s_smp_flight_lo_clus{nc,1});
            spk2plc_dist_tht{nc,1} = d2place(s_smp_flight_lo_clus{nc,1});   spk2plc_dist_tht{nc,1} = spk2plc_dist_tht{nc,1}(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
            spk2plc_dist_tmr{nc,1} = d2place(s_smp_flight_lo_clus{nc,1});   spk2plc_dist_tmr{nc,1} = spk2plc_dist_tmr{nc,1}(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
            spk2wbt_phase{nc,1} = spk_phase2wbt;    spk2wbt_phase_opt{nc,1} = spk_phase2wbt;
            spk2tht_phase{nc,1} = spk_phase2tht;    spk2tht_phase_opt{nc,1} = spk_phase2tht;
            spk2tmr_phase{nc,1} = spk_phase2tmr;    spk2tmr_phase_opt{nc,1} = spk_phase2tmr;
            
            if find_optimal_phase
                
                delta_phase = linspace(-pi,pi,20);
                pp_corr = NaN(numel(delta_phase),1);
                pp_pval = NaN(numel(delta_phase),1);
                
                %=== Wingbeat
                temp_phase = wbt_phase(s_smp_flight_hi_clus{nc,1});
                for dtp = 1:numel(delta_phase)
                    [pp_corr(dtp),pp_pval(dtp)] = corr(spk2plc_dist_wbt{nc,1},wrapToPi(temp_phase+delta_phase(dtp)),'type','Spearman');
                end
                [phase_wbt_place_corr(nc,1),bst_dtp] = max(abs(pp_corr));
                spk2wbt_phase_opt{nc,1} = wrapToPi(temp_phase+delta_phase(bst_dtp));
                
                %[pp_rho, pp_pval] = circ_corrcl(wbt_phase(s_smp_flight_hi_clus{nc,1}),spk2plc_dist_wbt{nc,1});
                NP_unit(sorted_tko_plcCells(nc)).f_clus(j).spk2wbt_rho = pp_corr(bst_dtp);
                NP_unit(sorted_tko_plcCells(nc)).f_clus(j).spk2wbt_pvl = pp_pval(bst_dtp);
                NP_unit(sorted_tko_plcCells(nc)).f_clus(j).spk2wbt_dtp = delta_phase(bst_dtp);
                
                %=== Shuffling for wingbeat
                shfl_corr = NaN(100,1);
                for ssf = 1:100
                    shfl_phase = -pi + (2*pi)*rand(1, numel(spk_phase2wbt));
                    shfl_corr(ssf) = corr(spk2plc_dist_wbt{nc,1},wrapToPi(shfl_phase+delta_phase(bst_dtp))');
                end
                NP_unit(sorted_tko_plcCells(nc)).f_clus(j).spk2wbt_pvl_sh = sum(shfl_corr<pp_corr(bst_dtp))/100;
                
                %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                %                 %=== Theta
                %                 temp_phase = spk_phase2tht;
                %                 for dtp = 1:numel(delta_phase)
                %                     pp_corr(dtp) = corr(spk2plc_dist_tht{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                %                 end
                %                 [phase_tht_place_corr(nc,1),bst_dtp] = min(pp_corr);
                %                 spk2tht_phase_opt{nc,1} = wrapToPi(temp_phase+delta_phase(bst_dtp));
                %
                %                 %=== Tamir
                %                 temp_phase = spk_phase2tmr;
                %                 for dtp = 1:numel(delta_phase)
                %                     pp_corr(dtp) = corr(spk2plc_dist_tmr{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                %                 end
                %                 [phase_tmr_place_corr(nc,1),bst_dtp] = min(pp_corr);
                %                 spk2tmr_phase_opt{nc,1} = wrapToPi(temp_phase+delta_phase(bst_dtp));
                %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
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
        
        %=== Visualize preferred phases
        pop_phase2wbt = mean(cell2mat(cellfun(@(x) histcounts(x,phase_bins,'Normalization','probability'),spk2wbt_phase,'UniformOutput',false)),1,'omitnan');
        pop_phase2tht = mean(cell2mat(cellfun(@(x) histcounts(x,phase_bins,'Normalization','probability'),spk2tht_phase,'UniformOutput',false)),1,'omitnan');
        pop_phase2tmr = mean(cell2mat(cellfun(@(x) histcounts(x,phase_bins,'Normalization','probability'),spk2tmr_phase,'UniformOutput',false)),1,'omitnan');
        ctrl_phase2wbt = histcounts(wbt_phase_ctrl,phase_bins,'Normalization','probability');
        ctrl_phase2tht = histcounts(tht_phase_ctrl,phase_bins,'Normalization','probability');
        ctrl_phase2tmr = histcounts(tmr_phase_ctrl,phase_bins,'Normalization','probability');
        
        figure('units','normalized','outerposition',[.2 .2 0.2 0.4],'Visible',fig_visibility);
        tiledlayout(3,3,'TileSpacing','tight');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[pop_phase2wbt,pop_phase2wbt],1,'k','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[pop_phase2tht,pop_phase2tht],1,'r','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Theta (F)');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[pop_phase2tmr,pop_phase2tmr],1,'g','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Tamir (F)');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[ctrl_phase2wbt,ctrl_phase2wbt],1,'k','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Control');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[ctrl_phase2tht,ctrl_phase2tht],1,'r','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Control');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],[ctrl_phase2tmr,ctrl_phase2tmr],1,'g','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Control');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],repmat(pop_phase2wbt-ctrl_phase2wbt,1,2),1,'k','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Diff');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],repmat(pop_phase2tht-ctrl_phase2tht,1,2),1,'r','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Diff');
        nexttile;   bar([phase_ctrs,phase_ctrs+2*pi],repmat(pop_phase2tmr-ctrl_phase2tmr,1,2),1,'g','EdgeColor','none');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Diff');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Average autocorrelogram
        cum_ac = [];
        for nc = 1:n_place_cells
            [ac_TMI_tmp,ac_TMI_bins] = cross_correlogram_AF_v0(s_flight{nc,1},s_flight{nc,1},max_lag_TMI,bin_TMI);
            cum_ac = [cum_ac,ac_TMI_tmp];
        end
        figure('units','normalized','outerposition',[0.3 .3 0.17 .4],'Visible',fig_visibility);
        area(ac_TMI_bins,mean(cum_ac,2),'FaceColor','k','EdgeColor','none','FaceAlpha',0.5);
        yticks([]); xlabel('Time lag (s)'); hold on;
        plot(repmat(1/f_wBeats*[-3:3],2,1),repmat(ylim,7,1)','k--');
        xlim([0 max_lag_TMI]);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Visualize accelerometer profile
        [~,tmp_idx] = sort(smp2_clus-smp1_clus);
        a_abs_all = zeros(max(smp2_clus_hi-smp1_clus_hi),n_fclus);
        for ii=1:n_fclus
            jj = tmp_idx(ii);
            a_abs_all(1:(smp2_clus_hi(jj)-smp1_clus_hi(jj)+1),ii) = a_flt_NP(smp1_clus_hi(jj):smp2_clus_hi(jj));
        end
        flight_pairs = [ones(n_fclus,1),[1:n_fclus]'];  opt_lag = zeros(n_fclus,1);
        for nn = 1:size(flight_pairs,1)
            [~,opt_lag(nn)] = max(xcorr(a_abs_all(:,flight_pairs(nn,1)),a_abs_all(:,flight_pairs(nn,2)),round(Fs_Hi/10)));
        end
        opt_lag=opt_lag-round(Fs_Hi/10);
        a_abs_all_sh = circshift(a_abs_all,opt_lag);
        figure('units','normalized','outerposition',[0.6 .3 0.3 .4],'Visible',fig_visibility);
        tiledlayout(1,2,'TileSpacing','tight');
        nexttile;   imagesc([0 size(a_abs_all,1)/Fs_Hi],[],a_abs_all');    colormap(gray);
        xlabel('Time (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
        nexttile;   imagesc([0 size(a_abs_all_sh,1)/Fs_Hi],[],a_abs_all_sh');    colormap(gray);
        xlabel('Time (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
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
        
        %=== Plot place cells around the selected cluster (+ Wingbeat Phase)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(6,4*ceil(n2show/6),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable_sorted.cell_id(ncc);
            fld_ctr = NP_subtable_sorted.phase_max(ncc)*0.01*NP_subtable_sorted.f_lenght(ncc);
            
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            scatter([spk2plc_dist_wbt{ncc,1};spk2plc_dist_wbt{ncc,1}],[spk2wbt_phase_opt{ncc,1};spk2wbt_phase_opt{ncc,1}+2*pi],6,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist_wbt{ncc,1},[1 99])*1);
            %title(num2str(NP_unit(sorted_tko_plcCells(ncc)).f_clus(j).spk2wbt_pvl_sh,2));
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Wingbeat Phase histogram)
        n_bins_phase = 10;
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            histogram([spk2wbt_phase{ncc,1};spk2wbt_phase{ncc,1}+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','k');
            xticks(pi*[-1 0 1 2 3]);  xticklabels({'-180', '0', '180', '360', '540'});
            xlabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Theta Phase)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(6,4*ceil(n2show/6),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable_sorted.cell_id(ncc);
            fld_ctr = NP_subtable_sorted.phase_max(ncc)*0.01*NP_subtable_sorted.f_lenght(ncc);
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            scatter([spk2plc_dist_tht{ncc,1};spk2plc_dist_tht{ncc,1}],[spk2tht_phase_opt{ncc,1};spk2tht_phase_opt{ncc,1}+2*pi],6,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist_tht{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Theta Phase histogram)
        n_bins_phase = 10;
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            histogram([spk2tht_phase{ncc,1};spk2tht_phase{ncc,1}+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','r');
            xticks(pi*[-1 0 1 2 3]);  xticklabels({'-180', '0', '180', '360', '540'});
            xlabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Tamir Phase)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(6,4*ceil(n2show/6),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable_sorted.cell_id(ncc);
            fld_ctr = NP_subtable_sorted.phase_max(ncc)*0.01*NP_subtable_sorted.f_lenght(ncc);
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            scatter([spk2plc_dist_tmr{ncc,1};spk2plc_dist_tmr{ncc,1}],[spk2tmr_phase_opt{ncc,1};spk2tmr_phase_opt{ncc,1}+2*pi],6,'filled','MarkerFaceColor',[0 .6 .3],'MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist_tmr{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Tamir Phase histogram)
        n_bins_phase = 10;
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
            xticks([]); yticks([]); title(['Unit ' num2str(nc), ', Pt ', num2str(ncc) ]);
            nexttile([1 2]);
            histogram([spk2tmr_phase{ncc,1};spk2tmr_phase{ncc,1}+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor',[0 .6 .3]);
            xticks(pi*[-1 0 1 2 3]);  xticklabels({'-180', '0', '180', '360', '540'});
            xlabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Autocorrelation)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(6,4*ceil(n2show/6),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable_sorted.cell_id(ncc);
            
            %=== Get the spikes happening during cluster j (add 0.5s at the flight tails)
            [ac_TMI,ac_TMI_bins] = cross_correlogram_AF_v0(s_flight{ncc,1},s_flight{ncc,1},max_lag_TMI,bin_TMI);
            
            %=== Plot
            nexttile([1 2]);   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],.5,[],1);
            nexttile([1 2]);    hold on;
            area(ac_TMI_bins,ac_TMI./max(ac_TMI),'FaceColor','k','EdgeColor','none','FaceAlpha',0.5);
            area(NP_unit(nc).AC_bins,(NP_unit(nc).AC)./max(NP_unit(nc).AC),'FaceColor',col_clus(j,:),'EdgeColor','k','FaceAlpha',0.5);
            %max_tmp = max( NP_unit(nc).AC(NP_unit(nc).AC_bins>bin_TMI)/max(NP_unit(nc).AC));
            yticks([]); xlabel('Time lag (s)');     title(['Unit ' num2str(nc)]);
            %title(['TMI: ',num2str(NP_unit(nc).TMI(1),2), ', p = ',num2str(NP_unit(nc).p_val_TMI,1)]);
            plot(repmat(1/f_wBeats*[-3:3],2,1),repmat(ylim,7,1)','k--');
            xlim([0 max_lag_TMI]);
            
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
        
    end
end

%% FLIGHT DECODING AND THETA-SWEEPS
while 0
%for hide=1
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
        
        %=== Get the subtable of cells for template matching and cells for decoding
        NP_table = struct2table(NP_unitOnClus{1, j});
        
        %=== Extract some variables for plotting
        smp1_clus = f_clus.strt_frame(f_clus.id == j)';                     % Takeoff sample
        smp2_clus = f_clus.stop_frame(f_clus.id == j)';                     % Landing sample
        cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
        smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
        n_fclus = numel(smp1_clus);                                         % Number of surviving flights in the cluster
        smp1_clus_hi = smp1_clus;   smp2_clus_hi = smp2_clus;               % Samples at higher sampling frequency
        for nn = 1:n_fclus
            smp1_clus_hi(nn) = knnsearch_fast_AF_v0(t_Hi',t(smp1_clus(nn)),0.1);
            smp2_clus_hi(nn) = knnsearch_fast_AF_v0(t_Hi',t(smp2_clus(nn)),0.1);
        end
        id = find(f_clus.id==j);                                            % ids of the flight clusters
        avg_takeoff = mean(squeeze(f_clus.ds_pos(:,1,id)),2);               % Average takeoff position
        avg_landing = mean(squeeze(f_clus.ds_pos(:,end,id)),2);             % Average landing position
        d_fclus = mean(t(smp2_clus)-t(smp1_clus));                          % Mean duration of the flight
        flight_bin_centers = NP_table.plc_ctr{1,1};                         % Positions at bin centers (useful later)
        
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
        
        %=====================================================================================================
        %=== For decoding FLIGHT ACTIVITY
        t_bin_dur = 0.030;                                                  % Time bin duration for decoding (def 25ms)
        t_bin_ovl = t_bin_dur-0.005;                                        % Time bin overlap for decoding (t_bin_dur-0.005)
        smooth_f = [1 .3];                                                  % [space bin,time bin]
        prc_lim = [1 98.5];
        
        %=== Cells used for decoding flights
        NP_table.use4decoding = NP_table.peakHz>0 & NP_table.fr<5 & NP_table.corr_m>0.6 & NP_table.spkPerflight>0;
        NP_subtable_dec = NP_table(NP_table.use4decoding,:);
        n_dec_cells = size(NP_subtable_dec,1);
        all_s_dec = sort(vertcat(s{NP_subtable_dec.cell_id,1}));                % Define rate by pooling spikes from cells used for decoding
        all_rate_dec = kernel_rate_AF_v1(all_s_dec,0.01,t)/n_dec_cells;         % Spike density from cells used for decoding
        
        %=== Ordering of cells, useful later for cross-correlation
        NP_subtable_sorted = sortrows(NP_subtable_dec,'phase_max');
        sorted_tko_plcCells = NP_subtable_sorted.cell_id;
        
        %=== DECODING Initialization
        n_bins = floor((t(T)-t_bin_dur)/(t_bin_dur-t_bin_ovl))+1;   % Number of time bins
        st_times = (0:n_bins-1)*(t_bin_dur - t_bin_ovl);            % Start times
        ed_times = st_times + t_bin_dur;                            % End times
        ct_times = st_times + t_bin_dur/2;                          % Center times
        t_d = zeros(1, 2 * n_bins);                                 % All times (initialize)
        t_d(1:2:end) = st_times;    t_d(2:2:end) = ed_times;        % All times (define)
        n_vector = zeros(1,n_dec_cells);                            % Initialize vector containing the number of spikes fired by each place cell
        p_x = NP_unitOnClus{1,j}(1).prob_x;                         % Probability of being on a given bin of the flight cluster
        f_x = zeros(numel(p_x),n_dec_cells);                        % Vector with the firing rates of place cells in a given bin of the flight cluster
        for nc=1:n_dec_cells, f_x(:,nc) = NP_subtable_dec.plc_map{nc,1};end
        tht_SWPS = struct();                                        % Initialize the tht_SWPS structure, for holding the posterior probability and other relevant variables
        
        %=== Generate mirrored spikes for smoothing (2 on the left and 2 on the right of each spike)
        s_smooth_1 = cellfun(@(x) x+0.005,s(:,1),'Uniformoutput',false);
        s_smooth_2 = cellfun(@(x) x+0.010,s(:,1),'Uniformoutput',false);
        s_smooth_3 = cellfun(@(x) x-0.005,s(:,1),'Uniformoutput',false);
        s_smooth_4 = cellfun(@(x) x-0.010,s(:,1),'Uniformoutput',false);
        s_smooth = cellfun(@(x,y,z,t,u) sort([x;y;z;t;u]),s(:,1),s_smooth_1,s_smooth_2,s_smooth_3,s_smooth_4,'Uniformoutput',false);
        %s_smooth = cellfun(@(x,y,z) sort([x;y;z]),s(:,1),s_smooth_1,s_smooth_3,'Uniformoutput',false);
        
        %=== Find peaks of the wingbeat phase (plotting)
        [~,wb_max_t] = findpeaks(wbt_phase,t_Hi,'MinPeakDistance',0.1);
        
        %=== Show decoded flights
        figure('units','normalized','outerposition',[0 0 1 1],'Visible',fig_visibility);
        tiledlayout(5,6,'TileSpacing','compact');
        for zz=1:n_fclus
            
            i_bin_strt = knnsearch_fast_AF_v0(st_times',t(smp1_clus(zz)),0.5);  % index of the start bin
            i_bin_stop = knnsearch_fast_AF_v0(st_times',t(smp2_clus(zz)),0.5);  % index of the stop bin
            n_bins_tmp = i_bin_stop-i_bin_strt+1;
            
            tht_SWPS(zz).unique_ID = unique_ID;
            tht_SWPS(zz).flight_id = j;
            tht_SWPS(zz).p_dec_flight = zeros(size(p_x,1),n_bins_tmp);  % Initialize posterior probability
            tht_SWPS(zz).p_dec_shifted = zeros(size(p_x,1),n_bins_tmp); % Shifted by the real position
            tht_SWPS(zz).n_spikes = zeros(n_dec_cells,n_bins_tmp);      % Initialize average number of spikes
            tht_SWPS(zz).pos_real = zeros(n_bins_tmp,1);                % Initialize real position
            tht_SWPS(zz).spk_dsty = zeros(n_bins_tmp,1);                % Initialize spike density
            tht_SWPS(zz).bin_time = zeros(n_bins_tmp,1);                % Initialize bin time
            tht_SWPS(zz).pos_bin_real = zeros(n_bins_tmp,1);            % Initialize real pos bin
            tht_SWPS(zz).wbt_phase = zeros(n_bins_tmp,1);               % Wingbeat phase
            tht_SWPS(zz).wbt_power = zeros(n_bins_tmp,1);               % Wingbeat power
            tht_SWPS(zz).tmr_phase = zeros(n_bins_tmp,1);               % Tamir's phase
            tht_SWPS(zz).tmr_power = zeros(n_bins_tmp,1);               % Tamir's power
            tht_SWPS(zz).LFP = zeros(n_bins_tmp,1);                     % LFP value
            tht_SWPS(zz).wbt = zeros(n_bins_tmp,1);                     % Accelerometer
            tht_SWPS(zz).clk = zeros(n_bins_tmp,1);                     % Number of clicks
            
            counter = 1;
            
            %=== Loop across bins
            for i_dec = i_bin_strt:i_bin_stop
                
                %for nc=1:n_dec_cells,n_vector(nc) = histcounts(s{NP_subtable_dec.cell_id(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
                for nc=1:n_dec_cells,n_vector(nc) = histcounts(s_smooth{NP_subtable_dec.cell_id(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
                p_x = double(p_x>-1);   % Uniform Prior
                tht_SWPS(zz).p_dec_flight(:,counter) = decode_1Dpos_AF_v1(n_vector,t_bin_dur,p_x,f_x);
                tht_SWPS(zz).n_spikes(:,counter) = n_vector';
                
                %=== Find current position (linearized)
                smpTotko = round((ct_times(i_dec)-t(smp1_clus(zz)))*f_clus.Fs);                     % Sample from takeoff (at bin center)
                if smpTotko==0,smpTotko=1;end                                                       % Keep it bounded between the trajectory start/stop
                if smpTotko>numel(fc_all{zz}),smpTotko=numel(fc_all{zz});end
                
                %=== Add a few more features
                tht_SWPS(zz).pos_real(counter) = fc_all{zz}(smpTotko);                                      % Corresponding position (distance from tko) at that sample
                tht_SWPS(zz).spk_dsty(counter) = all_rate_dec(knnsearch_fast_AF_v0(t,ct_times(i_dec),0.1)); % Spike density
                tht_SWPS(zz).bin_time(counter) = ct_times(i_dec);                                           % Time at the bin center
                tht_SWPS(zz).pos_bin_real(counter) = knnsearch(flight_bin_centers,fc_all{zz}(smpTotko));    % Real positional bin
                
                tmp_idx = knnsearch_fast_AF_v0(t_Hi',ct_times(i_dec),0.01);
                
                tht_SWPS(zz).wbt_phase(counter) = wbt_phase(tmp_idx);
                tht_SWPS(zz).wbt_power(counter) = wbt_power(tmp_idx);
                tht_SWPS(zz).tmr_phase(counter) = tmr_phase(tmp_idx);
                tht_SWPS(zz).tmr_power(counter) = tmr_power(tmp_idx);
                tht_SWPS(zz).LFP(counter) = LFP(tmp_idx);
                tht_SWPS(zz).wbt(counter) = a_flt_NP(tmp_idx);
                
                tht_SWPS(zz).clk(counter) = sum(Detected_Clicks.times>t_d(2*i_dec-1) & Detected_Clicks.times<t_d(2*i_dec));
                
                
                
                %=== Shift the posterior probability (correct for the actual position)
                tht_SWPS(zz).p_dec_shifted(:,counter)  = circshift(tht_SWPS(zz).p_dec_flight(:,counter) ,-tht_SWPS(zz).pos_bin_real(counter)+round(numel(p_x)/2));
                
                counter = counter+1;
                
            end
            
            %=== Quantify decoding error and the fraction decoded bins
            [~,dec_bin] = max(tht_SWPS(zz).p_dec_flight);
            dec_error = flight_bin_centers(dec_bin)-tht_SWPS(zz).pos_real;
            dec_error(~(sum(tht_SWPS(zz).p_dec_flight)>0)) = NaN;
            tht_SWPS(zz).rmsDec_error = rms(dec_error(~isnan(dec_error)));
            tht_SWPS(zz).prc_decoded = 1-sum(sum(tht_SWPS(zz).p_dec_flight)==0)/size(tht_SWPS(zz).p_dec_flight,2);
            
            %=== Plot the first 30 flights (or less)
            if zz<31
                nexttile;   imagesc([t(smp1_clus(zz)) t(smp2_clus(zz))],[flight_bin_centers(1) flight_bin_centers(end)],imgaussfilt(tht_SWPS(zz).p_dec_flight,smooth_f),prctile(tht_SWPS(zz).p_dec_flight, prc_lim,'all')');
                colormap(flipud(gray)); ylabel('Spatial bin'); axis off;  set(gca,'YDir','normal');
                xlabel('Temporal bin');
                hold on;   plot(t_Hi,normalize(a_flt_NP,'range',[0 flight_bin_centers(end)/4])); xlim([t(smp1_clus(zz)) t(smp2_clus(zz))]);
                plot(tht_SWPS(zz).bin_time,tht_SWPS(zz).pos_real,'r');    title(['RMS (m): ', num2str(tht_SWPS(zz).rmsDec_error,2)]);
                plot(repmat(wb_max_t,2,1),[zeros(size(wb_max_t));flight_bin_centers(end)*ones(size(wb_max_t))],'k--');
                plot(repmat(Detected_Clicks.times',2,1),ones(2,numel(Detected_Clicks.times)).*[0; flight_bin_centers(end)/4],'b');
                plot(tht_SWPS(zz).bin_time,normalize(tht_SWPS(zz).spk_dsty,'range',[0 flight_bin_centers(end)/4]),'g');
            end
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Add to the saved table
        SWP_table = [SWP_table; struct2table(tht_SWPS)];
        
        %======================================================================================================================================================================
        %=== Shuffling and modulation scores
        tht_SWPS_struct = struct2table(tht_SWPS);
        SWPs_sst = tht_SWPS_struct(tht_SWPS_struct.rmsDec_error<2 & tht_SWPS_struct.prc_decoded>0.5,:); % Select subset of good flights
        
        n_sweeps_sh = 3;
        n_reshape = 50;
        SWP_mod_score = zeros(n_sweeps_sh+1,1);
        avg_est_dist_sh = zeros(n_reshape,n_sweeps_sh+1);
        sem_est_dist_sh = zeros(n_reshape,n_sweeps_sh+1);
        rsz_bin_ctrs = linspace(-180,180,n_reshape);
        N_flights = size(SWPs_sst,1);
        max_int_shift = round(0.120/(t_bin_dur-t_bin_ovl));
        
        if ~isempty(SWPs_sst)
            for iii=1:n_sweeps_sh+1
                
                %=== Params and initialize relevant variables
                single_SWP = table();
                counter=1;
                warning('off');
                
                %=== Random shift
                if iii == 1,    sweep_shuffling = 0;
                else,           sweep_shuffling = 1;    end
                
                %=== Cut at wingbeat maxima
                for zz=1:N_flights
                    
                    zero_phs_idx = find(SWPs_sst.wbt_phase{zz,1}(1:end-1).* SWPs_sst.wbt_phase{zz,1}(2:end)<0 & diff(SWPs_sst.wbt_phase{zz,1})<0);  % Segment based on wingbeat
                    sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
                    
                    if sweep_shuffling
                        rand_shift =  -round(max_int_shift/2)+randi(max_int_shift,numel(sweep_strt),1);
                        sweep_strt_sh = sweep_strt+rand_shift; sweep_stop_sh = sweep_stop+rand_shift;
                        sweep_strt = sweep_strt_sh(sweep_strt_sh>1 & sweep_stop_sh<size(SWPs_sst.wbt_phase{zz,1},1));
                        sweep_stop = sweep_stop_sh(sweep_strt_sh>1 & sweep_stop_sh<size(SWPs_sst.wbt_phase{zz,1},1));
                    end
                    
                    N_spatial_bins = size(SWPs_sst.p_dec_flight{zz,1},1);
                    spt_bin_ids = [1:N_spatial_bins]';
                    
                    for ss=1:numel(sweep_strt)
                        
                        single_SWP.reg_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_flight{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
                        single_SWP.sft_posterior(counter) = {imgaussfilt(SWPs_sst.p_dec_shifted{zz,1}(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
                        single_SWP.rsp_posterior(counter) = {imresize(single_SWP.sft_posterior{counter,1},[size(single_SWP.sft_posterior{counter,1},1),n_reshape])};
                        
                        single_SWP.spk_dsty(counter) = {SWPs_sst.spk_dsty{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        single_SWP.rsz_spk_dsty(counter) = {interp1(single_SWP.spk_dsty{counter,1},linspace(1,numel(single_SWP.spk_dsty{counter,1}),n_reshape)')};
                        single_SWP.wbt_power(counter) = {SWPs_sst.wbt_power{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        single_SWP.LFP(counter) = {SWPs_sst.LFP{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        single_SWP.rsz_LFP(counter) = {interp1(single_SWP.LFP{counter,1},linspace(1,numel(single_SWP.LFP{counter,1}),n_reshape)')};
                        single_SWP.wbt(counter) = {SWPs_sst.wbt{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        single_SWP.rsz_wbt(counter) = {interp1(single_SWP.wbt{counter,1},linspace(1,numel(single_SWP.wbt{counter,1}),n_reshape)')};
                        single_SWP.fract_pos(counter) = {SWPs_sst.pos_real{zz,1}(sweep_strt(ss):sweep_stop(ss))/SWPs_sst.pos_real{zz,1}(end)};
                        
                        %single_SWP.clk(counter) = {SWPs_sst.clk{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        %single_SWP.rsz_clk(counter) = {interp1(single_SWP.clk{counter,1},linspace(1,numel(single_SWP.clk{counter,1}),n_reshape)')};
                        
                        single_SWP.wbt_phase(counter) = {SWPs_sst.wbt_phase{zz,1}(sweep_strt(ss):sweep_stop(ss))};
                        single_SWP.rsz_wbt_phase(counter) = {interp1(single_SWP.wbt_phase{counter,1},linspace(1,numel(single_SWP.wbt_phase{counter,1}),n_reshape)')};
                        
                        single_SWP.mean_spk_dsty(counter) = mean(single_SWP.spk_dsty{counter,1});
                        single_SWP.mean_wbt_power(counter) = mean(single_SWP.wbt_power{counter,1});
                        single_SWP.mean_fract_pos(counter) = mean(single_SWP.fract_pos{counter,1});
                        
                        [max_p,max_loc] = max(single_SWP.sft_posterior{counter,1},[],1);
                        %cnt_mass = spt_bin_ids'*single_SWP.sft_posterior{counter,1};
                        cnt_mass = max_loc;
                        single_SWP.med_jmp_distance(counter) = mean(abs(diff(cnt_mass)));
                        %                 pos_sprd = 0;
                        %                 for bb=1:numel(cnt_mass)
                        %                     pos_sprd = pos_sprd+sqrt((spt_bin_ids'-cnt_mass(bb)).^2*single_SWP.sft_posterior{counter,1}(:,bb));
                        %                 end
                        %single_SWP.avg_pos_spread(counter) = pos_sprd/numel(cnt_mass);
                        single_SWP.med_max_post(counter) = mean(max_p);
                        single_SWP.max_loc(counter) = {max_loc-N_spatial_bins/2};
                        single_SWP.cnt_mas(counter) = {cnt_mass-N_spatial_bins/2};
                        single_SWP.est_dist(counter) = {interp1((cnt_mass-N_spatial_bins/2)*bin_size_1D,linspace(1,numel(cnt_mass),n_reshape)')};
                        
                        counter=counter+1;
                    end
                    
                end
                warning('on');
                
                %=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
                SWP_sst = single_SWP(single_SWP.mean_fract_pos>0.15 & single_SWP.mean_fract_pos<0.85 & single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,10) & single_SWP.med_jmp_distance>0.1,:);
                N_sweeps = size(SWP_sst,1);
                
                %=== Calculate averages and STDs
                est_dist = zeros(size(SWP_sst.est_dist{1,1},1),size(SWP_sst,1));
                rsz_spk_dsty = zeros(size(SWP_sst.rsz_spk_dsty{1,1},1),size(SWP_sst,1));
                rsz_wbt = zeros(size(SWP_sst.rsz_wbt{1,1},1),size(SWP_sst,1));
                rsz_LFP = zeros(size(SWP_sst.rsz_LFP{1,1},1),size(SWP_sst,1));
                rsz_wbt_phase = zeros(size(SWP_sst.rsz_wbt{1,1},1),size(SWP_sst,1));
                for i=1:N_sweeps
                    est_dist(:,i) = SWP_sst.est_dist{i,1};
                    rsz_spk_dsty(:,i) = SWP_sst.rsz_spk_dsty{i,1};
                    rsz_wbt(:,i) = SWP_sst.rsz_wbt{i,1};
                    rsz_LFP(:,i) = SWP_sst.rsz_LFP{i,1};
                    rsz_wbt_phase(:,i) = SWP_sst.rsz_wbt_phase{i,1};
                end
                avg_est_dist = mean(est_dist,2);            sem_est_dist = std(est_dist,[],2)/sqrt(N_sweeps);
                avg_rsz_spk_dsty = mean(rsz_spk_dsty,2);    sem_rsz_spk_dsty = std(rsz_spk_dsty,[],2)/sqrt(N_sweeps);
                avg_rsz_wbt = mean(rsz_wbt,2);              sem_rsz_wbt = std(rsz_wbt,[],2)/sqrt(N_sweeps);
                avg_rsz_LFP = mean(rsz_LFP,2);              sem_rsz_LFP = std(rsz_LFP,[],2)/sqrt(N_sweeps);
                avg_rsz_wbt_phase = mean(rsz_wbt_phase,2);  sem_rsz_wbt_phase = std(rsz_wbt_phase,[],2)/sqrt(N_sweeps);
                
                %=== Calculate modulation score and averages
                SWP_mod_score(iii) = sum(abs(avg_est_dist-mean(avg_est_dist)));
                %SWP_mod_score(iii) = max(avg_est_dist);
                %SWP_mod_score(iii) = diff(prctile(avg_est_dist,[5 95]));
                avg_est_dist_sh(:,iii) = avg_est_dist;
                sem_est_dist_sh(:,iii) = sem_est_dist;
                
                %=== Show averages
                if iii == 1
                    rl_avg_est_dist = avg_est_dist;
                    rl_sem_est_dist = sem_est_dist;
                    figure('units','normalized','outerposition',[.3 .3 .45 .6],'Visible',fig_visibility);
                    tiledlayout(2,4,'TileSpacing','compact');
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_est_dist,avg_est_dist-sem_est_dist,avg_est_dist+sem_est_dist,'k');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    title('Average Decoding Error');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_spk_dsty,avg_rsz_spk_dsty-sem_rsz_spk_dsty,avg_rsz_spk_dsty+sem_rsz_spk_dsty,'r');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    title('Average Spike Density');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_LFP,avg_rsz_LFP-sem_rsz_LFP,avg_rsz_LFP+sem_rsz_LFP,'g');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    title('Average LFP');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt_phase,avg_rsz_wbt_phase-sem_rsz_wbt_phase,avg_rsz_wbt_phase+sem_rsz_wbt_phase,'m');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    title('Wingbeat Phase');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
                    nexttile;
                    plotWinterval_AF_v0(rsz_bin_ctrs,avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
                    hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
                    xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
                    sgtitle(['Average of ',num2str(N_sweeps),' cycles']);
                end
            end
            
            sh_avg_est_dist = mean(avg_est_dist_sh(:,2:end),2);
            sh_sem_est_dist = mean(sem_est_dist_sh(:,2:end),2);
            
            figure;
            plot_shuffled(SWP_mod_score);
            sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
            fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
            
            figure;
            plotWinterval_AF_v0(rsz_bin_ctrs,rl_avg_est_dist,rl_avg_est_dist-rl_sem_est_dist,rl_avg_est_dist+rl_sem_est_dist,'k'); hold on;
            plotWinterval_AF_v0(rsz_bin_ctrs,sh_avg_est_dist,sh_avg_est_dist-sh_sem_est_dist,sh_avg_est_dist+sh_sem_est_dist,'b');
            plot(0*[1 1],ylim,'k--'); xlim('tight');
            xticks([-180 0 180]);    title('Average Decoding Error');
            sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
            fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
            
        end
        
        %======================================================================================================================================================================
    end
    
end

%% SAVE THE DATA
for hide=1
    if options.savedata
        save([analysis_directory,'/Analyzed_NPs_', unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3},'.mat'],'SWP_table','LFP_PSD','NP_unit');
    end
end

