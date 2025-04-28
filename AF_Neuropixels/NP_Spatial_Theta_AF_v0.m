function NP_Spatial_Theta_AF_v0(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

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
    
    %=== Recluster Flights
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
    NP_unit_FS = NP_unit(NP_unit.fr>5,:);                   % Store high firing cells (putative interneurons) in a different variable (def 5 Hz)
    NP_unit = NP_unit(NP_unit.fr<5,:);                      % Exclude neurons with high-firing rate
    n_cells = size(NP_unit,1);                              % Number of units
    n_cells_FS = size(NP_unit_FS,1);                        % Number of units (putative interneurons)
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 10;                                             % Number of repetitions for shuffling (Place Cells)
    n_rep_mod_idx = 50;                                     % Number of repetitions for shuffling (Theta Sweeps)
    times_to_shift = t(randi([3*Fs T-3*Fs],1,n_rep));       % Time intervals for circshift
    smpls_to_shift = randi([3*Fs T-3*Fs],1,n_rep);          % Sample intervals for circshift
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
    t_Hi = NP_imu.t;                                        % Time for fast sampling (500 Hz, IMU and LFP)
    Fs_Hi = NP_imu.Fs;                                      % Sampling frequency for fast sampling (500 Hz, IMU and LFP)
    opt_int = best_prb_ch(2)+[-2:2];                        % Optimal channel interval for getting the LFP
    Tht_rng = [4 11];                                       % Theta range
    Swr_rng = [100 200];                                    % Sharp-wave-ripple range
    Dlt_rng = [1 4];                                        % Delta range
    Lwp_rng = [1 10];                                       % Low Pass Range (from Eliav et al., 2018)
    Sgm_rng = [30 60];                                      % Slow gamma 
    
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
    wbt_phase = angle(hilbert(a_flt_NP));   wbt_power = abs(hilbert(a_flt_NP)).^2;  wbt_freq = instfreq(a_flt_NP,Fs_Hi,'Method','hilbert');
    tht_phase = angle(hilbert(LFP_tht));    tht_power = abs(hilbert(LFP_tht)).^2;   tht_freq = instfreq(LFP_tht,Fs_Hi,'Method','hilbert');
    dlt_phase = angle(hilbert(LFP_dlt));    dlt_power = abs(hilbert(LFP_dlt)).^2;   dlt_freq = instfreq(LFP_dlt,Fs_Hi,'Method','hilbert');
    lwp_phase = angle(hilbert(LFP_lwp));    lwp_power = abs(hilbert(LFP_lwp)).^2;   lwp_freq = instfreq(LFP_lwp,Fs_Hi,'Method','hilbert');
    sgm_phase = angle(hilbert(LFP_sgm));    sgm_power = abs(hilbert(LFP_sgm)).^2;   sgm_freq = instfreq(LFP_sgm,Fs_Hi,'Method','hilbert');

    
    %=== Extract phase and power of Eliav's style LFP
    [~,min_locs] = findpeaks(-LFP_lwp,'MinPeakDistance',round(Fs_Hi*0.01)); % Find minima
    pi_tmp_phase = 2*pi*ones(numel(min_locs),1);                            % Assign phase 0 to minima
    pi_tmp_phase = cumsum(pi_tmp_phase);
    tmr_phase = interp1(min_locs,pi_tmp_phase,1:numel(t_Hi),'linear',0)';   % Interpolate
    tmr_phase = wrapTo2Pi(tmr_phase)-pi;                                    % Wrap to -pi:pi
    %plot(tmr_phase);    hold on;    plot(normalize(LFP_lwp));   plot(lwp_phase);
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
    
    %=== LFP and PSD
    LFP_r = LFP(~wBeats);      [PSD_r,f_PSD_r] = pwelch(LFP_r,[],[],[],Fs_Hi);  PSD_r = PSD_r./sum(PSD_r,1);  f_PSD_smpl_r = mean(diff(f_PSD_r)); % Rest
    LFP_f = LFP( wBeats);      [PSD_f,f_PSD_f] = pwelch(LFP_f,[],[],[],Fs_Hi);  PSD_f = PSD_f./sum(PSD_f,1);  f_PSD_smpl_f = mean(diff(f_PSD_f)); % Flight
    LFP_b = LFP( btheta);      [PSD_b,f_PSD_b] = pwelch(LFP_b,[],[],[],Fs_Hi);  PSD_b = PSD_b./sum(PSD_b,1);  f_PSD_smpl_b = mean(diff(f_PSD_b)); % Flight
    
    %=== PSD of the accelerometer signal during flight
    [PSD_a,f_PSD_a] = pwelch(a_abs_NP(wBeats),[],[],[],Fs_Hi);  PSD_a = PSD_a./sum(PSD_a,1);  f_PSD_smpl_a = mean(diff(f_PSD_a));
    
    %=== Average, smooth and normalize PSDs
    scaled_PSD_r = normalize(smoothdata(mean(PSD_r,2),'movmedian',1/f_PSD_smpl_r),'range');
    scaled_PSD_f = normalize(smoothdata(mean(PSD_f,2),'movmedian',1/f_PSD_smpl_f),'range');
    scaled_PSD_b = normalize(smoothdata(mean(PSD_b,2),'movmedian',1/f_PSD_smpl_b),'range');
    scaled_PSD_a = normalize(smoothdata(mean(PSD_a,2),'movmedian',1/f_PSD_smpl_a),'range');
    
    %=== Define pre-flight, flight and postflight time windows
    pre_flgt = false(size(wBeats)); for i=1:f_num_Hi,pre_flgt(f_smp_Hi(i,1)+[round(-2*Fs_Hi) :        0      ]) = 1;    end
    dur_flgt = false(size(wBeats)); for i=1:f_num_Hi,dur_flgt(f_smp_Hi(i,1)+[      0         : round(2*Fs_Hi)]) = 1;    end
    pst_flgt = false(size(wBeats)); for i=1:f_num_Hi,pst_flgt(f_smp_Hi(i,2)+[      0         : round(2*Fs_Hi)]) = 1;    end
    LFP_pre = LFP(pre_flgt);      [PSD_pre,f_PSD_pre] = pwelch(LFP_pre,[],[],[],Fs_Hi);  PSD_pre = PSD_pre./sum(PSD_pre,1);  f_PSD_smpl_pre = mean(diff(f_PSD_pre)); % Rest
    LFP_dur = LFP(dur_flgt);      [PSD_dur,f_PSD_dur] = pwelch(LFP_dur,[],[],[],Fs_Hi);  PSD_dur = PSD_dur./sum(PSD_dur,1);  f_PSD_smpl_dur = mean(diff(f_PSD_dur)); % Flight
    LFP_pst = LFP(pst_flgt);      [PSD_pst,f_PSD_pst] = pwelch(LFP_pst,[],[],[],Fs_Hi);  PSD_pst = PSD_pst./sum(PSD_pst,1);  f_PSD_smpl_pst = mean(diff(f_PSD_pst)); % Flight
    scaled_PSD_pre = normalize(smoothdata(mean(PSD_pre,2),'movmedian',1/f_PSD_smpl_pre),'range');
    scaled_PSD_dur = normalize(smoothdata(mean(PSD_dur,2),'movmedian',1/f_PSD_smpl_dur),'range');
    scaled_PSD_pst = normalize(smoothdata(mean(PSD_pst,2),'movmedian',1/f_PSD_smpl_pst),'range');
    
end

%% DEFINE SIMULATED SIGNAL WITH A GIVEN SPECTRUM AND RANDOM PHASE
while 0
    
    %=== Look at phase and amplitude of the low pass LFP
    N_lwp = length(LFP_lwp);              % Length of the signal
    Y_lwp = fft(LFP_lwp);                 % Compute the FFT of the signal
    f_lwp = (0:N_lwp-1)*(Fs_Hi/N_lwp);    % Frequency bins corresponding to FFT points
    
    % Calculate the amplitude spectrum (magnitude)
    amplitude = abs(Y_lwp) / N_lwp;      % Normalize by the number of points
    amplitude = amplitude(1:N_lwp/2+1);  % Only take positive frequencies
    amplitude(2:end-1) = 2*amplitude(2:end-1);  % Adjust for single-sided spectrum
    
    % Calculate the phase spectrum
    phase = angle(Y_lwp);
    phase_unwrapped = unwrap(phase);  % Optional: unwrap the phase
    
    % Frequency vector for positive frequencies only
    f_positive = f_lwp(1:N_lwp/2+1);
    phase_pos = phase(1:N_lwp/2+1);
    
    % Plot the amplitude spectrum
    figure('units','normalized','outerposition',[.1 .3 .4 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot(f_positive, smoothdata(amplitude,'movmean',100));
    xlabel('Frequency (Hz)');   ylabel('Amplitude');    title('Amplitude Spectrum of Signal y');    xlim([0 11]);
    nexttile;   plot(f_positive, phase(1:N_lwp/2+1));  % Phase for positive frequencies
    xlabel('Frequency (Hz)');   ylabel('Phase (radians)');  title('Phase Spectrum of Signal y');    xlim([0 11]);
    nexttile;   plot(f_positive,smoothdata(phase_pos,'movmean',1000))
    xlabel('Frequency (Hz)');   ylabel('Phase (radians)');  title('Phase Spectrum of Signal y');    xlim([0 11]);
    
    %     for i=1:10
    %       nexttile;  histogram(phase_pos(f_positive>i*1 & f_positive<(i+1)*1))
    %     end
    
    %=============================================================================================== SIMULATIONS
    %=== Generate a signal with random phase
    
    histogram(phase_pos(f_positive>1 & f_positive<10))
    
    N_sim = 2^16;  % Length of the simulated signal
    f_sim = (0:N_sim-1)*(Fs_Hi/N_sim);  % Frequency vector for the simulated signal
    
    % Generate random phase uniformly distributed between 0 and 2*pi for all frequencies
    random_phase = 2 * pi * rand(1, N_sim/2+1);
    %random_phase(f_sim(1:N_sim/2+1)<1 | f_sim(1:N_sim/2+1)>10) = 0;
    %random_phase = interp1(f_positive, smoothdata(phase_pos,'movmean',100), f_sim(1:N_sim/2+1), 'linear', 'extrap');
    amplitude_sim = interp1(f_positive, smoothdata(amplitude,'movmean',100), f_sim(1:N_sim/2+1), 'linear', 'extrap');
    Y_randphase = amplitude_sim .* exp(1i * random_phase);
    Y_randphase_full = [Y_randphase, conj(Y_randphase(end-1:-1:2))];
    
    % Perform inverse FFT to obtain the time-domain signal
    LFP_sim = ifft(Y_randphase_full, 'symmetric');
    t_sim = (0:N_sim-1)/Fs_Hi;
    
    %=== Plot the generated signal
    figure('units','normalized','outerposition',[.1 .3 .2 .3]);
    tiledlayout(2,2,'TileSpacing','tight');
    nexttile([1 2]);        plot(t_sim, LFP_sim);
    xlabel('Time (s)'); ylabel('Amplitude');    title('Generated Signal with Specified PSD');
    nexttile;   plot(f_positive, smoothdata(amplitude,'movmean',100));
    xlabel('Frequency (Hz)');   ylabel('Amplitude');    title('Amplitude Spectrum of Signal y (simulated)');    xlim([0 11]);
    nexttile;
    [pxx, f_sim] = periodogram(LFP_sim, [], N_sim, Fs_Hi);
    plot(f_sim, pxx);
    xlabel('Frequency (Hz)');   ylabel('Amplitude');    title('Amplitude Spectrum of Signal y');    xlim([0 11]);
    
    sim_phase = angle(hilbert(LFP_sim));    sim_power = abs(hilbert(LFP_sim)).^2;   sim_freq = instfreq(LFP_sim,Fs_Hi,'Method','hilbert');
    sim_phase = sim_phase(sim_power>prctile(sim_power,min_prctile_power));
    
    figure('units','normalized','outerposition',[.1 .3 .2 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    histogram([sim_phase;sim_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    nexttile;
    histogram([lwp_phase;lwp_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    
    
end

%% LFP PLOTTING
for hide=1
    
    %=== Plot the PSDs
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.4]);
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
    figure('units','normalized','outerposition',[0.1 0.3 0.3 0.4]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   plot(f_PSD_pre,scaled_PSD_pre,'LineWidth',2);    hold on;     plot(f_PSD_dur,scaled_PSD_dur,'LineWidth',2); plot(f_PSD_pst,scaled_PSD_pst,'LineWidth',2);
    rectangle('Position',[Tht_rng(1) 0 diff(Tht_rng) 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Pre)','LFP (Dur)','LFP (Post)');
    nexttile;   plot(f_PSD_pre,scaled_PSD_pre,'LineWidth',2);    hold on;     plot(f_PSD_dur,scaled_PSD_dur,'LineWidth',2); plot(f_PSD_pst,scaled_PSD_pst,'LineWidth',2);   plot(f_PSD_a,scaled_PSD_a,'LineWidth',2);
    rectangle('Position',[Tht_rng(1) 0 diff(Tht_rng) 1],'FaceColor',[0 0 0 0.2],'EdgeColor','none');
    xlim([0 20]);   xlabel('Frequency (Hz)');   ylabel('PSD (Norm.)');  legend('LFP (Pre)','LFP (Dur)','LFP (Post)','Accelerometer');
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
    figure('units','normalized','outerposition',[0.4 0.3 0.35 0.5]);
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
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(7,10,'TileSpacing','compact');
    int_t = [-2 4]; interval = [round(int_t(1)*Fs_Hi) : round(int_t(2)*Fs_Hi)];
    for i=1:min(f_num_Hi,70)
        zoom_int = f_smp_Hi(i,1)+interval;
        nexttile;   plot(t_Hi(zoom_int),normalize(LFP_tht(zoom_int)'),'r');
        hold on;    plot(t_Hi(zoom_int),normalize(a_abs_NP(zoom_int)')+5,'k');
        plot(t_Hi(zoom_int),2*click_trace(zoom_int)'-5,'b');
        plot(t_Hi(f_smp_Hi(i,1))*[1 1],ylim,'k--'); yticks([]);
        if i==min(f_num_Hi,70),xlabel('Time (s)');else,xticks([]); end
        if i==1,ylabel('Flight aligned LFP');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot theta bouts (70 max)
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(7,10,'TileSpacing','compact');
    for i=1:min(thtBt_num,70)
        zoom_int = thtBt_smp(i,1)+interval;
        nexttile;   plot(t_Hi(zoom_int),normalize(LFP_tht(zoom_int)'),'r');
        hold on;    plot(t_Hi(zoom_int),normalize(a_abs_NP(zoom_int)')+5,'k');
        plot(t_Hi(zoom_int),2*click_trace(zoom_int)'-5,'b');
        plot(t_Hi(thtBt_smp(i,1))*[1 1],ylim,'k--'); yticks([]);
        plot(t_Hi(thtBt_smp(i,2))*[1 1],ylim,'k--'); yticks([]);
        if i==min(thtBt_num,70),xlabel('Time (s)');else,xticks([]); end
        if i==1,ylabel('Theta Bouts');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Plot the accelerometer and the LFP
    figure('units','normalized','outerposition',[0 0.3 1 0.5]);
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
            NP_unitOnClus{1, j}(nc).map_interp = NP_unit(nc).f_clus(j).map_interp;
            
        end
    end
end

%% AUTOCORRELATION AND THETA MODULATION INDEX
for hide=1
    
    %=== Params
    bin_TMI = 0.01;
    max_lag_TMI = 0.5;
    
    %=== Autocorrelation and MI for principal cells
    for nc=1:n_cells
        
        %=== Initialize Autocorrelation, TMI and p value
        NP_unit(nc).AC = [];
        NP_unit(nc).AC_bins = [];
        NP_unit(nc).AC_PSD = [];
        NP_unit(nc).AC_PSD_freq = [];
        
        %=== Get empirical values
        [~,s_flight,~] = count_spikes_AF_v3(s{nc,1},t,[f_smp(:,1) f_smp(:,2)+1]);         % Keep only spikes emitted during flight
        [NP_unit(nc).AC,NP_unit(nc).AC_bins] = cross_correlogram_AF_v0(s_flight,s_flight,max_lag_TMI,bin_TMI);  % Autocorrelogram
        
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
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.7]);
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

%% PHASE LOCKING WITH LFP
for hide=1
    
    min_prctile_power = 25;
    n_bins_phase = 10;
    
    %=== Show phase locking to wingbeat
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        spk_phase2wbt = wbt_phase(s_smp_flight_hi{nc,1});   power_wbtAtspike =  wbt_power(s_smp_flight_hi{nc,1});  spk_phase2wbt = spk_phase2wbt(power_wbtAtspike>prctile(power_wbtAtspike,min_prctile_power));
        nexttile;
        histogram([spk_phase2wbt;spk_phase2wbt+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','k');
        yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'}); grid('on');
        %polarhistogram(spk_phase2wbt,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');  rticks([]);  thetaticks([]);
        if nc==1,title('Wingbeat Phase (Flight)');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Show phase locking to theta
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        spk_phase2tht = tht_phase(s_smp_flight_hi{nc,1});   power_thtAtspike =  tht_power(s_smp_flight_hi{nc,1});  spk_phase2tht = spk_phase2tht(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
        nexttile;
        histogram([spk_phase2tht;spk_phase2tht+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','r');
        yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'}); grid('on');
        %polarhistogram(spk_phase2tht,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');   rticks([]);  thetaticks([]);
        if nc==1,title('Theta Phase (Flight)');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Show phase locking to low pass
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        spk_phase2lwp = lwp_phase(s_smp_flight_hi{nc,1});   power_lwpAtspike =  lwp_power(s_smp_flight_hi{nc,1});  spk_phase2lwp = spk_phase2lwp(power_lwpAtspike>prctile(power_lwpAtspike,min_prctile_power));
        nexttile;
        histogram([spk_phase2lwp;spk_phase2lwp+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','b');
        yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'}); grid('on');
        %polarhistogram(spk_phase2lwp,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');   rticks([]);  thetaticks([]);
        if nc==1,title('Low Pass Phase (Flight)');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Show phase locking to Tamir style phase
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(10,20,'TileSpacing','none');
    for nc=1:min(n_cells,200)
        spk_phase2tmr = tmr_phase(s_smp_flight_hi{nc,1});   power_tmrAtspike =  tmr_power(s_smp_flight_hi{nc,1});  spk_phase2tmr = spk_phase2tmr(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
        nexttile;
        histogram([spk_phase2tmr;spk_phase2tmr+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','g');
        yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'}); grid('on');
        %polarhistogram(spk_phase2lwp,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');   rticks([]);  thetaticks([]);
        if nc==1,title('Tamir s Phase (Flight)');end
    end
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    n_bins_phase = 20;
    
    %=== Accumulate and show phase locking from all neurons (FLIGHT)
    spk_phase2wbt_all = []; spk_phase2tht_all = []; spk_phase2lwp_all = []; spk_phase2tmr_all = [];
    for nc=1:n_cells
        spk_phase2wbt = wbt_phase(s_smp_flight_hi{nc,1});   power_wbtAtspike =  wbt_power(s_smp_flight_hi{nc,1});  spk_phase2wbt = spk_phase2wbt(power_wbtAtspike>prctile(power_wbtAtspike,min_prctile_power));
        spk_phase2tht = tht_phase(s_smp_flight_hi{nc,1});   power_thtAtspike =  tht_power(s_smp_flight_hi{nc,1});  spk_phase2tht = spk_phase2tht(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
        spk_phase2lwp = lwp_phase(s_smp_flight_hi{nc,1});   power_lwpAtspike =  lwp_power(s_smp_flight_hi{nc,1});  spk_phase2lwp = spk_phase2lwp(power_lwpAtspike>prctile(power_lwpAtspike,min_prctile_power));
        spk_phase2tmr = tmr_phase(s_smp_flight_hi{nc,1});   power_tmrAtspike =  tmr_power(s_smp_flight_hi{nc,1});  spk_phase2tmr = spk_phase2tmr(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
        spk_phase2wbt_all = [spk_phase2wbt_all;  spk_phase2wbt];
        spk_phase2tht_all = [spk_phase2tht_all;  spk_phase2tht];
        spk_phase2lwp_all = [spk_phase2lwp_all;  spk_phase2lwp];
        spk_phase2tmr_all = [spk_phase2tmr_all;  spk_phase2tmr];
    end
    figure('units','normalized','outerposition',[.2 .2 0.3 0.5]);
    tiledlayout(2,4,'TileSpacing','none');
    nexttile;   polarhistogram(spk_phase2wbt_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');  title('Spikes in flight');
    nexttile;   polarhistogram(spk_phase2tht_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    nexttile;   polarhistogram(spk_phase2lwp_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');
    nexttile;   polarhistogram(spk_phase2tmr_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');
    nexttile;   histogram([spk_phase2wbt_all;spk_phase2wbt_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    nexttile;   histogram([spk_phase2tht_all;spk_phase2tht_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Theta');
    nexttile;   histogram([spk_phase2lwp_all;spk_phase2lwp_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Low pass');
    nexttile;   histogram([spk_phase2tmr_all;spk_phase2tmr_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Tamir');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Accumulate and show phase locking from all neurons (ALL)
    spk_phase2wbt_all = []; spk_phase2tht_all = []; spk_phase2lwp_all = []; spk_phase2tmr_all = [];
    for nc=1:n_cells
        spk_phase2wbt = wbt_phase(s_smp_hi{nc,1});   power_wbtAtspike =  wbt_power(s_smp_hi{nc,1});  spk_phase2wbt = spk_phase2wbt(power_wbtAtspike>prctile(power_wbtAtspike,min_prctile_power));
        spk_phase2tht = tht_phase(s_smp_hi{nc,1});   power_thtAtspike =  tht_power(s_smp_hi{nc,1});  spk_phase2tht = spk_phase2tht(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
        spk_phase2lwp = lwp_phase(s_smp_hi{nc,1});   power_lwpAtspike =  lwp_power(s_smp_hi{nc,1});  spk_phase2lwp = spk_phase2lwp(power_lwpAtspike>prctile(power_lwpAtspike,min_prctile_power));
        spk_phase2tmr = tmr_phase(s_smp_hi{nc,1});   power_tmrAtspike =  tmr_power(s_smp_hi{nc,1});  spk_phase2tmr = spk_phase2tmr(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
        spk_phase2wbt_all = [spk_phase2wbt_all;  spk_phase2wbt];
        spk_phase2tht_all = [spk_phase2tht_all;  spk_phase2tht];
        spk_phase2lwp_all = [spk_phase2lwp_all;  spk_phase2lwp];
        spk_phase2tmr_all = [spk_phase2tmr_all;  spk_phase2tmr];
    end
    figure('units','normalized','outerposition',[.2 .2 0.3 0.5]);
    tiledlayout(2,4,'TileSpacing','none');
    nexttile;   polarhistogram(spk_phase2wbt_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');  title('All Spikes');
    nexttile;   polarhistogram(spk_phase2tht_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    nexttile;   polarhistogram(spk_phase2lwp_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');
    nexttile;   polarhistogram(spk_phase2tmr_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');
    nexttile;   histogram([spk_phase2wbt_all;spk_phase2wbt_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    nexttile;   histogram([spk_phase2tht_all;spk_phase2tht_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Theta');
    nexttile;   histogram([spk_phase2lwp_all;spk_phase2lwp_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Low pass');
    nexttile;   histogram([spk_phase2tmr_all;spk_phase2tmr_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Tamir');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Accumulate and show phase locking from all neurons (CONTROL)
    spk_phase2wbt_all = []; spk_phase2tht_all = []; spk_phase2lwp_all = []; spk_phase2tmr_all = [];
    for nc=1:n_cells
        s_smp_hi_sh = randi([min(s_smp_hi{nc,1}),max(s_smp_hi{nc,1})],numel(s_smp_hi{nc,1}),1);
        spk_phase2wbt = wbt_phase(s_smp_hi_sh);  power_wbtAtspike =  wbt_power(s_smp_hi_sh);  spk_phase2wbt = spk_phase2wbt(power_wbtAtspike>prctile(power_wbtAtspike,min_prctile_power));
        spk_phase2tht = tht_phase(s_smp_hi_sh);  power_thtAtspike =  tht_power(s_smp_hi_sh);  spk_phase2tht = spk_phase2tht(power_thtAtspike>prctile(power_thtAtspike,min_prctile_power));
        spk_phase2lwp = lwp_phase(s_smp_hi_sh);  power_lwpAtspike =  lwp_power(s_smp_hi_sh);  spk_phase2lwp = spk_phase2lwp(power_lwpAtspike>prctile(power_lwpAtspike,min_prctile_power));
        spk_phase2tmr = tmr_phase(s_smp_hi_sh);  power_tmrAtspike =  tmr_power(s_smp_hi_sh);  spk_phase2tmr = spk_phase2tmr(power_tmrAtspike>prctile(power_tmrAtspike,min_prctile_power));
        spk_phase2wbt_all = [spk_phase2wbt_all;  spk_phase2wbt];
        spk_phase2tht_all = [spk_phase2tht_all;  spk_phase2tht];
        spk_phase2lwp_all = [spk_phase2lwp_all;  spk_phase2lwp];
        spk_phase2tmr_all = [spk_phase2tmr_all;  spk_phase2tmr];
    end
    figure('units','normalized','outerposition',[.2 .2 0.3 0.5]);
    tiledlayout(2,4,'TileSpacing','none');
    nexttile;   polarhistogram(spk_phase2wbt_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');  title('Control');
    nexttile;   polarhistogram(spk_phase2tht_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    nexttile;   polarhistogram(spk_phase2lwp_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');
    nexttile;   polarhistogram(spk_phase2tmr_all,linspace(-pi,pi,20),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');
    nexttile;   histogram([spk_phase2wbt_all;spk_phase2wbt_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Wingbeat');
    nexttile;   histogram([spk_phase2tht_all;spk_phase2tht_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Theta');
    nexttile;   histogram([spk_phase2lwp_all;spk_phase2lwp_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Low pass');
    nexttile;   histogram([spk_phase2tmr_all;spk_phase2tmr_all+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','g');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});   xlabel('Tamir');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Show phase distributions
    figure('units','normalized','outerposition',[.2 .2 .3 .3]);
    tiledlayout(1,4,'TileSpacing','none');
    nexttile;   histogram([wbt_phase;wbt_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','k');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    nexttile;   histogram([tht_phase;tht_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','r');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    nexttile;   histogram([lwp_phase;lwp_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','b');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    nexttile;   histogram([tmr_phase;tmr_phase+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'Normalization','probability','facealpha',.8,'edgecolor','none','FaceColor','g');    yticks([]); xticks(-pi:pi:3*pi);    xticklabels({'-180', '0', '180', '360', '540'});
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');   title('Phase distribution')
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% SPIKE TRIGGERED LFP
for hide=1
    
    interval = [round(-2*Fs_Hi):round(2*Fs_Hi)];
    
    %=== Spike triggered LFP (all spikes)
    figure('units','normalized','outerposition',[0 0 1 1]);
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
    figure('units','normalized','outerposition',[0 0 1 1]);
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
    
    %=== Accumulate (flight)
    spkTgdLfP_flight = [];
    spkTgdLfP_all = [];
    for nc=1:n_cells
        
        LFP_tmp = zeros(numel(interval),numel(s_smp_flight_hi{nc,1}));
        for i=1:numel(s_smp_flight_hi{nc,1})
            if all((s_smp_flight_hi{nc,1}(i)+interval)>1) && all((s_smp_flight_hi{nc,1}(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_flight_hi{nc,1}(i)+interval);
            end
        end
        spkTgdLfP_flight = [spkTgdLfP_flight,mean(LFP_tmp,2)];
        
        LFP_tmp = zeros(numel(interval),numel(s_smp_hi{nc,1}));
        for i=1:numel(s_smp_hi{nc,1})
            if all((s_smp_hi{nc,1}(i)+interval)>1) && all((s_smp_hi{nc,1}(i)+interval)<numel(t_Hi))
                LFP_tmp(:,i) = LFP(s_smp_hi{nc,1}(i)+interval);
            end
        end
        spkTgdLfP_all = [spkTgdLfP_all,mean(LFP_tmp,2)];
        
    end
    
    %=== Plot
    figure('units','normalized','outerposition',[.2 .1 .4 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    data_tmp = spkTgdLfP_flight';
    nexttile;
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('uV');
    data_tmp = spkTgdLfP_all';
    nexttile;
    plotWinterval_AF_v0(interval/Fs_Hi,mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    hold on;    plot([0 0],ylim,'k--'); xlabel('Time (s)'); ylabel('uV');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% SHOW EXAMPLE PLACE CELLS AND LOOK AT PHASE PRECESSION AND AUTOCORRELATIONS
for hide=1
    
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
        
        find_optimal_phase = 1;
        
        %=== Get the spike samples for this cluster
        s_smp_flight_lo_clus = cell(n_place_cells,1);% Spike samples during the flight cluster (sorted), low sampling
        s_smp_flight_hi_clus = cell(n_place_cells,1);% Spike samples during the flight cluster (sorted), hi sampling
        for nc = 1:n_place_cells
            [~,~,valid_samples] = count_spikes_AF_v2(s{sorted_tko_plcCells(nc),1},t_Hi,[smp1_clus_hi smp2_clus_hi]);
            s_smp_flight_hi_clus{nc,1} = s_smp_hi{sorted_tko_plcCells(nc),1}(valid_samples);
            s_smp_flight_lo_clus{nc,1} = round(m_slope * s_smp_flight_hi_clus{nc,1}+q_intcp);
        end
        
        %=== Initialize the cells containing, for each cell and spike the...
        spk2plc_dist = s_smp_flight_lo_clus;    %... distance to place field
        spk2wbt_phase = s_smp_flight_lo_clus;   %... wingbeat phase of the spike
        spk2tht_phase = s_smp_flight_lo_clus;   %... theta phase of the spike
        spk2lwp_phase = s_smp_flight_lo_clus;   %... lowpass phase of the spike
        spk2tmr_phase = s_smp_flight_lo_clus;   %... lowpass phase of the spike
        
        phase_wbt_place_corr = zeros(n_place_cells,1);
        phase_tht_place_corr = zeros(n_place_cells,1);
        phase_lwp_place_corr = zeros(n_place_cells,1);
        phase_tmr_place_corr = zeros(n_place_cells,1);
        
        %=== Extract a few variables for phase precession (NEW)
        for nc = 1:n_place_cells
            
            %d2place = lin_flight - NP_subtable_sorted.phase_max(nc)*0.01;
            d2place = len_flight;

            
            spk2plc_dist{nc,1} = d2place(s_smp_flight_lo_clus{nc,1});
            spk2wbt_phase{nc,1} = wbt_phase(s_smp_flight_hi_clus{nc,1});
            spk2tht_phase{nc,1} = tht_phase(s_smp_flight_hi_clus{nc,1});
            spk2lwp_phase{nc,1} = lwp_phase(s_smp_flight_hi_clus{nc,1});
            spk2tmr_phase{nc,1} = tmr_phase(s_smp_flight_hi_clus{nc,1});
            
            if find_optimal_phase
                
                delta_phase = linspace(-pi,pi,20);
                pp_corr = NaN(numel(delta_phase),1);
                
                %=== Wingbeat
                temp_phase = wbt_phase(s_smp_flight_hi_clus{nc,1});
                for dtp = 1:numel(delta_phase)
                    pp_corr(dtp) = corr(spk2plc_dist{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                end
                [phase_wbt_place_corr(nc,1),min_dtp] = min(pp_corr);
                spk2wbt_phase{nc,1} = wrapToPi(temp_phase+delta_phase(min_dtp));
                
                %=== Theta
                temp_phase = tht_phase(s_smp_flight_hi_clus{nc,1});
                for dtp = 1:numel(delta_phase)
                    pp_corr(dtp) = corr(spk2plc_dist{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                end
                [phase_tht_place_corr(nc,1),min_dtp] = min(pp_corr);
                spk2tht_phase{nc,1} = wrapToPi(temp_phase+delta_phase(min_dtp));
                
                %=== Lowpass
                temp_phase = lwp_phase(s_smp_flight_hi_clus{nc,1});
                for dtp = 1:numel(delta_phase)
                    pp_corr(dtp) = corr(spk2plc_dist{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                end
                [phase_lwp_place_corr(nc,1),min_dtp] = min(pp_corr);
                spk2lwp_phase{nc,1} = wrapToPi(temp_phase+delta_phase(min_dtp));
                
                %=== Tamir
                temp_phase = tmr_phase(s_smp_flight_hi_clus{nc,1});
                for dtp = 1:numel(delta_phase)
                    pp_corr(dtp) = corr(spk2plc_dist{nc,1},wrapToPi(temp_phase+delta_phase(dtp)));
                end
                [phase_tmr_place_corr(nc,1),min_dtp] = min(pp_corr);
                spk2tmr_phase{nc,1} = wrapToPi(temp_phase+delta_phase(min_dtp));
                
            end
            
        end
        
        %=== Visualize cluster
        figure('units','normalized','outerposition',[0.35 .3 0.2 .4]);
        plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   zlim(r_lim(3,:));  view(0,90);
        xlabel(['x, Cluster' num2str(j) ' (' num2str(n_fclus) ' flights),']);    ylabel('y');    hold on;
        for ii=1:n_fclus
            plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(j,:));
        end
        textscatter(avg_takeoff(1),avg_takeoff(2),"Take-off");     hold off;   axis equal;
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
        figure('units','normalized','outerposition',[0.6 .3 0.3 .4]);
        tiledlayout(1,2,'TileSpacing','tight');
        nexttile;   imagesc([0 size(a_abs_all,1)/Fs_Hi],[],a_abs_all');    colormap(gray); 
        xlabel('Time (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
        nexttile;   imagesc([0 size(a_abs_all_sh,1)/Fs_Hi],[],a_abs_all_sh');    colormap(gray); 
        xlabel('Time (s)');  ylabel('Sorted Flight #'); title('Accelerometer');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
        figure('units','normalized','outerposition',[0 0 1 1]);
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
            scatter([spk2plc_dist{ncc,1};spk2plc_dist{ncc,1}],[spk2wbt_phase{ncc,1};spk2wbt_phase{ncc,1}+2*pi],6,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            
            
            
            
            xlim(prctile(spk2plc_dist{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Wingbeat Phase histogram)
        n_bins_phase = 10;
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
        figure('units','normalized','outerposition',[0 0 1 1]);
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
            scatter([spk2plc_dist{ncc,1};spk2plc_dist{ncc,1}],[spk2tht_phase{ncc,1};spk2tht_phase{ncc,1}+2*pi],6,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Theta Phase histogram)
        n_bins_phase = 10;
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
        
        %=== Plot place cells around the selected cluster (+ Lowpass Phase)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
            scatter([spk2plc_dist{ncc,1};spk2plc_dist{ncc,1}],[spk2lwp_phase{ncc,1};spk2lwp_phase{ncc,1}+2*pi],6,'filled','MarkerFaceColor','b','MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Lowpass Phase histogram)
        n_bins_phase = 10;
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
            histogram([spk2lwp_phase{ncc,1};spk2lwp_phase{ncc,1}+2*pi],unique([linspace(-pi,pi,n_bins_phase),linspace(pi,3*pi,n_bins_phase)]),'facealpha',.5,'edgecolor','none','FaceColor','b');
            xticks(pi*[-1 0 1 2 3]);  xticklabels({'-180', '0', '180', '360', '540'});
            xlabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Tamir Phase)
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
            scatter([spk2plc_dist{ncc,1};spk2plc_dist{ncc,1}],[spk2tmr_phase{ncc,1};spk2tmr_phase{ncc,1}+2*pi],6,'filled','MarkerFaceColor',[0 .6 .3],'MarkerFaceAlpha',0.5);    hold on;
            plot([fld_ctr fld_ctr],ylim,'k--'); yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});
            xlim(prctile(spk2plc_dist{ncc,1},[1 99])*1);
            %xlim([-0.5 0.5]);
            axis square;
            xlabel('Traj. Lenght (m)');   ylabel('Spike Phase');
            
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells around the selected cluster (+ Tamir Phase histogram)
        n_bins_phase = 10;
        n2show = min(n_place_cells,60);
        figure('units','normalized','outerposition',[0 0 1 1]);
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
        figure('units','normalized','outerposition',[0 0 1 1]);
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
        figure('units','normalized','outerposition',[.3 .1 .1 .35]);
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
        figure('units','normalized','outerposition',[0 .2 .8 .3]);
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
        figure('units','normalized','outerposition',[.5 .3 .1 .4]);
        Snake = NP_subtable.map_interp;         Snake = Snake./max(Snake,[],2);
        [~,order] = sort(NP_subtable.phase_max,'descend');
        imagesc(Snake(order,:));    colormap(viridis);  ylabel('1D Field');   xlabel('Flight Phase (%)');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
    end
end

