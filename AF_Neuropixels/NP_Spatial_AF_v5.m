function NP_Spatial_AF_v5(folder_name)
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)
for hide=1
    %==========================%
    %=== RELEVANT VARIABLES ===%
    %==========================%
    
    
end

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    nrnFile = dir('SU_kilosort*');                                                      % Single Units
    unique_ID = options.unique_ID;                                                      % Session identifier
    load('LFP_probe2.mat');                                                             % This file contains LFP data from probe 1
    load('RPL_probe2.mat');                                                             % This file contains Ripple data from probe 1
    NP_unit = table();                  MUA_unit = table();                             % Single Units and MUA structures
    n_probes = size(nrnFile,1);                                                         % Number of probes
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        NP_unit = [NP_unit; out.good_units];                                            % Single units and MUA are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;                                                                          % Clear unnecessary variables
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);             % Add average firing rate
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);           % Add average firing rate
    
    %=== Create analysis folder for storing the results
    options.savedata = 0;                                   % If creating folder for saving the data or not
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
    NP_unit_FS = NP_unit(NP_unit.fr>1,:);                   % Store high firing cells (putative interneurons) in a different variable
    NP_unit = NP_unit(NP_unit.fr<1,:);                      % Exclude neurons with high-firing rate
    n_cells = size(NP_unit,1);                              % Number of units
    n_cells_FS = size(NP_unit_FS,1);                        % Number of units (putative interneurons)
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 100;                                            % Number of repetitions for shuffling (Place Cells)
    times_to_shift = t(randi([0.5*Fs T-0.5*Fs],1,n_rep));   % Time intervals for circshift
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
    n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
    col_clus = hsv(n_surv_clusters);                        % Colors for clusters
    min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
    t_Hi = NP_imu.t;                                        % Time for fast sampling (500 Hz, IMU and LFP)
    Fs_Hi = NP_imu.Fs;                                      % Sampling frequency for fast sampling (500 Hz, IMU and LFP)
    
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
    s_FS = cell(n_cells,n_rep);               % Single units' spikes
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
    LFP_rp = bandpass(RPL_out.LFP_trace,[100 200],RPL_out.Fs);                                           % Ripple trace
    LFP_th_power = smoothdata(zscore(abs(hilbert(LFP_th)),[],1),1,'gaussian',0.05*Fs_Hi);                % Theta power
    
    %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                         % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_Hi);              % Filtered at the wing-beat frequency
    [up,lo] = envelope(a_flt_NP,round(0.06*Fs_Hi),'peak');  % Upper and lower envelopes
    env = normalize(up - lo,'range');                           % Amplitude of the envelope
    env_th = otsuthresh(histcounts(env));                       % Threshold (based on Otsu method). Can be set at 0.35
    wBeats = movsum(env>env_th,2*Fs_Hi)>Fs_Hi/5;        % Euristic criterion for flight detection
    
    %=== Calculate Population Spike Density
    all_s = sort(vertcat(s{:,1}));
    all_rate = kernel_rate_AF_v1(all_s,0.1,t)/n_cells;
    
    %=== Calculate Population Spike Density for putative interneurons
    all_s_FS = sort(vertcat(s_FS{:,1}));
    all_rate_FS = kernel_rate_AF_v1(all_s_FS,0.1,t)/n_cells_FS;
    
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
            [~,locs,width] = findpeaks(NP_unit(nc).f_clus(j).field,linspace(0,mean(fc_lgt),numel(NP_unit(nc).f_clus(j).field)),'MinPeakDistance',1,'MinPeakHeight',0.5*NP_unit(nc).f_clus(j).peakHz,'SortStr','descend');
            NP_unit(nc).f_clus(j).n_fields = numel(locs);
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
            NP_unitOnClus{1, j}(nc).n_fields = NP_unit(nc).f_clus(j).n_fields;
            NP_unitOnClus{1, j}(nc).f_width = NP_unit(nc).f_clus(j).f_width;
            NP_unitOnClus{1, j}(nc).map_interp = NP_unit(nc).f_clus(j).map_interp;
            NP_unitOnClus{1, j}(nc).phase_max = NP_unit(nc).f_clus(j).phase_max;
            
        end
    end
end

%% LOOK AT POPULATION RASTER PLOT
for hide=1
    
    %=== Plot raster with sorted units (by firing frequency), LFP and Spike Density
    figure('units','normalized','outerposition',[0 .3 1 .6]);
    tiledlayout(11,1,'TileSpacing','none');
    cx(1) = nexttile(1,[4 1]);
    for nc= 1:n_cells
        plot(s{rate_sorted_ids(nc)}, nc*ones(size(s{rate_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(5,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(6,[2 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_th);  xticks([]); legend('All','Theta');  ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
    cx(4) = nexttile(8,[2 1]);
    LFP_lw = interp1(t_Hi,LFP_mn,t,'linear',median(LFP_mn));
    [PS_LFP,freq_SG] = cwt(LFP_lw,Fs);
    imagesc([t(1),t(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[2 10])); shading interp;  colormap(hot);
    set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 50]);  yticks([1 5 10 20 50]);    ylabel('Freq (Hz)');
    cx(5) = nexttile(10,[2 1]);
    plot(t,all_rate,'k');   ylabel('Spike Density');    ylim([0 prctile(all_rate,99.99)]);   set(gca,'TickLength',[0 0]);
    
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% LOOK AT SLOW-OSCILLATIONS
for hide=1
    
    %=== Find periods of high VS LOW activity
    [hiRate,~,hiRate_str,hiRate_stp,~]= BiLevel_Segm_AF_v0(all_rate,3*mad(all_rate,1),Fs,0.25,0.5);
    s_hiRate = s(:,1);
    for nc=1:n_cells
        [~,s_hiRate{nc,1}] = count_spikes_AF_v1(s{nc,1},t,[hiRate_str hiRate_stp+1]);
    end
    
    %=== Try to run PCA on Rate during rest
    Rate_rest = Rate;
    Rate_rest(bflying,:) = 0;
    Rate_rest(~hiRate,:) = 0;
    [coeff,score,latent] = pca(Rate_rest);
    angle_rate = cart2pol(coeff(:,1),coeff(:,2));
    [sorted_angle,correlation_sorted_ids] = sort(angle_rate);
    
    %=== Plot raster with sorted units (by firing frequency), LFP and Spike Density
    figure('units','normalized','outerposition',[0 .3 1 .6]);
    tiledlayout(14,1,'TileSpacing','none');
    cx = [];
    cx(1) = nexttile(1,[5 1]);
    for nc= 1:n_cells
        plot(s{correlation_sorted_ids(nc)}, nc*ones(size(s{correlation_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(6,[5 1]);
    for nc= 1:n_cells
        plot(s_hiRate{correlation_sorted_ids(nc)}, nc*ones(size(s_hiRate{correlation_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(11,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);   xticks([]);
    cx(4) = nexttile(12,[2 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_th);  legend('All','Theta');  ylabel('LFP (norm)');  set(gca,'TickLength',[0 0]);
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Calculate ISI within the indicated range and average across cells
    figure('units','normalized','outerposition',[.2 .3 .25 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    Cross_c = [];
    for nc=1:n_cells
        [cross_c,ctrs_theta] = cross_correlogram_AF_v0(s_hiRate{nc,1},s_hiRate{nc,1},10,0.05);
        Cross_c = [Cross_c;cross_c'];
    end
    nexttile;   area(ctrs_theta,smoothdata(nanmean(Cross_c,1),'movmean',1)); set(gca,'YScale', 'log');  xlabel('Time (s)'); ylabel('Average Autocorrelation');
    nexttile;   histogram(diff(t(hiRate_str)),[0:0.5:30]);  xlabel('Time between High Activity epochs');    ylabel('Count');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% LOOK AT SLOW-OSCILLATIONS (INTERNEURONS)
while 0
    
    %=== Find periods of high VS LOW activity
    [hiRate,~,hiRate_str,hiRate_stp,~]= BiLevel_Segm_AF_v0(all_rate_FS,5*mad(all_rate_FS,1),Fs,0.25,0.5);
    s_hiRate = s_FS(:,1);
    for nc=1:n_cells_FS
        [~,s_hiRate{nc,1}] = count_spikes_AF_v1(s_FS{nc,1},t,[hiRate_str hiRate_stp+1]);
    end
    
    %=== Try to run PCA on Rate during rest
    Rate_rest = Rate;
    Rate_rest(bflying,:) = 0;
    Rate_rest(~hiRate,:) = 0;
    [coeff,score,latent] = pca(Rate_rest);
    angle_rate = cart2pol(coeff(:,1),coeff(:,2));
    [sorted_angle,correlation_sorted_ids] = sort(angle_rate);
    
    %=== Plot raster with sorted units (by firing frequency), LFP and Spike Density
    figure('units','normalized','outerposition',[0 .3 1 .6]);
    tiledlayout(14,1,'TileSpacing','none');
    cx = [];
    cx(1) = nexttile(1,[5 1]);
    for nc= 1:n_cells_FS
        plot(s_FS{correlation_sorted_ids(nc)}, nc*ones(size(s_FS{correlation_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells_FS*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells_FS,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells_FS]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(6,[5 1]);
    for nc= 1:n_cells_FS
        plot(s_hiRate{correlation_sorted_ids(nc)}, nc*ones(size(s_hiRate{correlation_sorted_ids(nc)})), 'k|','MarkerSize', max(1,round(n_cells_FS*0.01)));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells_FS,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells_FS]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(11,[1 1]);  plot(t_Hi,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);   xticks([]);
    cx(4) = nexttile(12,[2 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_th);  legend('All','Theta');  ylabel('LFP (norm)');  set(gca,'TickLength',[0 0]);
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
    %=== Calculate ISI within the indicated range and average across cells
    figure('units','normalized','outerposition',[.2 .3 .25 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    Cross_c = [];
    for nc=1:n_cells_FS
        if NP_unit_FS(nc).fr<10
            nc
            %Cross_c = [Cross_c;numel(s_hiRate{nc,1})];
            [cross_c,ctrs_theta] = cross_correlogram_AF_v0(s_hiRate{nc,1},s_hiRate{nc,1},5,0.010);
            Cross_c = [Cross_c;cross_c'];
        end
    end
    nexttile;   area(ctrs_theta,smoothdata(nanmean(Cross_c,1),'movmean',1)); set(gca,'YScale', 'log');  xlabel('Time (s)'); ylabel('Average Autocorrelation');
    nexttile;   histogram(diff(t(hiRate_str)),[0:0.5:30]);  xlabel('Time between High Activity epochs');    ylabel('Count');
    sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
    fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
    
end

%% SHOW EXAMPLE PLACE CELLS AND LOOK AT REPLAY
for hide=1
    
    for j = setdiff(imp_bat_clusters,1)
        
        %=== Get the subtable of cells for template matching and cells for decoding
        NP_table = struct2table(NP_unitOnClus{1, j});
        
        %=== Some example conditions
        %NP_table.place_cond = NP_table.spkPerflight>1 & NP_table.peakHz>3 & NP_table.corr_m>0.4;    % DEFAULT
        %NP_table.place_cond = NP_table.peakHz>1 & NP_table.corr_m>0.3 & NP_table.n_fields == 1 & NP_table.f_width<0.5;
        %NP_table.place_cond = NP_table.peakHz>1 & NP_table.corr_m>0 & NP_table.n_fields == 1 & NP_table.f_width<0.5;
        
        %=== For template matching
        NP_table.place_cond = NP_table.spkPerflight>1 & NP_table.peakHz>3 & NP_table.stab_m>0.4;
        NP_subtable = NP_table(NP_table.place_cond,:);
        n_place_cells = size(NP_subtable,1);
        
        %=== For decoding
        NP_table.use4decoding = NP_table.peakHz>0;
        NP_subtable_dec = NP_table(NP_table.use4decoding,:);
        n_dec_cells = size(NP_subtable_dec,1);
        
        %=== Extract some variables for plotting
        smp1_clus = f_clus.strt_frame(f_clus.id == j)';                     % Takeoff sample
        smp2_clus = f_clus.stop_frame(f_clus.id == j)';                     % Landing sample
        cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
        smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
        n_fclus = numel(smp1_clus);                                         % Number of surviving flights in the cluster
        id = find(f_clus.id==j);                                            % ids of the flight clusters                                            
        avg_takeoff = mean(squeeze(f_clus.ds_pos(:,1,id)),2);               % Average takeoff position
        avg_landing = mean(squeeze(f_clus.ds_pos(:,end,id)),2);             % Average landing position
        [~,sorted_pk] = sort([NP_subtable.field_loc],'descend');            % Target sequence of activation (from first to last in time)
        d_fclus = mean(t(smp2_clus)-t(smp1_clus));                          % Mean duration of the flight
        all_s_plc = sort(vertcat(s{NP_subtable.cell_id,1}));                % Define rate by pooling spikes from place cells
        all_rate_plc = kernel_rate_AF_v1(all_s_plc,0.1,t)/n_place_cells;
        t_ic = [-1 1];                                                      % Time for plotting around the flight
        max_Replay_dur = 0.5;                                               % Maximal duration for Replay events
        Rate_treshold_SD = 2;                                               % Can try to decrease this one
        t_bin_dur = 0.010;                                                  % Time bin duration for decoding
        t_bin_ovl = 0.005;                                                  % Time bin overlap for decoding
        smooth_f = [.1 1];                                                  % [space bin,time bin]
        
        %=== Get the spike times happening around flight (plus t_ic) and calculate for each cell median time spike to takeoff
        s_flight1 = cell(n_cells,1);
        for nc = 1:n_cells,[~,s_flight1{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);end
        t2tko = NaN(n_cells,1);
        for nc =1:n_cells
            tgt_seq = knnsearch(s_flight1{nc},t(smp1_clus));
            if ~isempty(tgt_seq)
                t2tko(nc) = median(t(smp1_clus)-s_flight1{nc}(tgt_seq));
            end
        end
        [~,sorted_tko] = sort(t2tko,'descend');
        sorted_tko_plcCells = intersect(sorted_tko,NP_subtable.cell_id,'stable');   % Restrict to place cells only
        
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
        
        %=== Plot place cells around the selected cluster
        %NP_subtable = sortrows(NP_subtable,'f_width','ascend');
        n2show = min(n_place_cells,50);
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(5,4*ceil(n2show/5),'TileSpacing','tight');
        sgtitle('Example place cells');
        for ncc=1:n2show
            
            %=== Get the cell ID
            nc = NP_subtable.cell_id(ncc);
            
            %=== Get the spikes happening during cluster j
            [~,s_flight] = count_spikes_AF_v1(s{nc,1},t,[smp1_clus smp2_clus+1]);
            
            %=== Calculate linearized trajectories and average 3d path
            flight_pos = f_clus.pos(:,:,f_clus.id==j);
            flight_pos_rsh = reshape(flight_pos,3,[]);
            flight_pos_rsh = flight_pos_rsh(:,all(~isnan(flight_pos_rsh),1));
            
            %=== Find the closest positional sample to each spike (plotting purposes)
            if ~isempty(s_flight)
                k = round((s_flight-t(1))*Fs);    k = k(k<T);
                spikes_pos = r(k,:);
            else, spikes_pos = [nan nan nan];end
            
            %=== Plot
            nexttile([1 2]);   plot3(flight_pos_rsh(1,:),flight_pos_rsh(2,:),0*flight_pos_rsh(3,:),'.','MarkerSize',1,'MarkerEdgeColor',.9*[1 1 1]);
            hold on;    textscatter(flight_pos_rsh(1,1)-0.05,flight_pos_rsh(2,1),"**");
            scatter3(spikes_pos(:,1),spikes_pos(:,2),spikes_pos(:,3),6,'MarkerFaceColor', 'r','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
            hold off;    axis equal;
            view(0,90); xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
            xticks([]); yticks([]);
            nexttile([1 2]);   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],.5,[],1);
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot place cells sorted by peak location
        figure('units','normalized','outerposition',[.3 .1 .1 .35]);
        tiledlayout(n_place_cells,1,'TileSpacing','none');
        for ncc=1:n_place_cells
            nc = NP_subtable.cell_id((sorted_pk(ncc)));
            nexttile;   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],0.1,[],1);
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
        
        %=== REPLAY ANALYSIS STARTS HERE
        %=== Detect candidate events from spike density from place cells only
        trace_for_candidate = all_rate_plc;
        trace_for_candidate(bflying) = 0;
        [candidate_event_vector,~,candidate_event_str,candidate_event_stp,candidate_event_width] = BiLevel_Segm_AF_v0(zscore(trace_for_candidate)',Rate_treshold_SD,Fs,0.1,0.05);
        t_temp= t(candidate_event_str);
        candidate_event_str =     candidate_event_str(candidate_event_width<max_Replay_dur & t_temp>0);
        candidate_event_stp =     candidate_event_stp(candidate_event_width<max_Replay_dur & t_temp>0);
        candidate_event_width = candidate_event_width(candidate_event_width<max_Replay_dur & t_temp>0);
        
        %=== Plot sorted activity (entire session)
        figure('units','normalized','outerposition',[0 .2 1 .5]);
        tiledlayout(6,1,'TileSpacing','none');
        bx(1) = nexttile([2 1]);
        for nc= 1:n_place_cells
            plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'k|','MarkerSize', 5);  hold on;             % Raster for each cell
        end
        area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
        plot((t(f_clus.strt_frame(f_clus.id == j)')*[1 1])',(ones(size(f_clus.strt_frame(f_clus.id == j)'))*[0 n_cells])','r');
        ylim([0 n_place_cells]);  ylabel('Unit #');    xticks([]);
        bx(2) = nexttile([1 1]);  plot(t,all_rate_plc,'k');     ylim([0 prctile(all_rate_plc,99.99)]);   set(gca,'TickLength',[0 0]);
        bx(3) = nexttile([1 1]);  area(t,candidate_event_vector);                                        set(gca,'TickLength',[0 0]);
        bx(4) = nexttile([1 1]);  plot(RPL_out.time_vector,RPL_out.LFP_trace,'b');                       set(gca,'TickLength',[0 0]);
        bx(5) = nexttile([1 1]);   plot(t_Hi,LFP_mn); hold on;    plot(t_Hi,LFP_th);  legend('All','Theta');  ylabel('LFP (norm)');
        ylabel('Ave Population firing rate');   xlabel('Time (s)'); linkaxes(bx,'x');   xlim('tight');
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Find Candidate Replays
        cRt = [t(candidate_event_str),t(candidate_event_stp)];
        n_candidate_events = size(cRt,1);
        candidate_perc_act = zeros(n_candidate_events,1);
        candidate_rkcorr_c = zeros(n_candidate_events,1);
        candidate_rkcorr_p = NaN(n_candidate_events,1);
        candidate_rkcorr_shuffles = NaN(n_candidate_events,100);
        seq_min_spikes = 3;         % Minimum number of spikes in the sequence
        seq_min_fr_active = 0.1;    % Minimum fraction of active cells in the sequence
        
        %=== Assess event quality
        warning('off');
        for i=1:n_candidate_events
            tL = cRt(i,1);
            tR = cRt(i,2);
            candidate_s = cellfun(@(x) x(x>tL & x<tR), s(sorted_tko_plcCells,1),'UniformOutput', false);   % Extract the spikes happening during the event
            candidate_perc_act(i) =1-sum(cellfun(@isempty,candidate_s))/n_place_cells;
            if  numel(vertcat(candidate_s{:}))>seq_min_spikes && candidate_perc_act(i)>seq_min_fr_active
                [candidate_rkcorr_c(i),candidate_rkcorr_p(i),candidate_rkcorr_shuffles(i,:)] = spikeSeq_analysis_AF_v0(candidate_s,100);
            end
        end
        warning('on');
        
        %=== Define Replay Table
        Rpl_table = table();
        Rpl_table.tL = cRt(:,1);    Rpl_table.tR = cRt(:,2);                                            % Start Stop times
        Rpl_table.smpL = knnsearch(t,Rpl_table.tL); Rpl_table.smpR = knnsearch(t,Rpl_table.tR);         % Start Stop samples
        Rpl_table.f = candidate_perc_act;                                                               % Fraction of active cells
        Rpl_table.n = candidate_perc_act*n_place_cells;                                                 % Number of active cells
        Rpl_table.C = candidate_rkcorr_c;                                                               % Rank correlation
        Rpl_table.p = candidate_rkcorr_p;                                                               % p value
        Rpl_table.shuffles = mat2cell(candidate_rkcorr_shuffles,ones(1,numel(candidate_rkcorr_p)),100); % Shuffled Rank correlation
        Rpl_table.d1 = vecnorm(r(Rpl_table.smpL,:)-avg_takeoff',2,2);                                   % Distance from takeoff spot
        Rpl_table.d2 = vecnorm(r(Rpl_table.smpL,:)-avg_landing',2,2);                                   % Distance from landing spot
        [~,Rpl_table.t2tko] = knnsearch(t(smp1_clus),Rpl_table.tL);                                     % Time from nearest takeoff
        [~,Rpl_table.t2lnd] = knnsearch(t(smp2_clus),Rpl_table.tL);                                     % Time from nearest landing

        %=== Define subset of interest
        cRT_good = find(Rpl_table.n>4 & Rpl_table.p<0.05);
        n_event2show = min(100,size(cRT_good,1));
        cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
        
        %=== Plot a few features for all events, good and shuffled (when available)
        figure('units','normalized','outerposition',[.2 .3 .5 .35]);
        tiledlayout(1,5,'TileSpacing','tight');
        cRT_else = setdiff(1:size(Rpl_table,1),cRT_good);
        maxD = max([Rpl_table.d1;Rpl_table.d2]);
        C_edges = [-1:0.2:1];
        D_edges =[0,10.^[-.5:0.1:log10(maxD)]];
        T_edges = [0,10.^[0:0.1:log10(3600)]];
        nexttile;   
        histogram(Rpl_table.C(cRT_good),C_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table.C(cRT_else),C_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');   
        histogram(horzcat(Rpl_table.shuffles{:}),C_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','k');   hold off;
        ylabel('Fraction'); xlabel('Rank Correlation'); legend('Good','Excluded','Shuffled','Location','northoutside');
        nexttile;
        histogram(Rpl_table.d1(cRT_good),D_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table.d1(cRT_else),D_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Distance from Takeoff (m)'); 
        nexttile;
        histogram(Rpl_table.d2(cRT_good),D_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table.d2(cRT_else),D_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Distance from Landing (m)'); 
        nexttile;
        histogram(Rpl_table.t2tko(cRT_good),T_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table.t2tko(cRT_else),T_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Time from Nrst-Takeoff (s)'); 
        nexttile;
        histogram(Rpl_table.t2lnd(cRT_good),T_edges,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');    hold on;
        histogram(Rpl_table.t2lnd(cRT_else),T_edges,'Normalization','probability','facealpha',.1,'edgecolor','none','FaceColor','b');    hold off;
        ylabel('Fraction'); xlabel('Time from Nrst-Landing (s)'); 
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot replays
        figure('units','normalized','outerposition',[.2 .3 .6 .5]);
        tiledlayout(5,20,'TileSpacing','tight');
        for zz=cRT_subset'
            tL = cRt(zz,1);
            tR = cRt(zz,2);
            nexttile;
            for nc= 1:n_place_cells
                plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
            end
            ylim([0 n_place_cells+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[num2str(candidate_rkcorr_c(zz),1)],[num2str((tR-tL)*1000,3), ' ms']});
            title(['p ', num2str(candidate_rkcorr_p(zz),1)]);
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot replays with LFP
        figure('units','normalized','outerposition',[.05 .1 .9 .8]);
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
        
        %=== Plot replays with Position of the bat
        figure('units','normalized','outerposition',[.05 .1 .9 .8]);
        tiledlayout(10,20,'TileSpacing','tight');
        i=1;
        for zz=cRT_subset'
            tL = cRt(zz,1);
            tR = cRt(zz,2);
            sL = Rpl_table.smpL(zz);
            sR = Rpl_table.smpR(zz);
            
            nexttile(20*floor((i-1)/20)+i);
            for nc= 1:n_place_cells
                plot(s{sorted_tko_plcCells(nc)}, nc*ones(size(s{sorted_tko_plcCells(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
            end
            ylim([0 n_place_cells+1]); xlim([tL, tR]); yticks([]);
            xticks([tL, tR]);   xticklabels({'',[num2str((tR-tL)*1000,3), ' ms']});
            nexttile(20*(1+floor((i-1)/20))+i);
            plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
            xlim(r_lim(1,:)); ylim(r_lim(2,:));   zlim(r_lim(3,:));  view(0,90);    hold on;
            for ii=1:n_fclus, plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(j,:));end
            textscatter(avg_takeoff(1),avg_takeoff(2),">>");
            scatter3(r(sL:sR,1),r(sL:sR,2),r(sL:sR,3),50,'MarkerFaceColor', 'r','MarkerEdgeColor','none');
            hold off;   axis equal;  axis off;
            i=i+1;
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
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
        
        %=== DECODING Plotting
        figure('units','normalized','outerposition',[.2 .3 .6 .5]);
        tiledlayout(5,20,'TileSpacing','compact');
        for zz=cRT_subset'
            
            %             %=== Find closest bins for decoding start/stop (OLD VERSION, SLOW)
            %             i_bin_strt = knnsearch(st_times',cRt(zz,1));    % Find closest bin for starting
            %             i_bin_stop = knnsearch(st_times',cRt(zz,2));    % Find closest bin for stopping
            
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
                for nc=1:n_dec_cells,n_vector(nc) = histcounts(s{NP_subtable_dec.cell_id(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
                %=== Decoding
                p_x = double(p_x>-1);   % Uniform Prior
                p_dec_global = [p_dec_global, decode_1Dpos_AF_v0(n_vector,t_bin_dur,p_x,f_x)];
            end
            nexttile;   imagesc(imgaussfilt(p_dec_global,smooth_f),prctile(p_dec_global,[1 99],'all')');
            colormap('hot');  ylabel('Spatial bin'); axis off;  set(gca,'YDir','normal');
        end
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
        %=== Plot when Replays are happening
        figure('units','normalized','outerposition',[0 .4 1 .4]);
        tiledlayout(2,1,'TileSpacing','none');
        ax = [];    ax(1) = nexttile;   hold on;
        area(t,bflying,0,'FaceColor',[0 0 0],'FaceAlpha',0.2,'LineStyle','none');
        plot((t(f_clus.strt_frame(f_clus.id == j)')*[1 1])',(ones(size(f_clus.strt_frame(f_clus.id == j)'))*[0 1])','Color',col_clus(j,:));
        plot((cRt(cRT_subset,1)*[1 1])',(ones(size(cRt(cRT_subset,1)))*[0 -.5])','Color','r');    set(gca,'TickLength',[0 0]);  box('off'); ylabel('Replays/Flights');  yticks([]);
        ax(2) = nexttile; plot(t_Hi,LFP_mn);    ylabel('LFP');
        xlabel('Time (s)'); linkaxes(ax,'x');   xlim('tight');  box('off'); yticks([]);
        sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
        fig_count = saveFig(analysis_directory,[cell2mat(unique_ID)],fig_count,options.savefigures);
        
    end
end
