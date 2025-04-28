function NP_Clean_Duplicates_AF_v1()
%% SCRIPT FOR LOOKING AT NP RECORDINGS WITH MULTIPLE PROBES IN THE HIPPOCAMPUS DURING NAVIGATION
% Updated by A.F. on May 2024, based on previous versions
% BRIEF DESCRIPTION (TBD)

%% LOAD DATA and SAVING OPTIONS
for hide=1
    
    %=== Load data
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);    % Behavioral file
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);    % IMU data
    nrnFile = dir('SU_kilosort*');                                                      % Single Units
    unique_ID = options.unique_ID;                                                      % Session identifier
    NP_unit = table();                  MUA_unit = table();                             % Single Units and MUA structures
    n_probes = size(nrnFile,1);                                                         % Number of probes
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        out.good_units.probe_id = i*ones(size(out.good_units,1),1);
        out.mua_units.probe_id = i*ones(size(out.mua_units,1),1);
        NP_unit = [NP_unit; out.good_units];                                            % Single units and MUA are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;                                                                          % Clear unnecessary variables
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);             % Add average firing rate
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);           % Add average firing rate
  
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
    
    %=== Params
    if isrow(t),t=t';end                                    % Make sure t is a column vector
    NP_unit_FS = NP_unit(NP_unit.fr>2,:);                   % Store high firing cells (putative interneurons) in a different variable (def 1 Hz)
    NP_unit = NP_unit(NP_unit.fr<2,:);                      % Exclude neurons with high-firing rate
    n_cells = size(NP_unit,1);                              % Number of units
    n_cells_FS = size(NP_unit_FS,1);                        % Number of units (putative interneurons)
    bflying = logical(bflying);                             % Convert to logical
    n_rep = 10;                                            % Number of repetitions for shuffling (Place Cells)
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
    NP_unit.id = [1:n_cells]';
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
    
    %=== Extract flight periods from accelerometer
    a_abs_NP = vecnorm(NP_imu.acc,2,2);                     % Absolute acceleration
    a_flt_NP = bandpass(a_abs_NP,[7 9],Fs_Hi);              % Filtered at the wing-beat frequency
    
 
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

%% EVALUATE DUPLICATES
for hide=1
    
    merges = cell(n_probes,1);
    
    for i=1:n_probes
    
        NP_subunit = NP_unit([NP_unit.probe_id]==i,:);
        cell_pairs = nchoosek(1:size(NP_subunit,1),2);
        n_pairs = size(cell_pairs,1);
        
        %=== Calculate similarity between waveforms
        template_similarity = zeros(n_pairs,1);
        for nn=1:n_pairs
            template_similarity(nn) = corr(NP_subunit(cell_pairs(nn,1)).template(:),NP_subunit(cell_pairs(nn,2)).template(:));
        end
    
        %=== Calculate similarity between place fields
        mp_corr = zeros(n_pairs,n_surv_clusters);       % Correlation between place fields
        pk_dist = zeros(n_pairs,n_surv_clusters);       % Distance between peaks
        ph_dist = zeros(n_pairs,n_surv_clusters);       % Distance between phase max
        sf_minb = zeros(n_pairs,n_surv_clusters);       % Min number of spikes per flight
        pr_minb = zeros(n_pairs,n_surv_clusters);       % Min peak rate
        pk_crco = zeros(n_pairs,1);                     % Peak of cross correlation
        
        for j = setdiff(imp_bat_clusters,1)
            
            smp1_clus = f_clus.strt_frame(f_clus.id == j)';                     % Takeoff sample
            smp2_clus = f_clus.stop_frame(f_clus.id == j)';                     % Landing sample
            for nc = 1:n_cells,[~,s_flight2{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus smp2_clus]);end
            
            for nn=1:n_pairs
                mp_corr(nn,j) = corr(NP_subunit(cell_pairs(nn,1)).f_clus(j).map(1,:)',NP_subunit(cell_pairs(nn,2)).f_clus(j).map(1,:)');
                pk_dist(nn,j) = abs(NP_subunit(cell_pairs(nn,1)).f_clus(j).field_loc-NP_subunit(cell_pairs(nn,2)).f_clus(j).field_loc)*bin_size_1D;
                ph_dist(nn,j) = abs(NP_subunit(cell_pairs(nn,1)).f_clus(j).phase_max-NP_subunit(cell_pairs(nn,2)).f_clus(j).phase_max);
                sf_minb(nn,j) = min([NP_subunit(cell_pairs(nn,1)).f_clus(j).sum_spk,NP_subunit(cell_pairs(nn,2)).f_clus(j).sum_spk])/NP_subunit(cell_pairs(nn,1)).f_clus(j).n;
                pr_minb(nn,j) = min([NP_subunit(cell_pairs(nn,1)).f_clus(j).peakHz,NP_subunit(cell_pairs(nn,2)).f_clus(j).peakHz]);
                
                %=== Calculate cross correlation
                if j==2
                    [CC_cells,CC_bins] = cross_correlogram_AF_v0(s_flight2{NP_subunit(cell_pairs(nn,1)).id,1},s_flight2{NP_subunit(cell_pairs(nn,2)).id,1},0.1,0.001);
                    [~,tmp_idx] = max(CC_cells);
                    pk_crco(nn,1)= CC_bins(tmp_idx);
                end
            end
        end
        
        
%         %=== Plot cells sorted by flight phase around each cluster
%         for j = setdiff(imp_bat_clusters,1)
%             
%             %=== Extract some variables for plotting
%             smp1_clus = f_clus.strt_frame(f_clus.id == j)';                     % Takeoff sample
%             smp2_clus = f_clus.stop_frame(f_clus.id == j)';                     % Landing sample
%             cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
%             smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
%             
%             %=== Extract features on the cluster and select place cells
%             NP_table = struct2table(NP_unitOnClus{1, j});
%             NP_table.place_cond = NP_table.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
%             NP_table.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
%             NP_table.stab_m>0.4 &...      % Min Stability (DEF: .4)
%             NP_table.sff<0.7;             % Min Peakyness (DEF: .7)
%             NP_subtable = NP_table(NP_table.place_cond,:);
%             n_place_cells = size(NP_subtable,1);
%             NP_subtable_sorted = sortrows(NP_subtable,'phase_max');
%             sorted_tko_plcCells = NP_subtable_sorted.cell_id;
%         
%             %=== Plot place cells sorted by peak location
%             figure('units','normalized','outerposition',[0 0 .5 1]);
%             tiledlayout(n_place_cells,1,'TileSpacing','none');
%             for ncc=flip(1:n_place_cells)
%                 nexttile;   Raster_TimeWrap_AF_v1(s{sorted_tko_plcCells(ncc),1},t(smp1_clus),t(smp2_clus),[],[],0.1,[],0);    box('Off');
%                 if ncc<n_place_cells
%                     set(gca, 'YTick', [],'XLabel',[],'Box','off','XTick', []);
%                 else
%                     set(gca, 'YTick', [],'XLabel',[],'Box','off');
%                 end
%                 ylabel(num2str(sorted_tko_plcCells(ncc)));
%             end
%             xticks([0, 1]); xticklabels({'Takeoff','Landing'});
%             sgtitle([unique_ID{1},', Bat ',unique_ID{2},', ',unique_ID{3}],'Interpreter','none');
%             
%             %=== Check manually picked duplicates
%             n1 = 27;    
%             n2 = 38;
%             temp_idx = find(all(cell_pairs == sort([n1,n2]),2));
%             display(template_similarity(temp_idx))
%             display(mp_corr(temp_idx,j)); 
%             display(pk_dist(temp_idx,j));
%             display(ph_dist(temp_idx,j));
%             display(pk_crco(temp_idx,1));
%             [CC_cells,CC_bins] = cross_correlogram_AF_v0(s_flight2{n1,1},s_flight2{n2,1},0.1,0.001);
%             nexttile;   area(CC_bins,CC_cells,'FaceColor','k');
%             nexttile;   plot(NP_subunit(cell_pairs(temp_idx,1)).template,'r');  hold on;
%             plot(NP_subunit(cell_pairs(temp_idx,2)).template,'k');
% 
%         end
%         
%         
        %=== Plot with putative duplicates
        figure('units','normalized','outerposition',[.3 0 .3 1]);
        tiledlayout(10,3,'TileSpacing','tight');
        unique_nplet = {};
        for j = setdiff(imp_bat_clusters,1)
            
            %=== Identify potential duplicates: similar waveform and place fields OR identycal place fields (without necessarily similar waveform)
            pair_merges = [];
            sim_id = (template_similarity>0.7 & mp_corr(:,j)>0.5 & pk_dist(:,j)<1.5 & ph_dist(:,j)<15) | (mp_corr(:,j)>0.85 & pk_dist(:,j)<0.5 & ph_dist(:,j)<15 & sf_minb(:,j)>=1 & pr_minb(:,j)>=3 & abs(pk_crco(:))<0.005);
            %sim_id = (template_similarity>0.7 & mp_corr(:,j)>0.5 & pk_dist(:,j)<1.5 & ph_dist(:,j)<15) ;
            
            nexttile;   scatter(template_similarity,mp_corr(:,j),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);    
            hold on;    scatter(template_similarity(sim_id),mp_corr(sim_id,j),15,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   xlim([-1.1 1.1]);
            nexttile;   scatter(template_similarity,pk_dist(:,j),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);     xlim([-1 1]);
            hold on;    scatter(template_similarity(sim_id),pk_dist(sim_id,j),15,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   xlim([-1.1 1.1]);
            nexttile;   scatter(template_similarity,ph_dist(:,j),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3);    xlim([-1 1]);
            hold on;    scatter(template_similarity(sim_id),ph_dist(sim_id,j),15,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   xlim([-1.1 1.1]);
            
            %=== Find duplicates of each cell
            pair_merges = cell_pairs(sim_id,:);
            if ~isempty(pair_merges)
                unique_nplet = [unique_nplet;   find_equiv_classes_AF_v0(mat2cell(pair_merges, ones(size(pair_merges,1),1), 2))];
            end

        end
        
        %=== Reduce to n-plets for merging (save the corresponding cluster id
        temp_ids = find_equiv_classes_AF_v0(unique_nplet);
        merges(i,1) = {cellfun(@(x) [NP_subunit(x).cluster_id],temp_ids,'UniformOutput',false)};
        
    end
    
    %=== Save locally
    save(['merges_',unique_ID{1},'_',unique_ID{2},'_',unique_ID{3},'.mat'],'merges');
    close all
    
end
