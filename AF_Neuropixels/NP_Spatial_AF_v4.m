%% Exploratory Script for looking at NP recordings
% Data are obtained from recordings with NP probes implanted in the hippocampus of a freely flying bat.

%% BRIEF DESCRIPTION
for hide=1
    %==========================%
    %=== RELEVANT VARIABLES ===%
    %==========================%
    
    % This is the list of the relevant variables that will be loaded/created in the initial sections of this code:
    
    % a_abs(time samples, recorded tag): vector with the absolute acceleration of the tags
    % bflying(time samples, recorded tag): vector with the flight state for each tag (1 flight, 0 rest)
    % f_num: vector with the total number of flights of each tag
    % f_smp{recorded tag, start/stop}: cell with the samples of takeoff and landing
    % Fs: behavioral sampling frequency (Ciholas system, 100 Hz)
    % NP_unit: table with the curated single units
    % r(time samples, xyz, recorded tag): vector with the 3D position of each tag
    % t: behavioral time, sampled at 100 Hz and starting at the first Master-9 TTL
    % T: Total number of behavioral time samples
    % t_NP: Neuropixels time, synchronized at with the behavioral time (Ask KQ regarding sampling)
    % v_abs(time samples, recorded tag): vector with the absolute velocity of the tags
    % s{number of units, 1}: cell with the time (in s) of the spikes from each single unit (will be created below)
end

%% LOAD DATA
for hide=1
    
    disp('Loading data...');
    bhvFile = dir('Extracted_Behavior_*');  load([bhvFile.folder '/' bhvFile.name]);        
    imuFile = dir('IMU_data.mat');          load([imuFile.folder '/' imuFile.name]);
    nrnFile = dir('SU_kilosort*');
    NP_unit = table();
    MUA_unit = table();
    n_probes = size(nrnFile,1);
    
    for i=1:n_probes
        load([nrnFile(i).folder '/' nrnFile(i).name]);
        %out.good_units.group = NaN(size(out.good_units.group));                        % This is the structure containing the single units for each probe
        NP_unit = [NP_unit; out.good_units];                                            % Single units are assembled into a table
        MUA_unit = [MUA_unit; out.mua_units];
    end
    clear out;
    NP_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),NP_unit.spikeTimes_usec);
    MUA_unit.fr = cellfun(@(x) 1/mean(diff(x/1e6)),MUA_unit.spikeTimes_usec);
    
    %=== Load LFP from Probe1
    load('LFP_probe1.mat');                                                             % This file contains LFP data from probe 1
                                 
end

%% RECLUSTER FLIGHTS
for hide=1
    alpha_clus = .8;   clear f_clus    % Default is 0.7
    f_clus = FlightClus_AF_v3(r,bflying,Fs,'Alpha',alpha_clus,'Frechet',1,'Points',7);
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
    bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
    min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
    options.savefigures = 0;                                % Save Figures
    fig_count = 1;                                          % Id of the first figure
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

%% LOOK AT POPULATION RASTER PLOT (UNSORTED AND SORTED BASED ON MEDIAN TIME RELATIVE TO TAKEOFF)
for hide=1
    
    %=== Plot raster with unsorted units
    figure('units','normalized','outerposition',[0 0 1 1]);
    for nc= 1:n_cells
        plot(s{nc}, nc*ones(size(s{nc})), 'k|','MarkerSize', round(n_cells*0.01));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  xlim('tight');  xlabel('Time (s)'); ylabel('Unit #');
    
    %=== Calculate for each cell median time spike to takeoff
    t2tko = zeros(n_cells,1);
    for nc =1:n_cells
        tgt_seq = knnsearch(s{nc},t_tko);
        t2tko(nc) = median(t_tko-s{nc}(tgt_seq));
    end
    [~,sorted_tko] = sort(t2tko);
    
    %=== Calculate for each cell median time spike to takeoff for a given cluster
    t_tko_clus = t(f_clus.strt_frame(f_clus.id == 2)');   
    t2tko = zeros(n_cells,1);
    for nc =1:n_cells
        tgt_seq = knnsearch(s{nc},t_tko_clus);
        t2tko(nc) = median(t_tko_clus-s{nc}(tgt_seq));
    end
    [~,sorted_tko] = sort(t2tko);
    
    %=== Plot raster with sorted units, LFP and Spike Density
    figure('units','normalized','outerposition',[0 .3 1 .6]);
    tiledlayout(11,1,'TileSpacing','none');
    cx(1) = nexttile(1,[4 1]);
    for nc= 1:n_cells
        plot(s{sorted_tko(nc)}, nc*ones(size(s{sorted_tko(nc)})), 'k|','MarkerSize', round(n_cells*0.01));   hold on;          % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  ylabel('Unit #');   xticks([]); set(gca,'TickLength',[0 0]);
    cx(2) = nexttile(5,[1 1]);  plot(NP_imu.t,a_flt_NP);    ylabel('Accelerometer');    set(gca,'TickLength',[0 0]);
    cx(3) = nexttile(6,[2 1]);   plot(NP_imu.t,LFP_mn); hold on;    plot(NP_imu.t,LFP_th);  xticks([]); legend('All','Theta');  ylabel('LFP (norm)');    set(gca,'TickLength',[0 0]);
    cx(4) = nexttile(8,[2 1]);
    LFP_lw = interp1(NP_imu.t,LFP_mn,t,'linear',median(LFP_mn));
    [PS_LFP,freq_SG] = cwt(LFP_lw,Fs);
    %pcolor(linspace(NP_imu.t(1),NP_imu.t(end),size(PS_LFP,2)),freq_SG,smoothdata(abs(PS_LFP),1,'gaussian',10));    
    imagesc([t(1),t(end)],[freq_SG(1),freq_SG(end)],imgaussfilt(abs(PS_LFP),[2 10])); shading interp;  colormap(hot);
    set(gca, 'YScale', 'log','YDir', 'normal','TickLength',[0 0]);   ylim([1 50]);  yticks([1 5 10 20 50]);    ylabel('Freq (Hz)');
    cx(5) = nexttile(10,[2 1]);
    plot(t_rate,all_rate,'k');   ylabel('Spike Density');   set(gca,'TickLength',[0 0]);
    linkaxes(cx,'x');   xlim('tight');  xlabel('Time (s)');
    [x_range_theta,~] = ginput(2);
    
    %=== Calculate ISI within the indicated range and average across cells
    edges_theta = [0:0.010:60];
    ctrs_theta = edges_theta(1:end-1)+mean(diff(edges_theta))/2;
    ISI_counts = zeros(n_cells,numel(edges_theta)-1);
    Cross_c = [];
    for nc=1:n_cells
       %ISI_counts(nc,:) = histcounts(diff(s{nc,1}(s{nc,1}>x_range_theta(1) & s{nc,1}<x_range_theta(2))),edges_theta,'Normalization', 'countdensity');
       [cross_c,ctrs_theta] = cross_correlogram_AF_v0(s{nc,1}(s{nc,1}>x_range_theta(1) & s{nc,1}<x_range_theta(2)),s{nc,1}(s{nc,1}>x_range_theta(1) & s{nc,1}<x_range_theta(2)),1,0.010);
       Cross_c = [Cross_c;cross_c'];
    end
    figure;
    %plot(ctrs_theta,nanmean(ISI_counts,1)); set(gca,'YScale', 'log');
    plot(ctrs_theta,smoothdata(nanmean(Cross_c,1),'movmean',5)); set(gca,'YScale', 'log');  hold on;
    plot(([0.19 0.19].*[-3:1:3]')',(ylim.*[1 1 1 1 1 1 1]')','k--');
     
end

%% PLACE CELLS ON FLIGHT CLUSTERS
for hide=1
    for nc=1:n_cells
        
        disp(['Place Field Analysis for Cell ',num2str(nc)]);
        
        %=== Loop across fligth clusters
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
            [num_spikes_per_flight,s_flight] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus smp2_clus+1]);
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
                    
                    %                     %=== Use random redistribution of spikes for low number of flights
                    %                     if size(F_temp,1)<10
                    %                         spikes_pos_sh1(:,:,n) = datasample(flight_pos_1D,size(spikes_pos_1D,1));
                    %                     end
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
            NP_unit(nc).f_clus(j).corr_map = map_corrOdEv;
            NP_unit(nc).f_clus(j).dist_map = map_distOdEv;
            NP_unit(nc).f_clus(j).corr_pth = path_corrOdEv;
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
            
        end
        
        %=== Store features of the unit within a cluster of interest (limit place cell analysis to that cluster)
        NP_unit(nc).cell_id = nc;
        NP_unit(nc).SI_best = [NP_unit(nc).f_clus(clus_id).SI_value];
        NP_unit(nc).pv_best = [NP_unit(nc).f_clus(clus_id).SI_p_val];
        NP_unit(nc).sum_spk = [NP_unit(nc).f_clus(clus_id).sum_spk];
        NP_unit(nc).plc_map = {NP_unit(nc).f_clus(clus_id).map(1,:)'};
        NP_unit(nc).plc_ctr = {NP_unit(nc).f_clus(clus_id).binC'};
        NP_unit(nc).prob_x = {NP_unit(nc).f_clus(clus_id).prob_x'};
        NP_unit(nc).corr_m = [NP_unit(nc).f_clus(clus_id).corr_map];

    end
    
    %=== Create table for easier sorting and extract place cells
    NP_table = struct2table(NP_unit);
    NP_table = sortrows(NP_table,'pv_best','ascend');         % if sorting is helpful
    %place_cond = NP_table.pv_best<0.05 & NP_table.SI_best>1; % Criteria for place cells (minimal SI can be relaxed)
    %place_cond = NP_table.pv_best<0.05 & NP_table.SI_best>.5 & NP_table.corr_m>0.4;
    %place_cond = NP_table.pv_best<0.05;
    place_cond = NP_table.sum_spk>0 & NP_table.corr_m>0.4;
    cell_place_index = find(place_cond);                       % Indexes of place cells
    n_place_cells = numel(cell_place_index);                   % Number of place cells
end

%% LOOK AT THE SPIKES AROUND THE SAME FLIGHT CLUSTER (RESTRICT TO PLACE CELLS) + SNAKE PLOT
for hide=1
    
    %=== Redefine Place Cells
    NP_table = struct2table(NP_unit);
    NP_table = sortrows(NP_table,'pv_best','ascend');         % if sorting is helpful
    %place_cond = NP_table.pv_best<0.05 & NP_table.SI_best>1; 
    %place_cond = NP_table.pv_best<0.05 & NP_table.SI_best>.5 & NP_table.corr_m>0.5;
    place_cond = NP_table.sum_spk>0 & NP_table.corr_m>0.5;
    cell_place_index = find(place_cond);                       % Indexes of place cells
    n_place_cells = numel(cell_place_index);      
    
    %=== Params and initializations
    smp1_clus = f_clus.strt_frame(f_clus.id == clus_id)';
    smp2_clus = f_clus.stop_frame(f_clus.id == clus_id)';
    n_fclus = numel(smp1_clus);
    t_ic = [-2 2];
    d_fclus = mean(t(smp2_clus)-t(smp1_clus));
    
    %=== Get the spike times happening around flight
    s_flight = cell(n_cells,1);
    for nc = 1:n_cells
        [~,s_flight{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);
    end
    
    %=== Calculate for each cell median time spike to takeoff
    t2tko = NaN(n_cells,1);
    for nc =1:n_cells
        tgt_seq = knnsearch(s_flight{nc},t(smp1_clus));
        if ~isempty(tgt_seq)
            t2tko(nc) = median(t(smp1_clus)-s_flight{nc}(tgt_seq));
        end
    end
    [~,sorted_tko] = sort(t2tko);
    
    %=== Restrict to place cells only
    sorted_tko_ok = sorted_tko;
    sorted_tko_ok = intersect(sorted_tko_ok,NP_table.cell_id(cell_place_index),'stable');
    n_cells_ok = numel(sorted_tko_ok);
    
    %=== Define rate by pooling spikes from place cells
    all_s_plc = sort(vertcat(s{sorted_tko_ok,1}));
    all_rate_plc = kernel_rate_AF_v1(all_s_plc,0.1,t_rate)/n_cells_ok;
    
    %=== Visualize cluster
    figure('units','normalized','outerposition',[0.35 .5 0.3 .4]);
    id = find(f_clus.id==clus_id);
    plot3(r(:,1),r(:,2),r(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
    xlim([r_lim(1,1) r_lim(1,2)]); ylim([r_lim(2,1) r_lim(2,2)]);   zlim([r_lim(3,1) r_lim(3,2)]);  view(0,90);
    xlabel('x');    ylabel('y');    hold on;
    avg_take_off = [];
    for ii=1:n_fclus
        title(['Cluster' num2str(clus_id) ' (' num2str(n_fclus) ' flights),'])
        plot3(f_clus.pos(1,:,id(ii)),f_clus.pos(2,:,id(ii)),f_clus.pos(3,:,id(ii)),'-','LineWidth',1,'Color', col_clus(clus_id,:));
        avg_take_off = [avg_take_off f_clus.pos(:,1,id(ii))];
    end
    take_off = mean(avg_take_off,2);    textscatter(take_off(1),take_off(2),"Take-off");     hold off;   axis equal;
    
    %=== Plot sorted activity around each flight
    figure('units','normalized','outerposition',[0 .2 1 .35]);
    tiledlayout(1,n_fclus,'TileSpacing','none');
    for i=1:n_fclus
        ax(i) = nexttile;
        for nc = 1:n_cells_ok
            plot(s_flight{sorted_tko_ok(nc)}, nc*ones(size(s_flight{sorted_tko_ok(nc)})), 'k|','MarkerSize', max(round(n_cells_ok*0.01),1));    hold on;
            xlim([t(smp1_clus(i))+t_ic(1) t(smp1_clus(i))+d_fclus+t_ic(2)]);
        end
        area(t,bflying*n_cells_ok,0,'FaceColor',col_clus(clus_id,:),'FaceAlpha',.1,'LineStyle','none');
        yticks([]); xticks([]); ylim([0 n_cells_ok]);   box off;
        if i==1, ylabel('Unit #'); yticks([20:20:n_cells_ok]); end
        %else, ax(i).YAxis.Visible = 'off';end
    end
    linkaxes(ax,'y');   sgtitle(['Activity during cluster ',num2str(clus_id), ', Time to takeoff: [', num2str(t_ic(1)), ',' num2str(d_fclus+t_ic(2),2), '] s']);
    
    %=== Detect candidate events from spike density from place cells only
    [candidate_event_vector,~,candidate_event_str,candidate_event_stp,candidate_event_width] = BiLevel_Segm_AF_v0(zscore(all_rate_plc)',3,1e3,0.1,0.05);
    candidate_event_str = candidate_event_str(candidate_event_width<0.5);
    candidate_event_stp = candidate_event_stp(candidate_event_width<0.5);
    candidate_event_width = candidate_event_width(candidate_event_width<0.5);
    
    %=== Plot sorted activity (entire session)
    figure('units','normalized','outerposition',[0 .2 1 .5]);
    tiledlayout(4,1,'TileSpacing','none');
    bx(1) = nexttile([2 1]);
    for nc= 1:n_cells_ok
        plot(s{sorted_tko_ok(nc)}, nc*ones(size(s{sorted_tko_ok(nc)})), 'k|','MarkerSize', 5);  hold on;             % Raster for each cell
    end
    area(t,bflying*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    plot((t(f_clus.strt_frame(f_clus.id == clus_id)')*[1 1])',(ones(size(f_clus.strt_frame(f_clus.id == clus_id)'))*[0 n_cells])','r');
    ylim([0 n_cells_ok]);  ylabel('Unit #');    xticks([]);
    bx(2) = nexttile([1 1]);  plot(t_rate,all_rate_plc,'k');
    bx(3) = nexttile([1 1]);  area(t_rate,candidate_event_vector);
    ylabel('Ave Population firing rate');   xlim('tight');  xlabel('Time (s)'); linkaxes(bx,'x');
    
    %=== Make snake plot
    Field1D = table();  c = 1;  warning('off','all');
    for nc = 1:numel(cell_place_index)
        %=== Extract map and bin centers
        map_tmp = NP_table.plc_map{cell_place_index(nc),1};
        ctr_tmp = NP_table.plc_ctr{cell_place_index(nc),1};
        
        %=== Interpolate between takeoff and landing, 100 pts
        Field1D.map_interp(c,:) = {interp1(ctr_tmp,map_tmp,linspace(ctr_tmp(1),ctr_tmp(end),100)','linear')'};
        
        %=== Get the location of the maximum, in terms of flight phase
        [~,Field1D.phase_max(c,:)] = max(Field1D.map_interp{c,:});
        
        c = c+1;
    end
    warning('on','all');
    
    %=== Extract snake plot for significant cells and plot
    Snake = cell2mat(Field1D.map_interp);
    Snake = Snake./max(Snake,[],2);
    [~,loc] = cellfun(@max,Field1D.map_interp);
    [~,order] = sort(loc);
    figure('units','normalized','outerposition',[.5 .3 .1 .4]);
    imagesc(Snake(order,:));    colormap(viridis);  ylabel('1D Field');   xlabel('Flight Phase (%)');   colorbar;
    
end

%% PLOT DETECTED REPLAYS
for hide=1
    
    cRt = [t_rate(candidate_event_str)',t_rate(candidate_event_stp)'];
    n_candidate_events = size(cRt,1);
    candidate_perc_act = zeros(n_candidate_events,1);
    candidate_rkcorr_c = zeros(n_candidate_events,1);
    candidate_rkcorr_p = NaN(n_candidate_events,1);
    seq_min_spikes = 3;         % Minimum number of spikes in the sequence        
    seq_min_fr_active = 0.1;    % Minimum fraction of active cells in the sequence
    
    %=== Assess event quality
    for i=1:n_candidate_events
        tL = cRt(i,1);
        tR = cRt(i,2);
        candidate_s = cellfun(@(x) x(x>tL & x<tR), s(sorted_tko_ok,1),'UniformOutput', false);   % Extract the spikes happening during the event
        candidate_perc_act(i) =1-sum(cellfun(@isempty,candidate_s))/n_cells_ok;
        if  numel(vertcat(candidate_s{:}))>seq_min_spikes && candidate_perc_act(i)>seq_min_fr_active
            [candidate_rkcorr_c(i),candidate_rkcorr_p(i)] = spikeSeq_analysis_AF_v0(candidate_s,100);
        end
    end
    
    %=== Get best events
    Rpl_table = table();
    Rpl_table.tL = cRt(:,1);
    Rpl_table.tR = cRt(:,2);
    Rpl_table.f = candidate_perc_act;
    Rpl_table.C = candidate_rkcorr_c;
    Rpl_table.p = candidate_rkcorr_p;

    cRT_subset = find(Rpl_table.f>0.1 & Rpl_table.p<0.05);
    n_event2show = min(100,size(cRT_subset,1));
    [p_rw,p_cl] = optimal_arrangement_AF_v0(n_event2show);
    cRT_subset = datasample(cRT_subset,n_event2show,'Replace',false);
    
    %=== Plot replays
    figure('units','normalized','outerposition',[.1 0 .6 1]);
    %tiledlayout(p_rw,p_cl,'TileSpacing','compact');
    for zz=cRT_subset'
        tL = cRt(zz,1);
        tR = cRt(zz,2);
        nexttile;
        for nc= 1:n_cells_ok
            plot(s{sorted_tko_ok(nc)}, nc*ones(size(s{sorted_tko_ok(nc)})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
        end
        ylim([0 n_cells_ok+1]); xticks(tR); xlim([tL, tR]); yticks([]); xticklabels([num2str((tR-tL)*1000), ' ms']);
        title(['C = ', num2str(candidate_rkcorr_c(zz),2)]);
    end
    
%     %=== Find significant events
%     disp(['Fraction of significant events: ' num2str(sum(candidate_rkcorr_p<0.05)/sum(~isnan(candidate_rkcorr_p)),2)]);
%     fw_replay_idx = find(candidate_rkcorr_c>0 & candidate_rkcorr_p<0.05);   [~,tmp_fw] = sort(candidate_perc_act(fw_replay_idx),'descend'); fw_replay_idx = fw_replay_idx(tmp_fw);
%     rw_replay_idx = find(candidate_rkcorr_c<0 & candidate_rkcorr_p<0.05);   [~,tmp_rv] = sort(candidate_perc_act(rw_replay_idx),'descend'); rw_replay_idx = rw_replay_idx(tmp_rv);
    
end

%% RUN DECODING ON CANDIDATE EVENTS 
for hide=1
    
    %=== Redefine Place Cells
    NP_table = struct2table(NP_unit);                           % Reformatting
    place_cond = NP_table.sum_spk>-1;                           % Select cells for decoding (NP_table.corr_m>0.4);
    cell_place_index = find(place_cond);                        % Indexes of place cells
    n_place_cells = numel(cell_place_index);                    % Number of place cells
    
    %=== Params
    t_bin_dur = 0.010;                                          % Time bin duration
    t_bin_ovl = 0.005;                                          % Time bin overlap
    smooth_f = [.01 .01];                                       % Smoothing sigmas [space bin,time bin]
    
    %=== Initialization
    n_bins = floor((t(T)-t_bin_dur)/(t_bin_dur-t_bin_ovl))+1;   % Number of time bins
    st_times = (0:n_bins-1)*(t_bin_dur - t_bin_ovl);            % Start times
    ed_times = st_times + t_bin_dur;                            % End times
    t_d = zeros(1, 2 * n_bins);                                 % All times (initialize)
    t_d(1:2:end) = st_times;    t_d(2:2:end) = ed_times;        % All times (define)
    t_c = st_times+0.5*(t_bin_dur - t_bin_ovl);                 % Bin centers
    n_vector = zeros(1,n_place_cells);                          % Initialize vector containing the number of spikes fired by each place cell
    p_x = NP_unit(1).prob_x{:};                                 % Probability of being on a given bin of the flight cluster
    c_x = NP_unit(1).plc_ctr{:};                                % Centers of the positional bins
    f_x = zeros(numel(p_x),n_place_cells);                      % Vector with the firing rates of place cells in a given bin of the flight cluster
    for nc=1:n_place_cells, f_x(:,nc) = NP_unit(cell_place_index(nc)).plc_map{:};end
    
    %=== Run decoding on selected epochs
    figure('units','normalized','outerposition',[.1 0 .6 1]);
    tiledlayout(p_rw,p_cl,'TileSpacing','compact');
    for zz=cRT_subset'
        i_bin_strt = knnsearch(st_times',cRt(zz,1));    % Find closest bin for starting
        i_bin_stop = knnsearch(st_times',cRt(zz,2));    % Find closest bin for stopping
        p_dec_global = [];                                          % Initialize
        
        %=== Loop across bins
        for i_dec = i_bin_strt:i_bin_stop
            
            %=== Count spikes of each unit on the time bin
            for nc=1:n_place_cells
                n_vector(nc) = histcounts(s{cell_place_index(nc),1},[t_d(2*i_dec-1),t_d(2*i_dec)]);
            end
            %=== Decoding
            p_x = double(p_x>-1);   % Uniform Prior
            p_dec_global = [p_dec_global, decode_1Dpos_AF_v0(n_vector,t_bin_dur,p_x,f_x)];
            
        end
        nexttile;   imagesc(imgaussfilt(p_dec_global,smooth_f),prctile(p_dec_global,[1 99],'all')');  colormap('hot');  ylabel('Spatial bin'); axis off;   
    end
    
end