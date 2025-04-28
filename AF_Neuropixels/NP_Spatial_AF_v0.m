%% Exploratory Script for looking at NP recordings
% Data are obtained from recordings with NP probes implanted in the hippocampus of a freely flying bat.
% Positional data are pre-processed by AF
% Neural data are pre-processed by KQ

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
% rec_bat: tag of the implanted bat
% t: behavioral time, sampled at 100 Hz and starting at the first Master-9 TTL
% T: Total number of behavioral time samples
% t_NP: Neuropixels time, synchronized at with the behavioral time (Ask KQ regarding sampling)
% v_abs(time samples, recorded tag): vector with the absolute velocity of the tags
% s{number of units, 1}: cell with the time (in s) of the spikes from each single unit (will be created below)

%% LOAD DATA

%=== Concatenate Single Units after Manual Curation
probe_idx = [1,2,3];
NP_unit = table();
for i=1:numel(probe_idx)
    load(['Sorted_units_AF\SU_kilosort_outdir_probe',num2str(probe_idx(i)),'.mat']);
    out.good_units.group = NaN(size(out.good_units.group));                         % This is the structure containing the single units for each probe
    NP_unit = [NP_unit; out.good_units];                                            % Single units are assembled into a table
end

%=== Get the time associated with NP recordings (Time is already synchronized with behavioral time t)
t_NP = out.allGlobalSpikeTimes_usec/1e6;    clear out;                              % This is the time recorded by the Neuropixels system

%=== Behavioral data
BHV1_file = dir(fullfile(cd,'Ext_Beh*','Extracted_Behavior*'));                                        % Preprocessed Behavioral Data
load([BHV1_file.folder,'\',BHV1_file.name],'a','a_abs','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
    'batdate','bflying','Fs','f_num','n_tags','r','r_lim','t','T','v','v_abs','v_th','wBeats');        % Processed Behavioral Data
BHV2_file = dir(fullfile(cd,'Ext_Beh*','Analysis*/Analyzed_DATA*'));
load([BHV2_file.folder,'\',BHV2_file.name]);
load('rec_bat.mat');                                                                                   % id of the recorded bat

%% RECLUSTER FLIGHTS

alpha_clus = 1.4;   clear f_clus
f_clus = FlightClus_AF_v2(squeeze(r(:,:,rec_bat)),bflying(:,rec_bat),'Alpha',alpha_clus,'Frechet',1,'Points',7);


%% BASIC PARAMS, ASSIGNMENTS AND SANITY CHECKS

%=== Params
NP_unit = NP_unit(NP_unit.fr<1,:);                      % Exclude neurons with high-firing rate
%NP_unit = NP_unit(NP_unit.Amplitude>1e3,:);            % Exclude neurons with low amplitude
n_cells = size(NP_unit,1);                              % Number of units
r_lim = [-4.6 5.5; -2.3 1.8; 0 2.8]*1.1;                % Override Room boundaries
bflying = logical(bflying);                             % Convert to logical
n_rep = 100;                                            % Number of repetitions for shuffling (Place Cells)
times_to_shift = t(randi([0.5*Fs T-0.5*Fs],1,n_rep));   % Time intervals for circshift
bin_size_1D = 0.15;                                     % Bin Size for 1D Firing Maps
min_flights_with_spikes = 3;                            % Minimum number of flights with spikes to compute spatial info
options.savefigures = 1;                                % Save Figures
fig_count = 1;                                          % Id of the first figure
rstr_int = [-3 3];                                      % Interval for looking at rasters
t_tko = t(f_smp{rec_bat,1});                            % Times of takeoff
t_lnd = t(f_smp{rec_bat,2});                            % Times of landing
imp_bat_clusters =  unique(f_clus.id);                  % Surviving fligth clusters
n_surv_clusters = numel(imp_bat_clusters);              % Number of surviving flight clusters
col_clus = hsv(n_surv_clusters);                        % Colors for clusters
min_time_2D_fly = 0.2;                                  % Minimum time for 2D maps (flight)
plot_field1D = 0;                                       % If plotting the 1D fields

%=== Sanity Checks
figure('units','normalized','outerposition',[.3 .2 .5 .4]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;   plot(t,ones(size(t)),'.');  hold on;    plot(t_NP,1.1*ones(size(t_NP)),'.');  ylim([0.9 1.2]);
xlabel('Time (s)'); legend('Ciholas','Neuropixels');    yticks([]); title('Sampled Time');  xlim('tight');
nexttile;   plot(t(end)*[1 1],[0 1]);   hold on;    plot(t_NP(end)*[1 1],[1 2]);
xlabel('Time (s)'); title('End Points');

%=== Populate the s cell with the single units and calculate firing rate
s = cell(n_cells,n_rep);            % Single units' spikes
Rate = zeros(length(t),n_cells);    % Smoothed Firing Rate
for nc = 1:n_cells
    s{nc,1} = NP_unit.spikeTimes_usec{nc,1}/1e6;                    % Each row contains the time (in s) of the spikes from a single unit
    Rate(2:end,nc) = histcounts(s{nc,1},t)*Fs;                      % Calculate firing rate at behavioral samples
    Rate(:,nc) = smoothdata(Rate(:,nc),1,'gaussian',round(Fs*0.1)); % Smooth firing rate
    
    for n = 1:n_rep;    s{nc,n+1} = mod(s{nc,1}+times_to_shift(n),t(T)); end        % Shuffling
end
Rate_norm = normalize(Rate,1,'range');

%=== Convert NP unit to structure
NP_unit = table2struct(NP_unit);

%% LOOK AT POPULATION RASTER PLOT (UNSORTED AND SORTED BASED ON MEDIAN TIME RELATIVE TO TAKEOFF)
for hide=1
    
    %=== Plot raster with unsorted units
    figure('units','normalized','outerposition',[0 0 1 1]);
    for nc= 1:n_cells
        plot(s{nc}, nc*ones(size(s{nc})), 'k|','MarkerSize', round(n_cells*0.01));   hold on;          % Raster for each cell
    end
    area(t,bflying(:,rec_bat)*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  xlim('tight');  xlabel('Time (s)'); ylabel('Unit #');
    
    %=== Calculate for each cell median time spike to takeoff
    t2tko = zeros(n_cells,1);
    for nc =1:n_cells
        tmp_idx = knnsearch(s{nc},t_tko);
        t2tko(nc) = median(t_tko-s{nc}(tmp_idx));
    end
    [~,sorted_tko] = sort(t2tko);
    
    %=== Plot raster with sorted units
    figure('units','normalized','outerposition',[0 0 1 1]);
    for nc= 1:n_cells
        plot(s{sorted_tko(nc)}, nc*ones(size(s{sorted_tko(nc)})), 'k|','MarkerSize', round(n_cells*0.01));   hold on;          % Raster for each cell
    end
    area(t,bflying(:,rec_bat)*n_cells,0,'FaceColor',[0 0 1],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 n_cells]);  xlim('tight');  xlabel('Time (s)'); ylabel('Unit #');
    
end

%% LOOK AT THE SPIKES AROUND THE SAME FLIGHT CLUSTER
for hide=1
    
    %=== Params and initializations
    clus_id = 3;
    smp1_clus = f_clus.strt_frame(f_clus.id == clus_id)';
    smp2_clus = f_clus.stop_frame(f_clus.id == clus_id)';
    n_fclus = numel(smp1_clus);
    t_ic = [-2 2];
    d_fclus = mean(t(smp2_clus-smp1_clus));
    
    %=== Get the spike times happening around flight
    s_flight = cell(n_cells,1);
    for nc = 1:n_cells
        [~,s_flight{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);
    end
    
    %=== Calculate for each cell median time spike to takeoff
    t2tko = NaN(n_cells,1);
    for nc =1:n_cells
        tmp_idx = knnsearch(s_flight{nc},t(smp1_clus));
        if ~isempty(tmp_idx)
            t2tko(nc) = median(t(smp1_clus)-s_flight{nc}(tmp_idx));
        end
    end
    [~,sorted_tko] = sort(t2tko);
    
    %cond_tmp = all([sort(t2tko)>t_ic(1),sort(t2tko)<(t_ic(2)+mean(t(smp2_clus-smp1_clus)))],2);
    cond_tmp = true(n_cells,1);
    n_cells_ok = sum(cond_tmp);
    sorted_tko_ok = sorted_tko(cond_tmp);
    
    %=== Visualize cluster
    figure('units','normalized','outerposition',[0.35 .5 0.3 .4]);
    id = find(f_clus.id==clus_id);
    plot3(r(:,1,rec_bat),r(:,2,rec_bat),r(:,3,rec_bat),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
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
        area(t,bflying(:,rec_bat)*n_cells_ok,0,'FaceColor',col_clus(clus_id,:),'FaceAlpha',.1,'LineStyle','none');
        yticks([]); xticks([]); ylim([0 n_cells_ok]);   box off;
        if i==1, ylabel('Unit #'); yticks([20:20:n_cells_ok]); end
        %else, ax(i).YAxis.Visible = 'off';end
    end
    linkaxes(ax,'y');   sgtitle(['Activity during cluster ',num2str(clus_id), ', Time to takeoff: [', num2str(t_ic(1)), ',' num2str(d_fclus+t_ic(2),2), '] s'])
    
end

%% PLACE CELLS ON FLIGHT CLUSTERS
for hide=1
    for nc=1:n_cells
        
        disp(['Place Field Analysis for Cell ',num2str(nc)])
        for j = imp_bat_clusters
            
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
            
            %=== Calculate linearized trajectories and average 3d path
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
            s_temp=histcounts(s_flight,[0:1:T]/Fs);                                         % Number of above spikes in each sample [0,T]
            s_cut = s_temp(batinclus==1)';                                                  % Number of above spikes in each valid flight j sample
            s_pos = [];
            for nr=1:max(s_cut)                                                             % Get the linearized positions associated with those spikes
                s_pos = [s_pos; flight_pos_1D(s_cut>0,:)];
                s_cut = s_cut-1;
            end
            spikes_pos_1D = s_pos;                                                          % Store the data
            
            %=== Initialize and Calculate place properties
            map_1 = NaN(n_rep+1,numel(bin_size_1D:bin_size_1D:mean(fc_lgt)));SI_1 = NaN(1,n_rep+1);SP_1 = NaN(1,n_rep+1);
            [map_1(1,:),SI_1(1),SP_1(1)] = RateMap_AF_v2(spikes_pos_1D,flight_pos_1D,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
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
                    [map_1(n,:),SI_1(n),SP_1(n)] = RateMap_AF_v2(spikes_pos_sh1(:,:,n),flight_pos_1D,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
                end
            end
            
            %=== Split trials into even and odd (or just random half)
            idx_lo = [1:2:n_fclus]';                 % sort(datasample([1:size(F_temp,1)],ceil(size(F_temp,1)/2),'Replace',false))';
            idx_hi = setdiff([1:n_fclus]',idx_lo);
            
            %=== Calculate correlation between paths falling into those 2 categories
            path_corrOdEv = trace(corr(mean(flight_rsmp(:,:,idx_lo),3),mean(flight_rsmp(:,:,idx_hi),3)))/3;
            
            %=== Map Calculation for odd/even trials
            [spikes_pos_1DOd,flight_pos_1DOd] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_lo',t);
            [map_1Od,~,~] = RateMap_AF_v2(spikes_pos_1DOd,flight_pos_1DOd,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            [spikes_pos_1DEv,flight_pos_1DEv] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_hi',t);
            [map_1Ev,~,~] = RateMap_AF_v2(spikes_pos_1DEv,flight_pos_1DEv,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
            %=== Calculate correlation between maps
            map_corrOdEv = corr(map_1Ev',map_1Od','Type','Spearman');
            
            %=== Calculate distance between maps
            map_distOdEv = pdist2(map_1Ev,map_1Od)/size(oned_ctrs,2);
            
            %=== Find the closest positional sample to each spike (plotting purposes)
            if ~isempty(s_flight)
                k = round(s_flight*100);    k = k(k<T);
                spikes_pos = r(k,:,rec_bat);
            else, spikes_pos = [nan nan nan];end
            
            %=== Save 1D field information
            NP_unit(nc).f_clus(j).id = j;
            NP_unit(nc).f_clus(j).n = size(flight_pos,3);
            NP_unit(nc).f_clus(j).SI = SI_1;
            NP_unit(nc).f_clus(j).SP = SP_1;
            NP_unit(nc).f_clus(j).map = map_1;
            NP_unit(nc).f_clus(j).binC = oned_ctrs;
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
            
        end
        
        %=== Store features of the unit within a cluster of interest
        clus_id = 3;
        NP_unit(nc).cell_id = nc;
        NP_unit(nc).SI_best = [NP_unit(nc).f_clus(clus_id).SI_value];
        NP_unit(nc).pv_best = [NP_unit(nc).f_clus(clus_id).SI_p_val];
        NP_unit(nc).sum_spk = [NP_unit(nc).f_clus(clus_id).sum_spk];
        NP_unit(nc).plc_map = {NP_unit(nc).f_clus(clus_id).map(1,:)'};
        NP_unit(nc).plc_ctr = {NP_unit(nc).f_clus(clus_id).binC'};
        
    end
    
    %=== Create table for easier sorting and extract place cells
    NP_table = struct2table(NP_unit);
    NP_table = sortrows(NP_table,'pv_best','ascend');
    place_cond = NP_table.pv_best<0.05 & NP_table.SI_best>1.5; % & NP_table.sum_spk>10;
    cell_place_index = find(place_cond);
end

%% PLOT THE BEST 10 CELLS
for hide=1
    
    %=== Plot
    for ncc=1:10
        
        nc = NP_table.cell_id(ncc);
        figure('units','normalized','outerposition',[.3 0 .6 1]);
        tiledlayout(numel(imp_bat_clusters),9,'TileSpacing','tight');
        
        for j = imp_bat_clusters
            
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
            
            %=== Calculate linearized trajectories and average 3d path
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
            s_temp=histcounts(s_flight,[0:1:T]/Fs);                                         % Number of above spikes in each sample [0,T]
            s_cut = s_temp(batinclus==1)';                                                  % Number of above spikes in each valid flight j sample
            s_pos = [];
            for nr=1:max(s_cut)                                                             % Get the linearized positions associated with those spikes
                s_pos = [s_pos; flight_pos_1D(s_cut>0,:)];
                s_cut = s_cut-1;
            end
            spikes_pos_1D = s_pos;                                                          % Store the data
            
            %=== Initialize and Calculate place properties
            map_1 = NaN(n_rep+1,numel(bin_size_1D:bin_size_1D:mean(fc_lgt)));SI_1 = NaN(1,n_rep+1);SP_1 = NaN(1,n_rep+1);
            [map_1(1,:),SI_1(1),SP_1(1)] = RateMap_AF_v2(spikes_pos_1D,flight_pos_1D,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
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
                    
                    [map_1(n,:),SI_1(n),SP_1(n)] = RateMap_AF_v2(spikes_pos_sh1(:,:,n),flight_pos_1D,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
                end
            end
            
            %=== Split trials into even and odd (or just random half)
            idx_lo = [1:2:n_fclus]';                 % sort(datasample([1:size(F_temp,1)],ceil(size(F_temp,1)/2),'Replace',false))';
            idx_hi = setdiff([1:n_fclus]',idx_lo);
            
            %=== Calculate correlation between paths falling into those 2 categories
            path_corrOdEv = trace(corr(mean(flight_rsmp(:,:,idx_lo),3),mean(flight_rsmp(:,:,idx_hi),3)))/3;
            
            %=== Map Calculation for odd/even trials
            [spikes_pos_1DOd,flight_pos_1DOd] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_lo',t);
            [map_1Od,~,~] = RateMap_AF_v2(spikes_pos_1DOd,flight_pos_1DOd,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            [spikes_pos_1DEv,flight_pos_1DEv] = getSpikePlace1D_AF_v0(s{nc,1},smp1_clus,smp2_clus,fc_lin,idx_hi',t);
            [map_1Ev,~,~] = RateMap_AF_v2(spikes_pos_1DEv,flight_pos_1DEv,0,'1d',bin_size_1D/mean(fc_lgt),min_time_2D_fly);
            
            %=== Calculate correlation between maps
            map_corrOdEv = corr(map_1Ev',map_1Od','Type','Spearman');
            
            %=== Calculate distance between maps
            map_distOdEv = pdist2(map_1Ev,map_1Od)/size(oned_ctrs,2);
            
            %=== Find the closest positional sample to each spike (plotting purposes)
            if ~isempty(s_flight)
                k = round(s_flight*100);    k = k(k<T);
                spikes_pos = r(k,:,rec_bat);
            else, spikes_pos = [nan nan nan];end
            
            %=== Plot
            nexttile([1 2]);   plot3(flight_pos_rsh(1,:),flight_pos_rsh(2,:),0*flight_pos_rsh(3,:),'.','MarkerSize',1,'MarkerEdgeColor',.9*[1 1 1]);
            hold on;    textscatter(flight_pos_rsh(1,1)-0.5,flight_pos_rsh(2,1),">>");
            scatter3(spikes_pos(:,1),spikes_pos(:,2),spikes_pos(:,3),6,'MarkerFaceColor', 'r','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
            hold off;    axis equal;
            view(0,90); xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));   title(['Cluster', num2str(j),', N = ', num2str(n_fclus)]);
            xticks([]); yticks([]);
            nexttile([1 2]);   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],1,'T-Wrp (ALL)',1);
            nexttile([1 2]);   plot(bin_size_1D:bin_size_1D:mean(fc_lgt),mean(map_1(2:end,:),1),'LineWidth',0.3,'Color',[.6 .6 .6]);   xlabel('Distance from tko (m)');    ylabel('Firing Rate (Hz)');
            hold on; plot(oned_ctrs,map_1(1,:),'r','LineWidth',3);   hold off;  xlim([0 mean(fc_lgt)]);
            nexttile;   plot_shuffled(SI_1);  xlabel('Spatial Information');
            nexttile([1 2]);   ave_trajectory = mean(flight_rsmp,3)';
            map_tmp = interp1(oned_ctrs,map_1(1,:)',mean(length_rsmp,1)*mean(fc_lgt),'linear','extrap')';
            map_color = uint8(round(normalize(map_tmp,'range',[2 100])));   cmap = jet(100);
            p = plot3(ave_trajectory(1,:),ave_trajectory(2,:),ave_trajectory(3,:),'r','LineWidth',4);   axis equal;
            grid on;    cdata = [uint8(cmap(map_color,:)*255) uint8(ones(100,1))].';    drawnow();
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cdata);  view(0,90); xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));
            title(['C (Evn VS Odd): ',num2str(map_corrOdEv,2)]);    axis off;
            
        end
    end
    
end

%% PLOT PLACE CELLS (REDUCED VISUALIZATION)
for hide=1
    
    n2show = 5;
    
    %=== Plot place cells around the selected cluster
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(5,4*n2show,'TileSpacing','tight');
    sgtitle('Example place cells');
    
    for ncc=1:(5*n2show)%sum(place_cond)
        
        nc = NP_table.cell_id(cell_place_index(ncc));
        j = clus_id;
        
        smp1_clus = f_clus.strt_frame(f_clus.id == j)';   % Takeoff sample
        smp2_clus = f_clus.stop_frame(f_clus.id == j)';   % Landing sample
        cond = smp1_clus>1*Fs & smp2_clus<(T-1*Fs);                         % Exclude flights close to edges of recording...
        smp1_clus = smp1_clus(cond);    smp2_clus = smp2_clus(cond);        % ...
        n_fclus = numel(smp1_clus);                                         % Number of surviving flights in the cluster
        
        %=== Linearize flights from start to stop, stack them
        fc_lgt = f_clus.length(f_clus.id == j);       % Length in m of the paths
        fc_lin = {f_clus.lin_tr{1,f_clus.id == j}}';  % Flights linearized between 0 and 1
        flight_pos_1D = []; batinclus = zeros(T,1);
        for m = 1:n_fclus
            flight_pos_1D = [flight_pos_1D; fc_lin{m, 1}];          %=== Stack together flights, parametrized between 0 and 1
            batinclus(smp1_clus(m,:):smp2_clus(m,:)) = 1;           %=== Define logical vector when bat is in cluster [0,T]
        end
        
        %=== Get the spikes happening during cluster j
        [num_spikes_per_flight,s_flight] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus smp2_clus+1]);
        
        %=== Calculate linearized trajectories and average 3d path
        flight_pos = f_clus.pos(:,:,f_clus.id==j);
        flight_pos_rsh = reshape(flight_pos,3,[]);
        flight_pos_rsh = flight_pos_rsh(:,all(~isnan(flight_pos_rsh),1));
        
        %=== Find the closest positional sample to each spike (plotting purposes)
        if ~isempty(s_flight)
            k = round(s_flight*100);    k = k(k<T);
            spikes_pos = r(k,:,rec_bat);
        else, spikes_pos = [nan nan nan];end
        
        %=== Plot
        nexttile([1 2]);   plot3(flight_pos_rsh(1,:),flight_pos_rsh(2,:),0*flight_pos_rsh(3,:),'.','MarkerSize',1,'MarkerEdgeColor',.9*[1 1 1]);
        hold on;    textscatter(flight_pos_rsh(1,1)-0.5,flight_pos_rsh(2,1),">>");
        scatter3(spikes_pos(:,1),spikes_pos(:,2),spikes_pos(:,3),6,'MarkerFaceColor', 'r','MarkerEdgeColor','none','MarkerFaceAlpha',.2);
        hold off;    axis equal;    
        view(0,90); xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  
        xticks([]); yticks([]);
        nexttile([1 2]);   Raster_TimeWrap_AF_v1(s{nc,1},t(smp1_clus),t(smp2_clus),[],[],1,[],1);
        
    end
    
    %=== Plot average firing rate (restricted on place cells) and flights
    figure('units','normalized','outerposition',[0 .3 1 .5]);
    plot(t,mean(Rate_norm,2),'k');  hold on;
    area(t,bflying(:,rec_bat),0,'FaceColor',[1 0 0],'FaceAlpha',0.5,'LineStyle','none');    % Plot flights
    ylim([0 max(mean(Rate_norm,2))]);  xlim('tight');  xlabel('Time (s)'); ylabel('Unit #');
    
end

%% LOOK AT THE SPIKES AROUND THE SAME FLIGHT CLUSTER (RESTRICT TO PLACE CELLS)
for hide=1
    
    %=== Params and initializations
    clus_id = 3;
    smp1_clus = f_clus.strt_frame(f_clus.id == clus_id)';
    smp2_clus = f_clus.stop_frame(f_clus.id == clus_id)';
    n_fclus = numel(smp1_clus);
    t_ic = [-2 2];
    d_fclus = mean(t(smp2_clus-smp1_clus));
    
    %=== Get the spike times happening around flight
    s_flight = cell(n_cells,1);
    for nc = 1:n_cells
        [~,s_flight{nc,1}] = count_spikes_AF_v0(s{nc,1},t,[smp1_clus+t_ic(1)*Fs smp2_clus+t_ic(2)*Fs]);
    end
    
    %=== Calculate for each cell median time spike to takeoff
    t2tko = NaN(n_cells,1);
    for nc =1:n_cells
        tmp_idx = knnsearch(s_flight{nc},t(smp1_clus));
        if ~isempty(tmp_idx)
            t2tko(nc) = median(t(smp1_clus)-s_flight{nc}(tmp_idx));
        end
    end
    [~,sorted_tko] = sort(t2tko);
    
    %=== Restrict to place cells only
    sorted_tko_ok = sorted_tko;
    sorted_tko_ok = intersect(sorted_tko_ok,NP_table.cell_id(cell_place_index),'stable');
    n_cells_ok = numel(sorted_tko_ok);
    
    %=== Visualize cluster
    figure('units','normalized','outerposition',[0.35 .5 0.3 .4]);
    id = find(f_clus.id==clus_id);
    plot3(r(:,1,rec_bat),r(:,2,rec_bat),r(:,3,rec_bat),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
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
        area(t,bflying(:,rec_bat)*n_cells_ok,0,'FaceColor',col_clus(clus_id,:),'FaceAlpha',.1,'LineStyle','none');
        yticks([]); xticks([]); ylim([0 n_cells_ok]);   box off;
        if i==1, ylabel('Unit #'); yticks([20:20:n_cells_ok]); end
        %else, ax(i).YAxis.Visible = 'off';end
    end
    linkaxes(ax,'y');   sgtitle(['Activity during cluster ',num2str(clus_id), ', Time to takeoff: [', num2str(t_ic(1)), ',' num2str(d_fclus+t_ic(2),2), '] s']);
    
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