%% POOL REPLAYS AND QUANTIFY INTERESTING FEATURES

%=== Load Data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Params for cross-correlation Replay SWR
    min_n = 5;                      % Minimum number of active cells
    min_f = 0.2;                    % Minimum fraction of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.2;                    % Minimum Correlation
    
    %=== Sessions to load (comment sessions you don't want to include)
    sessions2include = {...
        %'Dataset_1','32622','231006';...
        'Dataset_1','32622','231007';...
        'Dataset_1','32622','231008';...
        'Dataset_1','32622','231009';...
        'Dataset_1','32622','231010';...
        %'Dataset_2','14445','231208';...
        %'Dataset_2','14445','231209';...
        %'Dataset_2','14445','231210';...
        %'Dataset_2','14445','231211';...
        %'Dataset_2','14445','231212';...
        %'Dataset_2','14611','231208';...
        %'Dataset_2','14611','231209';...
        %'Dataset_2','14611','231210';...
        %'Dataset_2','14611','231211';...
        %'Dataset_2','14611','231212';...
        'Dataset_3','00000','240402';...
        'Dataset_3','00000','240403';...
        'Dataset_3','00000','240404';...
        'Dataset_3','00000','240405';...
        'Dataset_3','14543','240419';...
        'Dataset_3','14543','240420';...
        'Dataset_3','14543','240421';...
        'Dataset_3','14543','240422';...
        'Dataset_3','29959','240402';...
        'Dataset_3','29959','240403';...
        'Dataset_3','29959','240404';...
        'Dataset_3','29959','240405';...
        ...
        };
    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   SB = []; ids = [];  NP_units = []; MS = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(Rpl_table) && any(cellfun(@(row) isequal(row, Rpl_table.unique_ID(1,:)), num2cell(sessions2include, 2)))
            
            %=== Accumulate the replay table and single units
            SB = [SB; Rpl_table];
            NP_units = [NP_units;   NP_unit];
            MS = [MS; MSC];
            
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except SB RPSW_CrosCorr SWSW_AutoCorr RS_bin_centers NP_units MS
    
end

%% ADD SOME FEATURES TO EACH UNIT
for hide=1
    
    NP_units = struct2table(NP_units);

    N_cells = size(NP_units,1); % Number of cells
    bin_size_1D = 0.15;
    f_l = [];  f_s = []; s_c = []; p_H = []; p_M = [];  s_I = [];
    
    for nc=1:N_cells
        
        %=== Rearrange Place Cell Data
        n_clus = size(NP_units.f_clus{nc,1},2)-1;
        NP_unitOnClus = cell(1,n_clus);
        for j = 1:n_clus
            NP_unitOnClus{1, j} = struct(); % Initialize as struct array
            
            %=== Store features of the unit within a cluster of interest
            NP_unitOnClus{1, j}.f_lenght = NP_units.f_clus{nc,1}(j+1).f_lenght;
            NP_unitOnClus{1, j}.f_duration = NP_units.f_clus{nc,1}(j+1).f_duration;
            NP_unitOnClus{1, j}.fr = NP_units.fr;
            NP_unitOnClus{1, j}.SI = NP_units.f_clus{nc,1}(j+1).SI_value;
            NP_unitOnClus{1, j}.p_val = NP_units.f_clus{nc,1}(j+1).SI_p_val;
            NP_unitOnClus{1, j}.spkPerflight = NP_units.f_clus{nc,1}(j+1).sum_spk/NP_units.f_clus{nc,1}(j+1).n;
            NP_unitOnClus{1, j}.plc_map = NP_units.f_clus{nc,1}(j+1).map(1,:)';
            NP_unitOnClus{1, j}.plc_ctr = NP_units.f_clus{nc,1}(j+1).binC';
            NP_unitOnClus{1, j}.prob_x = NP_units.f_clus{nc,1}(j+1).prob_x';
            NP_unitOnClus{1, j}.corr_m = NP_units.f_clus{nc,1}(j+1).corr_map_OdEv;
            NP_unitOnClus{1, j}.dist_m = NP_units.f_clus{nc,1}(j+1).dist_map_OdEv;
            NP_unitOnClus{1, j}.stab_m = NP_units.f_clus{nc,1}(j+1).corr_map_1h2h;
            NP_unitOnClus{1, j}.stbd_m = NP_units.f_clus{nc,1}(j+1).dist_map_1h2h;
            NP_unitOnClus{1, j}.peakHz = NP_units.f_clus{nc,1}(j+1).peakHz;
            NP_unitOnClus{1, j}.field_loc = NP_units.f_clus{nc,1}(j+1).field_loc;
            NP_unitOnClus{1, j}.field_loc_m = NP_units.f_clus{nc,1}(j+1).field_loc*bin_size_1D;
            NP_unitOnClus{1, j}.n_fields = NP_units.f_clus{nc,1}(j+1).n_fields;
            NP_unitOnClus{1, j}.f_width = NP_units.f_clus{nc,1}(j+1).f_width;
            NP_unitOnClus{1, j}.phase_max = NP_units.f_clus{nc,1}(j+1).phase_max;
            NP_unitOnClus{1, j}.sff = NP_units.f_clus{nc,1}(j+1).sff;
            NP_unitOnClus{1, j}.map_interp = NP_units.f_clus{nc,1}(j+1).map_interp;
            NP_unitOnClus{1, j}.asymm_frct = sum(NP_units.f_clus{nc,1}(j+1).map(1,1:round(NP_units.f_clus{nc,1}(j+1).field_loc)))./sum(NP_units.f_clus{nc,1}(j+1).map(1,:));
            
            NP_unitOnClus{1, j}.place_cond =    NP_unitOnClus{1, j}.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
                NP_unitOnClus{1, j}.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
                NP_unitOnClus{1, j}.stab_m>0.4 &...      % Min Stability (DEF: .4)
                NP_unitOnClus{1, j}.sff<0.7;             % Min Peakyness (DEF: .7)
            
        end
        
        %=== Classify as place cell or not
        place_cond = zeros(n_clus,1);
        spkPerflight = zeros(n_clus,1);
        field_size = zeros(n_clus,1);
        stability = zeros(n_clus,1);
        peakHz = zeros(n_clus,1);
        phaseMax = zeros(n_clus,1);
        spatialInfo = zeros(n_clus,1);
        flight_l = zeros(n_clus,1);
        for j = 1:n_clus
            place_cond(j) = NP_unitOnClus{1, j}.place_cond;
            spkPerflight(j) = NP_unitOnClus{1, j}.spkPerflight;
            field_size(j) = NP_unitOnClus{1, j}.f_width;
            flight_l(j) = NP_unitOnClus{1, j}.f_lenght;
            stability(j) = NP_unitOnClus{1, j}.stab_m;
            peakHz(j)  = NP_unitOnClus{1, j}.peakHz;
            phaseMax(j)  = NP_unitOnClus{1, j}.phase_max;
            spatialInfo(j)  = NP_unitOnClus{1, j}.SI;
        end
        NP_units.place_cond(nc) = any(place_cond);
        NP_units.analz_cond(nc) = any(spkPerflight>1);
        if any(place_cond)
            NP_units.field_size(nc) = mean(field_size(place_cond>0));
            f_l = [f_l; flight_l(place_cond>0)];
            f_s = [f_s; field_size(place_cond>0)];
            s_c = [s_c; stability(place_cond>0)];
            p_H = [p_H; peakHz(place_cond>0)];
            p_M = [p_M; phaseMax(place_cond>0)];
            s_I = [s_I; spatialInfo(place_cond>0)];
        else
            NP_units.field_size(nc) = NaN;
        end
       
    end
    
    disp(['N place cells: ',num2str(sum(NP_units.place_cond)),' out of ',num2str(sum(NP_units.analz_cond)),',',num2str(sum(NP_units.place_cond)/sum(NP_units.analz_cond)*100,2),' %'])
    
    NP_place_cells = NP_units(NP_units.place_cond,:);
    
    figure('units','normalized','outerposition',[.3 .3 .35 .25]);
    tiledlayout(1,5,'TileSpacing','tight');
    nexttile;   histogram(f_s,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');                     xlabel('m');    ylabel('Fraction');  title('Field Size');
    nexttile;   histogram(s_c,[0:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');           xlabel('Correlation');    ylabel('Fraction');  title('Stability');  
    nexttile;   histogram(p_H,[1:2:100],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');   xlabel('Hz');    ylabel('Fraction');  title('Peak Firing');  set(gca, 'XScale', 'log');
    nexttile;   histogram(p_M,[0:10:100],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');          xlabel('Flight Phase');    ylabel('Fraction');  title('Preferred Phase');
    nexttile;   histogram(s_I,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');          xlabel('Bits/Spike');    ylabel('Fraction');  title('Spatial Info');

    figure('units','normalized','outerposition',[.3 .3 .3 .3]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;   scatter(f_l,f_s,10,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Flight Length (m)');    ylabel('Number of participating cells');
    nexttile;   plot_distr_AF_v0(f_s(f_l>4&f_l<8), f_s((f_l>8&f_l<14)), {'Short Flights', 'Long Flights'}, 'SEM', 'Place-Field Size(m)');
    
end

%% LOOK AT ACCELEROMETER AND ECHOLOCATION
for hide=1
   
    abs_acc = vertcat(MS.abs_acc);
    echo_rate = vertcat(MS.echo_rate);
    
    %=== Compare accelerometer and echolocation during flight vs rest
    figure('units','normalized','outerposition',[.3 .3 .2 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   plot_distr_AF_v0(abs_acc(:,1), abs_acc(:,2), {'Flight', 'Rest'}, 'SEM', 'Mean Absolute Deviation from g');    ylim([0 1]);
    nexttile;   plot_distr_AF_v0(echo_rate(:,1), echo_rate(:,2), {'Flight', 'Rest'}, 'SEM', 'Echolocation Rate (Hz)');        ylim([0 15]);
    nexttile;   boxchart(echo_rate);    ylim([0 16]);
    
    
    mean(echo_rate(:,1))
    std(echo_rate(:,1))/sqrt(numel(echo_rate(:,1)))
    mean(echo_rate(:,2))
    std(echo_rate(:,2))/sqrt(numel(echo_rate(:,2)))
    
  
    signrank(echo_rate(:,1),echo_rate(:,2))
    
    %=== Plot time course around takeoff and landing
    figure('units','normalized','outerposition',[.6 .3 .15 .6]);
    tiledlayout(2,2,'TileSpacing','tight');
    nexttile;   data_tmp = horzcat(MS.aa1)';
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');   ylim([0 1]);
    xlabel('Time from takeoff (s)');    ylabel('Absolute deviation from g');
    nexttile;   data_tmp = horzcat(MS.aa2)';
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');    ylim([0 1]);
    xlabel('Time from landing (s)');    yticks([]);    
    nexttile;   data_tmp = horzcat(MS.cr1)';
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');   ylim([0 15]);
    xlabel('Time from takeoff (s)');    ylabel('Eholocation Rate');
    nexttile;   data_tmp = horzcat(MS.cr2)';
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');   ylim([0 15]);
    xlabel('Time from landing (s)');    yticks([]);        
    
end

%% LOOK AT FORWARD VS REVERSE REPLAYS
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 5;                      % Minimum number of active cells
    min_f = 0.3;                    % Minimum fraction of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.2;                    % Minimum Correlation
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');
    [~,~,SB.sessionID] =  unique(string(SB.unique_ID),'rows');% Id for each cluster from each session
    SB.GoodFReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z> min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.GoodRReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z<-min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.Not__Replay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z>-inf  ) && (t<=inf ), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;
    SB.compression = SB.flight_dur./SB.dt;
    SB.cvrg =  cellfun(@(x) x(end)-x(1),SB.place_phases); 
    SB.local1 =  SB.d1<.3;
    SB.local2 =  SB.d2<.3;
    
    session_info_summary = groupsummary(SB,'sessionID','max',{'clus_id'});
    mean(session_info_summary.max_clus_id-1)
    
    %=== Select good events
    FRpl_table = SB(SB.GoodFReplay,:);   nFR = size(FRpl_table,1);  % Forward Replays
    RRpl_table = SB(SB.GoodRReplay,:);   nRR = size(RRpl_table,1);  % Reverse Replays
    NRpl_table = SB(SB.Not__Replay,:);   nNR = size(NRpl_table,1);  % Reverse Replays
    ARpl_table = SB(SB.GoodFReplay | SB.GoodRReplay,:);    nAR = size(ARpl_table,1);  % All Replays
    
    session_info_summary_all = groupsummary(ARpl_table,'sessionID','median',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','d1','d2','local1','local2'});
    
    median(ARpl_table.dt)
    std(ARpl_table.dt)
    median(ARpl_table.f)
    std(ARpl_table.f)

    %=== Display fractions of events
    disp(['Fraction of Forward Replays: ', num2str(nFR/(nFR+nRR),3)]);
    disp(['Total number of good Replays: ', num2str(nFR+nRR,4)]);
    disp(['Median fraction of active cells: ', num2str(median([FRpl_table.f;RRpl_table.f]))]);
    disp(['Median duration (s): ', num2str(median([FRpl_table.dt;RRpl_table.dt])*1e3,3),' p-m',num2str(std([FRpl_table.dt;RRpl_table.dt])*1e3./sqrt(nFR+nRR),3)]);
    disp(['Median correlation (s): ', num2str(median([FRpl_table.C;RRpl_table.C]),3),' p-m',num2str(std([FRpl_table.C;RRpl_table.C])./sqrt(nFR+nRR),3)]);
    
    %=================================================================================================================================================================================================================================
    %=== Plot time course of echolocatio rate and accelerometer around replays !!! REVIEWER !!!
    figure('units','normalized','outerposition',[.3 .3 .2 .6]);
    tiledlayout(2,1,'TileSpacing','tight');
    nexttile;   
    data_tmp = vertcat(ARpl_table.abs_acc{:,1});  
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');   hold on;
    data_tmp = vertcat(FRpl_table.abs_acc{:,2});  
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k'); 
    xlim('tight');
    plot([0 0], ylim,'k--');    plot(xlim,mean(abs_acc(:,1))*[1 1],'k--');  plot(xlim,mean(abs_acc(:,2))*[1 1],'k--');
    xlabel('Time from Replay (s)'); ylabel('Absolute deviation from g');
    nexttile;   
    data_tmp = vertcat(ARpl_table.ech_rat{:,1});  
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');   hold on;
    data_tmp = vertcat(FRpl_table.ech_rat{:,2});  
    plotWinterval_AF_v0(MS(1).tme',mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'k');
    xlim('tight');
    plot([0 0], ylim,'k--');    plot(xlim,mean(echo_rate(:,1))*[1 1],'k--');  plot(xlim,mean(echo_rate(:,2))*[1 1],'k--');
    xlabel('Time from Replay (s)'); ylabel('Echolocation Rate (Hz)');   
    
    abs_acc_R = vertcat(ARpl_table.abs_acc{:,1});
    abs_acc_C = vertcat(ARpl_table.abs_acc{:,2});
    ech_rat_R = vertcat(ARpl_table.ech_rat{:,1});
    ech_rat_C = vertcat(ARpl_table.ech_rat{:,2});
    
    disp(['p at time of replay (accelerometer): ', num2str(signrank(abs_acc_R(:,knnsearch(MS(1).tme,0)),abs_acc_C(:,knnsearch(MS(1).tme,0))))])
    disp(['p at time of replay (echolocation): ', num2str(signrank(ech_rat_R(:,knnsearch(MS(1).tme,0)),ech_rat_C(:,knnsearch(MS(1).tme,0))))])
    
    %=================================================================================================================================================================================================================================
    
    %=== Plot Rank Correlations
    figure('units','normalized','outerposition',[.3 .3 .3 .27]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    histogram(FRpl_table.C,[-1:0.1:1],'Normalization','count','facealpha',.5,'edgecolor','none','FaceColor','b');  hold on;
    histogram(RRpl_table.C,[-1:0.1:1],'Normalization','count','facealpha',.5,'edgecolor','none','FaceColor','r');
    histogram(NRpl_table.C,[-1:0.1:1],'Normalization','count','facealpha',.5,'edgecolor','none','FaceColor','k');
    xlabel('Rank-Correlation');   ylabel('Counts');
    nexttile;
    histogram([FRpl_table.C;RRpl_table.C],[-1:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  hold on;
    histogram(NRpl_table.C,[-1:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');
    xlabel('Rank-Correlation');   ylabel('Fraction');
    
    %=== Plot some stats
    figure('units','normalized','outerposition',[.2 .3 .6 .25]);
    tiledlayout(1,8,'TileSpacing','tight');
    nexttile;   histogram(FRpl_table.C,[0:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Rank Correlation');   ylabel('Fraction');
    nexttile;   histogram(FRpl_table.f,[0:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Fraction of participating cells');   ylabel('Fraction');
    nexttile;   histogram(FRpl_table.n,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Number of participating cells');   ylabel('Fraction');
    nexttile;   histogram(FRpl_table.dt,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Duration (s)');   ylabel('Fraction');
    nexttile;   histogram(FRpl_table.derived_dt,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Derived Duration (s)');   ylabel('Fraction');
    nexttile;   histogram(FRpl_table.compression,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','b');  
    xlabel('Compression Factor');   ylabel('Fraction');
    nexttile;   scatter(FRpl_table.flight_len,FRpl_table.n,10,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Flight Length (m)');    ylabel('Number of participating cells');
    nexttile;   scatter(FRpl_table.flight_len,FRpl_table.C,10,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Flight Length (m)');    ylabel('Rank Correlation');
    
%     scatter(FRpl_table.dt,FRpl_table.derived_dt,10,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);    xlim([0 1]);    ylim([0 1]);
    
    %=== Plot some stats
    figure('units','normalized','outerposition',[.2 .1 .6 .25]);
    tiledlayout(1,8,'TileSpacing','tight');
    nexttile;   histogram(RRpl_table.C,[-1:0.1:0],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Rank Correlation');   ylabel('Fraction');
    nexttile;   histogram(RRpl_table.f,[0:0.1:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Fraction of participating cells');   ylabel('Fraction');
    nexttile;   histogram(RRpl_table.n,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Number of participating cells');   ylabel('Fraction');
    nexttile;   histogram(RRpl_table.dt,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Duration (s)');   ylabel('Fraction');
    nexttile;   histogram(RRpl_table.derived_dt,[0:0.05:1],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Derived Duration (s)');   ylabel('Fraction');
    nexttile;   histogram(RRpl_table.compression,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');  
    xlabel('Compression Factor');   ylabel('Fraction');
    nexttile;   scatter(RRpl_table.flight_len,RRpl_table.n,10,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Flight Length (m)');    ylabel('Number of participating cells');
    nexttile;   scatter(RRpl_table.flight_len,RRpl_table.C,10,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Flight Length (m)');    ylabel('Rank Correlation');
    
    %=== Look at temporal and spatial distribution of replays
    figure('units','normalized','outerposition',[.3 .2 .3 .3]);
    tiledlayout(1,4,'TileSpacing','compact');
    edges_dist = [0:0.2:10];
    nexttile;   histogram(FRpl_table.d1,edges_dist,'Normalization','cumcount','edgecolor','none','FaceColor','b','FaceAlpha',1); ylim('tight');  hold on; plot(0.6*[1 1],ylim,'b'); hold off; set(gca, 'XScale', 'log');    yticks([]); title('All');   xlabel('Distance from takeoff (m)');
    nexttile;   histogram(FRpl_table.d2,edges_dist,'Normalization','cumcount','edgecolor','none','FaceColor','b','FaceAlpha',1); ylim('tight');  hold on; plot(0.6*[1 1],ylim,'b'); hold off; set(gca, 'XScale', 'log');    yticks([]); title('All');   xlabel('Distance from landing (m)'); 
    nexttile;   histogram(RRpl_table.d1,edges_dist,'Normalization','cumcount','edgecolor','none','FaceColor','r','FaceAlpha',1); ylim('tight');  hold on; plot(0.6*[1 1],ylim,'b'); hold off; set(gca, 'XScale', 'log');    yticks([]); title('All');   xlabel('Distance from takeoff (m)'); 
    nexttile;   histogram(RRpl_table.d2,edges_dist,'Normalization','cumcount','edgecolor','none','FaceColor','r','FaceAlpha',1); ylim('tight');  hold on; plot(0.6*[1 1],ylim,'b'); hold off; set(gca, 'XScale', 'log');    yticks([]); title('All');   xlabel('Distance from landing (m)');   
 
    %=== Session and cluster averages
    FRpl_session = groupsummary(FRpl_table,'sessionID','mean',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','d1','d2','local1','local2','issleep','sleep_fraction'});
    RRpl_session = groupsummary(RRpl_table,'sessionID','mean',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','d1','d2','local1','local2','issleep','sleep_fraction'});
    FRpl_cluster = groupsummary(FRpl_table,'groupID','mean',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','d1','d2','local1','local2','issleep','sleep_fraction'});
    RRpl_cluster = groupsummary(RRpl_table,'groupID','mean',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','d1','d2','local1','local2','issleep','sleep_fraction'});

    %=== Display fraction of local replays
    figure('units','normalized','outerposition',[.3 .5 .15 .3]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v0(FRpl_session.mean_local1(FRpl_session.GroupCount>10), RRpl_session.mean_local2(RRpl_session.GroupCount>10), {'Forward', 'Reverse'}, 'SEM', 'Fraction Local Replays');  hold on; plot(xlim,[.5 .5],'k--');  ylim([0 1]);  
    nexttile;   plot_distr_AF_v0(FRpl_cluster.mean_local1(FRpl_cluster.GroupCount>10), RRpl_cluster.mean_local2(RRpl_cluster.GroupCount>10), {'Forward', 'Reverse'}, 'SEM', 'Fraction Local Replays');  hold on; plot(xlim,[.5 .5],'k--');  ylim([0 1]);  
%     signrank(FRpl_session.mean_local1(FRpl_session.GroupCount>50)-1)
%     title(signrank(FRpl_cluster.mean_local1-0.5));

    %=== Fraction of replays happening during sleep (AKA bouts in between flights)
    Fsa_ratio = ((FRpl_session.mean_issleep.*FRpl_session.GroupCount)./FRpl_session.mean_sleep_fraction)./(((1-FRpl_session.mean_issleep).*FRpl_session.GroupCount)./(1-FRpl_session.mean_sleep_fraction));
    Rsa_ratio = ((RRpl_session.mean_issleep.*RRpl_session.GroupCount)./RRpl_session.mean_sleep_fraction)./(((1-RRpl_session.mean_issleep).*RRpl_session.GroupCount)./(1-RRpl_session.mean_sleep_fraction));
    figure('units','normalized','outerposition',[.5 .5 .1 .3]);
    plot_distr_AF_v0(Fsa_ratio(Fsa_ratio<inf),Rsa_ratio(Rsa_ratio<inf), {'Forward', 'Reverse'}, 'SEM', 'Sleep/Awake Ratio');  hold on; plot(xlim,[1 1],'k--'); 
    
    %=== Fraction of replays happening close to tko in space and time (landing for reverse)
    FFraction_local = [];   RFraction_local = [];
    for i=1:10
        d_locRpl = 0.1*i;     t_locRpl = 6*i;
        FRpl_local = FRpl_table(FRpl_table.d1<d_locRpl & FRpl_table.t2tko<t_locRpl,:);  nFR_local = size(FRpl_local,1); FFraction_local= [FFraction_local;nFR_local/nFR];  
        RRpl_local = RRpl_table(RRpl_table.d2<d_locRpl & RRpl_table.t2lnd<t_locRpl,:);  nRR_local = size(RRpl_local,1); RFraction_local= [RFraction_local;nRR_local/nRR];
    end
    
    figure('units','normalized','outerposition',[.2 .2 .2 .6]);
    tiledlayout(2,2,'TileSpacing','compact');
    ax= nexttile;      radialBar(ax, FFraction_local, 1, [pi/2 -2*pi],1);
    ax= nexttile;      radialBar(ax, RFraction_local, 1, [pi/2 -2*pi],1);   
    nexttile;   pie([nFR_local/nFR,1-nFR_local/nFR]);   lgd = legend('Fraction Local','Location','East');    title([num2str(nFR_local), '/',num2str(nFR), ' Replays']);
    nexttile;   pie([nRR_local/nRR,1-nRR_local/nRR]);   lgd = legend('Fraction Local','Location','East');    title([num2str(nRR_local), '/',num2str(nRR), ' Replays']);
    
    
    %%
    %=== Plot Forward Replays
    n_event2show = min(200,nFR);
    cRT_subset = sort(datasample(1:nFR,n_event2show,'Replace',false));
    figure('units','normalized','outerposition',[0 0 .5 1]);
    tiledlayout(10,20,'TileSpacing','tight');
    for i=1:n_event2show
        s = FRpl_table.spikes{cRT_subset(i),1};
        N = size(s,1);  tgt_seq = [1:N]';
        tL = min(vertcat(s{:}));    tR = max(vertcat(s{:}));
        nexttile;   hold on;
        for nc= 1:N,plotSpikes_AF_v0(s{nc}, nc, 'r', 1);  end
        ylim([0 N+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[],[num2str((tR-tL)*1000,3)]});
    end
    
    %=== Plot Reverse Replays
    n_event2show = min(200,nRR);
    cRT_subset = sort(datasample(1:nRR,n_event2show,'Replace',false));
    figure('units','normalized','outerposition',[.5 0 .5 1]);
    tiledlayout(10,20,'TileSpacing','tight');
    for i=1:n_event2show
        s = RRpl_table.spikes{cRT_subset(i),1};
        N = size(s,1);  tgt_seq = [1:N]';
        tL = min(vertcat(s{:}));    tR = max(vertcat(s{:}));
        nexttile;   hold on;
        for nc= 1:N,plotSpikes_AF_v0(s{nc}, nc, 'k', 1);  end
        ylim([0 N+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[],[num2str((tR-tL)*1000,3)]});
    end
    
    
end

%% EXTRACT EVENTS OF INTEREST (FORWARD) AND LOOK AT SPEED
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 7;                      % Minimum number of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.4;                    % Minimum Correlation
    min_f = 0.3;                    % Minimum fraction of active cells
    max_t = 1;                      % Max duration
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                           % Id for each cluster from each session
    SB.GoodReplay = rowfun(@(t,x,y,z,w) (t>min_f) && (x>min_n) && (y>min_C) && (z<min_p) && (w<max_t), SB, 'InputVariables', {'f','n','C','p','dt'},'OutputFormat', 'uniform');    % Good Replay or not
       
    %=== Select good events
    Rpl_table = SB(SB.GoodReplay,:);
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','ascend');
    disp(['Average Replay Duration: ', num2str(median(Rpl_table.dt),3),' s']);
    disp(['Average Derived Replay Duration: ', num2str(median(Rpl_table.derived_dt),3),' s']);
    
    %=== Compare velocities for short vs long flights (all replays)
    Rpl_srt = Rpl_table(Rpl_table.flight_len>3 & Rpl_table.flight_len<7,:);
    Rpl_lng = Rpl_table(Rpl_table.flight_len>7 & Rpl_table.flight_len<13,:);
    
    figure('units','normalized','outerposition',[0.2 0.4 0.65 0.3]);
    tiledlayout(1,9,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v0(Rpl_srt.flight_len, Rpl_lng.flight_len, {'Short Flights', 'Long Flights'}, 'SEM', 'Flight Length (m)');
    nexttile;   plot_distr_AF_v0(Rpl_srt.derived_dt, Rpl_lng.derived_dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Corrected Duration(s)');
    nexttile;   plot_distr_AF_v0(Rpl_srt.dt, Rpl_lng.dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Duration(s)');
    nexttile;   plot_distr_AF_v0(Rpl_srt.C, Rpl_lng.C, {'Short Flights', 'Long Flights'}, 'SEM', 'Rank Correlation');
    nexttile;   plot_distr_AF_v0(Rpl_srt.f, Rpl_lng.f, {'Short Flights', 'Long Flights'}, 'SEM', 'Fraction Active');
    nexttile;   plot_distr_AF_v0(Rpl_srt.n, Rpl_lng.n, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Active');
    nexttile;   plot_distr_AF_v0(Rpl_srt.n_place_cells, Rpl_lng.n_place_cells, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Place');
    nexttile;   plot_distr_AF_v0(Rpl_srt.v, Rpl_lng.v, {'Short Flights', 'Long Flights'}, 'SEM', 'Velocity (m/s)');
    nexttile;   plot_distr_AF_v0(Rpl_srt.cvrg, Rpl_lng.cvrg, {'Short Flights', 'Long Flights'}, 'SEM', 'Coverage');
    
    %=== Compare velocities for short vs long flights (flight clusters)
    Rpl_Clus_srt = groupsummary(Rpl_srt,'groupID','median',{'flight_len','flight_dur','f','dt','v','c_ratio','derived_dt','C','n','cvrg','n_place_cells'});
    Rpl_Clus_lng = groupsummary(Rpl_lng,'groupID','median',{'flight_len','flight_dur','f','dt','v','c_ratio','derived_dt','C','n','cvrg','n_place_cells'});
    
    figure('units','normalized','outerposition',[0.2 0.1 0.65 0.3]);
    tiledlayout(1,9,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_flight_len, Rpl_Clus_lng.median_flight_len, {'Short Flights', 'Long Flights'}, 'SEM', 'Flight Length (m)');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_derived_dt, Rpl_Clus_lng.median_derived_dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Corrected Duration(s)');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_dt, Rpl_Clus_lng.median_dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Duration(s)');  ylim([0 max([Rpl_Clus_srt.median_dt;Rpl_Clus_lng.median_dt])]);
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_C, Rpl_Clus_lng.median_C, {'Short Flights', 'Long Flights'}, 'SEM', 'Rank Correlation');   ylim([0 1]);
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_f, Rpl_Clus_lng.median_f, {'Short Flights', 'Long Flights'}, 'SEM', 'Fraction Active');    ylim([0 1]);
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_n, Rpl_Clus_lng.median_n, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Active');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_n_place_cells, Rpl_Clus_lng.median_n_place_cells, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Place');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_v, Rpl_Clus_lng.median_v, {'Short Flights', 'Long Flights'}, 'SEM', 'Speed (m/s)');    ylim([0 max([Rpl_Clus_srt.median_v;Rpl_Clus_lng.median_v])]);
    nexttile;   plot_distr_AF_v0(Rpl_Clus_srt.median_cvrg, Rpl_Clus_lng.median_cvrg, {'Short Flights', 'Long Flights'}, 'SEM', 'Coverage');
    
    %=== Violin plot version
    figure('units','normalized','outerposition',[0.2 0.1 0.65 0.3]);
    tiledlayout(1,11,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_flight_len, Rpl_Clus_lng.median_flight_len, {'Short Flights', 'Long Flights'}, 'SEM', 'Flight Length (m)', 'violin');      currentYLim = ylim;  ylim([0, currentYLim(2)]); 
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_flight_dur, Rpl_Clus_lng.median_flight_dur, {'Short Flights', 'Long Flights'}, 'SEM', 'Flight Duration (s)', 'violin');    currentYLim = ylim;  ylim([0, currentYLim(2)]);    
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_derived_dt, Rpl_Clus_lng.median_derived_dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Corrected Duration(s)', 'violin');
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_c_ratio, Rpl_Clus_lng.median_c_ratio, {'Short Flights', 'Long Flights'}, 'SEM', 'Compression', 'violin');  
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_dt, Rpl_Clus_lng.median_dt, {'Short Flights', 'Long Flights'}, 'SEM', 'Duration(s)', 'violin');  
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_C, Rpl_Clus_lng.median_C, {'Short Flights', 'Long Flights'}, 'SEM', 'Rank Correlation', 'scatter');    currentYLim = ylim;  ylim([0, 1]);  
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_f, Rpl_Clus_lng.median_f, {'Short Flights', 'Long Flights'}, 'SEM', 'Fraction Active', 'scatter');      currentYLim = ylim;  ylim([0, 1]);
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_n, Rpl_Clus_lng.median_n, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Active', 'violin');
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_n_place_cells, Rpl_Clus_lng.median_n_place_cells, {'Short Flights', 'Long Flights'}, 'SEM', 'Number Place', 'violin');
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_v, Rpl_Clus_lng.median_v, {'Short Flights', 'Long Flights'}, 'SEM', 'Speed (m/s)', 'violin');    %ylim([0 max([Rpl_Clus_srt.median_v;Rpl_Clus_lng.median_v])]);
    nexttile;   plot_distr_AF_v1(Rpl_Clus_srt.median_cvrg, Rpl_Clus_lng.median_cvrg, {'Short Flights', 'Long Flights'}, 'SEM', 'Coverage', 'violin');
    sgtitle([mean(Rpl_Clus_srt.median_flight_len),mean(Rpl_Clus_lng.median_flight_len)]);
    
    mean(Rpl_Clus_srt.median_flight_dur)
    mean(Rpl_Clus_lng.median_flight_dur)
    
   %%
    
%     %=== Plot a few trajectory types
%     figure('units','normalized','outerposition',[0 0.4 1 0.3]);
%     tiledlayout(1,15,'TileSpacing','compact');
%     for jj=round(linspace(1,size(Rpl_table_ascending,1),15))
%         nexttile;
%         clus_path = Rpl_table_ascending.mean_path3D{jj,1};
%         plot3(clus_path(:,1),clus_path(:,2),clus_path(:,3),'-','LineWidth',1,'Color', 'k');
%         view(0,90);    axis equal;  xlim([-6 6]); ylim([-3 3]);   zlim([0 3]);  title([num2str(Rpl_table_ascending.flight_len(jj),3),' m']);
%     end

    %=== Summary statistics for replays from a given bat & cluster
    RPL_Clus = groupsummary(Rpl_table,'groupID','median',{'sessionID','flight_len','flight_dur','f','dt','v','c_ratio','derived_dt','C','n'});
    
    %=== Plot cumulative distribution for flight lengths and example trajectories
    figure('units','normalized','outerposition',[.2 0.2 .15 0.3]);
    plot(sort(RPL_Clus.median_flight_len)); ylim([0 max(RPL_Clus.median_flight_len)]);  xlabel('Flight Cluster Id');    ylabel('Lenght (m)');
    
    %=== Fit
    ft1 = fittype('a*x^2+b*x');
    fit1 = fit(RPL_Clus.median_derived_dt,RPL_Clus.median_flight_len,ft1,'start',[100 0],'Lower',[0,0],'Upper',[Inf,0.001]);
    
    %=== Plot summary statistics
    figure('units','normalized','outerposition',[.2 .3 .65 .3]);
    tiledlayout(1,8,'TileSpacing','compact');
%     nexttile;   scatter(Rpl_table.dt,Rpl_table.flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
%     xlabel('Replay Duration (s)'); ylabel('Trajectory lenght (m)'); xlim([0 1]);  ylim([0 max(Rpl_table.flight_len)]);
    nexttile;   scatter(RPL_Clus.median_v,RPL_Clus.median_flight_len,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Replay Speed (m/s)'); ylabel('Trajectory lenght (m)');lsline; xlim([0 max(RPL_Clus.median_v)]);  ylim([0 max(RPL_Clus.median_flight_len)]); axis square;
    nexttile;   scatter(RPL_Clus.median_dt,RPL_Clus.median_flight_len,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);  
    xlabel('Replay duration (s)'); ylabel('Trajectory lenght (m)'); xlim([0 1]);  ylim([0 max(RPL_Clus.median_flight_len)]);    axis square; 
    nexttile;   scatter(RPL_Clus.median_dt,RPL_Clus.median_flight_dur,20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);  
    xlabel('Replay duration (s)'); ylabel('Trajectory Duration (s)'); xlim([0 1]);  ylim([0 max(RPL_Clus.median_flight_dur)]);    axis square; 
    nexttile;   scatter(RPL_Clus.median_derived_dt,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);  hold on;
    %plot(fit1); legend('off');
    xlabel('Replay DRV. duration (s)'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_derived_dt)]);  ylim([0 max(RPL_Clus.median_flight_len)]); axis square;
    nexttile;   histogram(Rpl_table.v,[0:100],'Normalization','probability','facealpha',1,'edgecolor','none','FaceColor','k');
    xlabel('Replay Speed (m/s)');   ylabel('Fraction'); axis square;
    nexttile;   scatter(RPL_Clus.median_C,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Rank-Correlation'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_C)]);  ylim([0 max(RPL_Clus.median_flight_len)]);  axis square;
    nexttile;   scatter(RPL_Clus.median_n,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('N. Active Cells'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_n)]);  ylim([0 max(RPL_Clus.median_flight_len)]);   axis square;
    nexttile;   scatter(RPL_Clus.median_f,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Fraction Active Cells'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_f)]);  ylim([0 max(RPL_Clus.median_flight_len)]);   axis square;
    
    [c_tmp,p_tmp] = corr(RPL_Clus.median_dt,RPL_Clus.median_flight_len)
    
    %%
    
    min_n = 5;                      % Minimum number of active cells
    min_p = 0.1;                   % Minimum p value
    min_C = 0.4;                    % Minimum Correlation
    min_f = 0.3;                    % Minimum fraction of active cells
    max_t = 1;                      % Max duration
    
    %=== Add a few features
    SB.GoodReplay = rowfun(@(t,x,y,z,w) (t>min_f) && (x>min_n) && (y>min_C) && (z<min_p) && (w<max_t), SB, 'InputVariables', {'f','n','C','p','dt'},'OutputFormat', 'uniform');    % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v; % Derived duration
    SB.compression = SB.flight_dur./SB.dt;
    Rpl_table = SB(SB.GoodReplay,:);
    
    %=== Plot a few examples, sorted by trajectory lenght and using the same timescale
    %=== Plot replays
    rng(5);
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','descend');
    cRT_good = find(Rpl_table_ascending.n>7 & Rpl_table_ascending.p<0.05 & Rpl_table_ascending.C>0.4 & Rpl_table_ascending.f>0.3);
    n_event2show = min(13*17,size(cRT_good,1));
    cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
    cRT_short = find(Rpl_table_ascending.n>6 & Rpl_table_ascending.p<0.1 & Rpl_table_ascending.C>0.4 & Rpl_table_ascending.f>0.3 & Rpl_table_ascending.flight_len<3);
    cRT_subset(end-numel(cRT_short)+1:end) = cRT_short;
    Rpl_subtable = Rpl_table_ascending(cRT_good,:);
    figure('units','normalized','outerposition',[.2 0 .6 1]);
    tiledlayout(13,17,'TileSpacing','compact');
    for i=1:n_event2show
        s = Rpl_table_ascending.spikes{cRT_subset(i),1};
        N = size(s,1);              % Total number of cells
        tgt_seq = [1:N]';           % Target sequence of activation
        
        %=== Left and right extremes
        tL = min(vertcat(s{:}));
        tR = max(vertcat(s{:}));
        nexttile;   hold on; 
        for nc= 1:N
            plotSpikes_AF_v0(s{nc}, nc, 'r', 1);              
        end
        ylim([0 N+1]); xticks([]); xlim([tL, tL+max(Rpl_table_ascending.dt)]); yticks([]);
        xlabel([num2str(Rpl_table_ascending.flight_len(cRT_subset(i)),2),' m'])
    end
    
    %%
    %=== Plot replays (vertical version)
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','descend');
    counter =1;
    wR = 0.07;  n_events2disp = 13;
    cRT_good = find(Rpl_table_ascending.n>6 & Rpl_table_ascending.p<0.1 & Rpl_table_ascending.C>0.6 & Rpl_table_ascending.f>0.1 & Rpl_table_ascending.dt<0.5);
    Rpl_subtable = Rpl_table_ascending(cRT_good,:);
    n_event2show = min(n_events2disp,size(cRT_good,1));
    for sd = round(linspace(1,15,15))
        
        rng(sd);
        cRT_subset = sort(datasample(1:size(Rpl_subtable,1),n_event2show,'Replace',false));
        
        figure('units','normalized','outerposition',[0.9*wR*(counter-1) 0 wR 1]);
        tiledlayout(n_events2disp,1,'TileSpacing','compact');
        for i=1:n_event2show
            s = Rpl_subtable.spikes{cRT_subset(i),1};
            N = size(s,1);              % Total number of cells
            tgt_seq = [1:N]';           % Target sequence of activation
            
            %=== Left and right extremes
            tL = min(vertcat(s{:}));
            tR = max(vertcat(s{:}));
            nexttile;   hold on;
            for nc= 1:N
                plotSpikes_AF_v0(s{nc}, nc, 'r', 1);             
            end
            ylim([0 N+1]); xticks([]); xlim([tL, tL+max(Rpl_subtable.dt)]); yticks([]);
            legend([num2str(Rpl_subtable.flight_len(cRT_subset(i)),2),' m'],'Location','eastoutside');
        end
        sgtitle(['Seed: ', num2str(sd)]);
        counter=counter+1;
    end
    
   %%
    %=== Pairwise timings vs field distances
    %Rpl_table_sst = Rpl_table(Rpl_table.n>7 & Rpl_table.p<0.01 & Rpl_table.C>0.6,:);
    Rpl_table_sst = Rpl_table(Rpl_table.n>7 & Rpl_table.p<0.05 & Rpl_table.C>0.4 & Rpl_table.f>0.3,:);
    pw_timing = vertcat(Rpl_table_sst.pw_timing{:});
    pw_distance = vertcat(Rpl_table_sst.pw_distance{:});
    pw_timing(pw_distance<0) = -pw_timing(pw_distance<0);
    pw_distance = abs(pw_distance);
    [pw_distance,id_tmp] = sort(pw_distance);
    pw_timing = pw_timing(id_tmp);
    
    %=== Subsample distances in a balanced way
    pwD_cell = {};
    pwT_cell = {};
    [D_counts,~,D_bin_id] = histcounts(pw_distance,[0:2:11]);
    bin_ids = setdiff(unique(D_bin_id),0);
    for i=bin_ids'
        pwD_cell{i} = pw_distance(D_bin_id==i);
        pwT_cell{i} = pw_timing(D_bin_id==i);
    end
    min_size = min(cellfun(@numel,pwD_cell));
    cell_sst = cellfun(@(x) datasample(1:numel(x),min_size,'Replace',false),pwD_cell,'UniformOutput',false);
    pwD_cell_sst = cellfun(@(x,y) x(y),pwD_cell,cell_sst,'UniformOutput',false);
    pwT_cell_sst = cellfun(@(x,y) x(y),pwT_cell,cell_sst,'UniformOutput',false);
    xData = vertcat(pwD_cell_sst{:});
    yData = vertcat(pwT_cell_sst{:});
    goodEntries = ~isnan(yData);
    xData = xData(goodEntries);
    yData = yData(goodEntries);
    
    % Define the fitting model
    model = fittype('a*(1-exp(-b*x))', 'independent', 'x', 'coefficients', {'a', 'b'});
    fit_result = fit(xData, yData, model, 'StartPoint', [1, 1],'Exclude', yData < 0.00);
    % model = fittype('a*x','independent', 'x', 'coefficients', {'a'});
    % fit_result = fit(xData, yData, model, 'StartPoint', 1,'Exclude', abs(xData) > 5);
    
    figure('units','normalized','outerposition',[.2 0.2 .15 0.5]);
    scatter(xData,yData,40,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.2);   set(gca,'YScale','log');
    hold on;plot(fit_result);   legend('off');  ylim(prctile(yData,[5 99]));
    
end

%% LOOK AT AVERAGE DECODED SHAPE OF REPLAYS AND FLIGHT TYPE ANALYSIS
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 7;                      % Minimum number of active cells (6)
    min_f = 0.3;                    % Minimum fraction of active cells (0.1)
    min_p = 0.05;                   % Minimum p value (0.05)
    min_C = 0.4;                    % Minimum Correlation (0.6)
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                                               % Id for each cluster from each session
    SB.GoodReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z> min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');    % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                                        % Derived duration
    SB.compression = SB.flight_dur./SB.dt;
    SB.cvrg =  cellfun(@(x) x(end)-x(1),SB.place_phases);   % Flight phase of last minus first cell
   
    %=== Select good events
    Rpl_table = SB(SB.GoodReplay,:);
    n_good = size(Rpl_table,1);
    
    %=== Calculate weighted correlation and percentage of decoded bins
    for i=1:n_good
        
        posterior = Rpl_table.posterior{i,1}(3:end-2,:);                % Get the posterior probability (cut the edges)
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial');
        [max_P,max_L] = max(pstP);
        invalidP = (max_P==0 | isnan(max_P));                           % Invalidate NaNs or zeros
        
       Rpl_table.fract_good(i) = 1-sum(sum(posterior,1)<0.99 | sum(posterior)==NaN)/nT;
        
        if ~all(invalidP)
            max_P(invalidP)=[]; max_L(invalidP)=[];
            corrmat = weightedcorrs([[1:numel(max_L)]',max_L'], max_P');
            Rpl_table.weightcorr(i) = corrmat(1,2);
        else
            Rpl_table.weightcorr(i) = NaN;
        end
        
    end
    
    %=== Select good events
    Rpl_table = Rpl_table(~isnan(Rpl_table.weightcorr),:);
    n_good = size(Rpl_table,1);
    
    %=== 'Flight type' analysis
    [~,tmp_idx,tmp_idx1] = unique(Rpl_table.groupID);
    num_clus = numel(tmp_idx);
    med_curv = NaN(num_clus,1);
    gaussR2 = NaN(num_clus,1);
    mean_v1D_rsz = NaN(10,num_clus);
    for jj = 1:num_clus
        
        %=== Get mean path 3D
        mean_path3D = Rpl_table.mean_path3D{tmp_idx(jj),1};
        
        %=== Calculate curvature
        norm_grad3D = diff(mean_path3D,1,1)./vecnorm(diff(mean_path3D,1,1),2,2);
        curv3D = vecnorm(diff(norm_grad3D,1,1),2,2);
        [~,R_tmp,~] = curvature(mean_path3D); % from Are Mjaavatten
        curv3D_bis = 1./R_tmp;
        med_curv(jj) = median(curv3D,'omitnan');
        
        %=== Get mean paths 1D
        mean_path1D = Rpl_table.mean_path1D{tmp_idx(jj),1};
        mean_time1D = Rpl_table.mean_time1D{tmp_idx(jj),1};
        
        %=== Calculate fit R2 with 1 or 2 gaussians
        mean_v1D = diff(mean_path1D)./diff(mean_time1D);
        [fg1,gof1] = fit(mean_time1D(2:end)',mean_v1D','gauss1','StartPoint',[6 1 1]);
        [fg2,gof2] = fit(mean_time1D(2:end)',mean_v1D','gauss2','StartPoint',[6 1 1 6 1 1]);
        gaussR2(jj) = gof1.rsquare;
        
        %=== Normalize velocity profile for dim-reduction
        %mean_v1D_rsz(:,jj) = interp1(mean_v1D./max(mean_v1D),linspace(1,numel(mean_v1D),size(mean_v1D_rsz,1)),'linear')';
        mean_v1D_rsz(:,jj) = interp1(mean_v1D./max(mean_v1D),linspace(1,numel(mean_v1D),10),'linear')';
        
    end
    Rpl_table.med_curv =  med_curv(tmp_idx1);   % Assign to Rpl table entries
    Rpl_table.gaussR2 =  gaussR2(tmp_idx1);     % Assign to Rpl table entries
    
    %=== Cluster flight types
    for i=5
    rng(i); warning('off');
    flight_type = kmeans([med_curv,gaussR2],2);
    [~,score,~] = pca(mean_v1D_rsz');
    tSNE = tsne(mean_v1D_rsz','Distance','chebychev','Exaggeration',5);
    figure('units','normalized','outerposition',[.2 .3 .2 .25]);
    tiledlayout(1,3,'TileSpacing','compact');
    nexttile;   gscatter(med_curv,gaussR2,flight_type);     legend('off');  %set(gca,'XScale','log');
    axis square;    xlabel('Median Curvature'); ylabel('R square (fit)');
    nexttile;   gscatter(score(:,1),score(:,2),flight_type); legend('off'); %set(gca,'XScale','log');
    axis square;    xlabel('PCA 1'); ylabel('PCA 2');
    nexttile;   gscatter(tSNE(:,1),tSNE(:,2),flight_type); legend('off');   %set(gca,'XScale','log');
    axis square;    xlabel('t-SNE 1'); ylabel('t-SNE 2');   xlim(xlim*2); ylim(ylim*2);
    warning('on');
    end
    
    %=== Assign cluster ID, resize flight shape and store speed around the center of the flight
    Rpl_table.flight_type = double(Rpl_table.gaussR2<0.6)+1;
    Rpl_table.rsz_mean_path1D = cellfun(@(x) interp1(x./max(x),linspace(1,numel(x),100)),Rpl_table.mean_path1D,'UniformOutput',false);
    Rpl_table.rsz_mean_velo1D = cellfun(@(x,y) interp1(diff(x)./diff(y),linspace(1,numel(x)-1,100)),Rpl_table.mean_path1D,Rpl_table.mean_time1D,'UniformOutput',false);
    Rpl_table.local_speed1D = cellfun(@(x) mean(x(30:70)),Rpl_table.rsz_mean_velo1D);
    
    %=== Quantify the slope of replay close to the center of the event
    %figure('units','normalized','outerposition',[.3 .6 .15 .3]);
    logisEqn = '1/(1+exp(-b*(x-c)))';
    lineEqn = 'm*x+q';
    smooth_f = [.01 1];
    n_nbh = 10;
    for i=1:n_good
        
        posterior = Rpl_table.posterior{i,1};                           % Get the posterior probability
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        spt_bin_ids = [1:nS];                                           % Spatial bin ids
        invalid_bins = (sum(posterior)<0.99 | isnan(sum(posterior)));   % 1 where invalid bin
        posterior(:,invalid_bins) = zeros(nS,sum(invalid_bins));        % Make sure invalid bins are at zero
        pstP = imgaussfilt(posterior,smooth_f,'FilterDomain','spatial');% Smooth the posterior
        [~,cnt_mass] = max(pstP);                                       % Calculate the max location
        %cnt_mass = spt_bin_ids*pstP;                                    % Calculate the center of mass
        cnt_mass(invalid_bins) = NaN;                                   % NaN invalid bins
        cnt_mass(1) = 1;    cnt_mass(end) = nS;                         % Force the tails 
        cnt_mass = fillmissing(cnt_mass,'linear');                      % Linearly interpolate missing entries
        y1 = normalize(cnt_mass,'range')';                              % Normalize between 0 and 1 
        x1 = [1:nT]';                                                   % Time bins  
        w = max(posterior); w(1)=1; w(end)=1;                           % Weight is proportional to the posterior
        [f,gof] = fit(x1,y1,logisEqn,'Start',[3 nT/2],'Weights',w);     % Fit with logistic function
        fitcoeff = coeffvalues(f);                                      % Get the coefficients 
        Rpl_table.central_bin(i) = fitcoeff(2);                         % Get the center of the logistics f
        Rpl_table.slope(i) = fitcoeff(1);                               % Get the slope of the logistics f
        Rpl_table.lrs(i) = gof.rsquare;                                 % Get the R square
        
        %         imagesc([],[0 1],posterior);     hold on; plot(f,x1,y1);
        %         colormap(flipud(gray));  set(gca,'YDir','normal');  legend('off');  axis('off');
        %         title(gof.rsquare);
        %
        %         choice = questdlg('Is this image of the correct class?', 'Class Verification','Yes', 'No', 'Stop','Yes');
        %         switch choice
        %             case 'Stop'
        %                 break;
        %         end
        
        %=== Fit slope around the center of the replay when possible
        if round(fitcoeff(2))-n_nbh>0 && round(fitcoeff(2))+n_nbh<nT
            y1 = cnt_mass(round(fitcoeff(2))+[-n_nbh:+n_nbh])';
            x1 = [1:numel(y1)]';
            w = max(posterior(:,round(fitcoeff(2))+[-n_nbh:+n_nbh]));
            [f,gof] = fit(x1,y1,lineEqn,'Start',[1 y1(1)],'Weights',w);
            fitcoeff = coeffvalues(f);
            Rpl_table.local_slope(i) = fitcoeff(1)*(bin_size_1D/Rpl_table.bin1Dt(1));
        
        else
            Rpl_table.local_slope(i) = NaN;
        end

    end

    %%
    %=== Adjust slope for temporal compression factor (Adjusted Local Slope)
    Rpl_table.als1 = Rpl_table.local_slope./Rpl_table.c_ratio;
    Rpl_table.als2 = Rpl_table.local_slope./Rpl_table.compression;
    
    %=== Restrict analysis to good events
    RPL_sst = Rpl_table(Rpl_table.weightcorr>0.5 & Rpl_table.fract_good>0.5 & Rpl_table.lrs>0.1,:);
    RPL_sst1 = RPL_sst(RPL_sst.flight_type == 1,:);
    RPL_sst2 = RPL_sst(RPL_sst.flight_type == 2,:);
    
    %=== Look at single replays
    figure('units','normalized','outerposition',[0.2 0.4 0.5 0.3]);
    tiledlayout(1,9,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v0(RPL_sst1.local_slope, RPL_sst2.local_slope, {'Type 1', 'Type 2'}, 'SEM', 'Local Slope');
    nexttile;   plot_distr_AF_v0(RPL_sst1.als1, RPL_sst2.als1, {'Type 1', 'Type 2'}, 'SEM', 'Adjuste Local Slope 1');
    nexttile;   plot_distr_AF_v0(RPL_sst1.als2, RPL_sst2.als2, {'Type 1', 'Type 2'}, 'SEM', 'Adjuste Local Slope 2');
    nexttile;   plot_distr_AF_v0(RPL_sst1.dt, RPL_sst2.dt, {'Type 1', 'Type 2'}, 'SEM', 'Duration');
    nexttile;   plot_distr_AF_v0(RPL_sst1.weightcorr, RPL_sst2.weightcorr, {'Type 1', 'Type 2'}, 'SEM', 'W. Corr');
    nexttile;   plot_distr_AF_v0(RPL_sst1.fract_good, RPL_sst2.fract_good, {'Type 1', 'Type 2'}, 'SEM', 'Fraction Good');
    nexttile;   plot_distr_AF_v0(RPL_sst1.c_ratio, RPL_sst2.c_ratio, {'Type 1', 'Type 2'}, 'SEM', 'C Ratio');
    nexttile;   plot_distr_AF_v0(RPL_sst1.compression, RPL_sst2.compression, {'Type 1', 'Type 2'}, 'SEM', 'Compression');
    nexttile;   plot_distr_AF_v0(RPL_sst1.local_speed1D, RPL_sst2.local_speed1D, {'Type 1', 'Type 2'}, 'SEM', 'Flight local speed');

    %=== Look at averages across flight clusters
    Rpl_Clus_sst1 = groupsummary(RPL_sst1,'groupID','median',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','local_slope','als1','als2','weightcorr','fract_good','compression','local_speed1D'});
    Rpl_Clus_sst2 = groupsummary(RPL_sst2,'groupID','median',{'flight_len','f','dt','v','c_ratio','derived_dt','C','n','cvrg','local_slope','als1','als2','weightcorr','fract_good','compression','local_speed1D'});
    
    figure('units','normalized','outerposition',[0.2 0.1 0.5 0.3]);
    tiledlayout(1,9,'TileSpacing','compact');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_local_slope, Rpl_Clus_sst2.median_local_slope, {'Type 1', 'Type 2'}, 'SEM', 'Local Slope');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_als1, Rpl_Clus_sst2.median_als1, {'Type 1', 'Type 2'}, 'SEM', 'Adjuste Local Slope 1');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_als2, Rpl_Clus_sst2.median_als2, {'Type 1', 'Type 2'}, 'SEM', 'Adjuste Local Slope 2');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_dt, Rpl_Clus_sst2.median_dt, {'Type 1', 'Type 2'}, 'SEM', 'Duration');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_weightcorr, Rpl_Clus_sst2.median_weightcorr, {'Type 1', 'Type 2'}, 'SEM', 'W. Corr');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_fract_good, Rpl_Clus_sst2.median_fract_good, {'Type 1', 'Type 2'}, 'SEM', 'Fraction Good');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_c_ratio, Rpl_Clus_sst2.median_c_ratio, {'Type 1', 'Type 2'}, 'SEM', 'C Ratio');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_compression, Rpl_Clus_sst2.median_compression, {'Type 1', 'Type 2'}, 'SEM', 'Compression');
    nexttile;   plot_distr_AF_v0(Rpl_Clus_sst1.median_local_speed1D, Rpl_Clus_sst2.median_local_speed1D, {'Type 1', 'Type 2'}, 'SEM', 'Flight local speed');
    sgtitle(['N clusters: ',num2str(size(Rpl_Clus_sst1,1)),'-',num2str(size(Rpl_Clus_sst2,1))]);

    %=== Look at average distance and speed profiles
    [~,tmp_idx1,~] = unique(RPL_sst1.groupID);
    [~,tmp_idx2,~] = unique(RPL_sst2.groupID);
    figure('units','normalized','outerposition',[.2 .3 .2 .25]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   
    data_tmp = vertcat(RPL_sst1.rsz_mean_path1D{tmp_idx1,1});
    plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');   hold on;
    data_tmp = vertcat(RPL_sst2.rsz_mean_path1D{tmp_idx2,1});
    plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');
    xlabel('Flight %'); ylabel('Normalized Distance');
    nexttile;   
    data_tmp = vertcat(RPL_sst1.rsz_mean_velo1D{tmp_idx1,1});
    plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');   hold on;
    data_tmp = vertcat(RPL_sst2.rsz_mean_velo1D{tmp_idx2,1});
    plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');
    xlabel('Flight %'); ylabel('Bat Speed (m/s)');
  
%     %=== Plot some example replays for straight and loop flights
%     n_evt = 50;
%     ssT1 = datasample(1:size(RPL_sst1,1),n_evt,'Replace',false);
%     ssT2 = datasample(1:size(RPL_sst2,1),n_evt,'Replace',false);
%     figure('units','normalized','outerposition',[0 .3 1 0.3]);
%     tiledlayout(2,50,'TileSpacing','compact');
%     for i=ssT1
%         pstP = imgaussfilt(RPL_sst1.posterior{i,1},[.01 1],'FilterDomain','spatial'); % Filter
%         nexttile;   imagesc(pstP,prctile(pstP,[1 99],'all')')
%         colormap('hot');  set(gca,'YDir','normal');                                     xticks([]); yticks([]);
%     end
%     for i=ssT2
%         pstP = imgaussfilt(RPL_sst2.posterior{i,1},[.01 1],'FilterDomain','spatial'); % Filter
%         nexttile;   imagesc(pstP,prctile(pstP,[1 99],'all')')
%         colormap('hot');  set(gca,'YDir','normal');                                     xticks([]); yticks([]);
%     end
    
end
