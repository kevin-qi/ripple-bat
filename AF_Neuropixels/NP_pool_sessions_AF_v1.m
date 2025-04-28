%% POOL REPLAYS AND QUANTIFY INTERESTING FEATURES

%=== Load Data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Sessions to load (comment sessions you don't want to include)
    sessions2include = {...
        'Dataset_1','32622','231006';...
        'Dataset_1','32622','231007';...
        'Dataset_1','32622','231008';...
        'Dataset_1','32622','231009';...
        'Dataset_1','32622','231010';...
        'Dataset_2','14445','231208';...
        'Dataset_2','14445','231209';...
        'Dataset_2','14611','231208';...
        %'Dataset_2','14611','231209';...
        'Dataset_2','14611','231210';...
        'Dataset_2','14611','231211';...
        'Dataset_2','14611','231212';...
        'Dataset_3','00000','240402';...
        'Dataset_3','00000','240403';...
        'Dataset_3','00000','240404';...
        'Dataset_3','00000','240405';...
        %'Dataset_3','00000','240406';...
        %'Dataset_3','00000','240407';...
        'Dataset_3','14543','240419';...
        'Dataset_3','14543','240420';...
        'Dataset_3','14543','240421';...
        %'Dataset_3','14543','240422';...
        'Dataset_3','29959','240402';...
        'Dataset_3','29959','240403';...
        'Dataset_3','29959','240404';...
        'Dataset_3','29959','240405';...
        %'Dataset_3','29959','240406';...
        %'Dataset_3','29959','240407';...
        ...
        };
    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   SB = []; ids = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(Rpl_table) && any(cellfun(@(row) isequal(row, Rpl_table.unique_ID(1,:)), num2cell(sessions2include, 2)))
            SB = [SB; Rpl_table];
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except SB
end

%% EXTRACT EVENTS OF INTEREST (FORWARD) AND LOOK AT SPEED
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 5;                      % Minimum number of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.7;                    % Minimum Correlation
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                           % Id for each cluster from each session
    SB.GoodReplay = rowfun(@(x,y,z) (x>min_n) && (y>min_C) && (z<min_p), SB, 'InputVariables', {'n','C','p'},'OutputFormat', 'uniform');    % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                    % Derived duration
    
    %=== Select good events
    Rpl_table = SB(SB.GoodReplay,:);
    
    %=== Summary statistics for replays from a given bat & cluster
    RPL_Clus = groupsummary(Rpl_table,'groupID','median',{'flight_len','dt','v','c_ratio','derived_dt','C','n'});
    
    %=== Fit
    ft1 = fittype('a*x^2+b*x');
    fit1 = fit(RPL_Clus.median_derived_dt,RPL_Clus.median_flight_len,ft1,'start',[100 0],'Lower',[0,0],'Upper',[Inf,0.001]);
    
    %=== Plot summary statistics
    figure('units','normalized','outerposition',[.3 .3 .5 .3]);
    tiledlayout(1,5,'TileSpacing','compact');
    nexttile;   scatter(RPL_Clus.median_v,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Replay Speed (m/s)'); ylabel('Trajectory lenght (m)');lsline; xlim([0 max(RPL_Clus.median_v)]);  ylim([0 max(RPL_Clus.median_flight_len)]);
    nexttile;   scatter(RPL_Clus.median_derived_dt,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);  hold on;
    plot(fit1); legend('off');
    xlabel('Replay duration (s)'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_derived_dt)]);  ylim([0 max(RPL_Clus.median_flight_len)]);
    nexttile;   histogram(Rpl_table.v,[0:100],'Normalization','probability','facealpha',1,'edgecolor','none','FaceColor','k');
    xlabel('Replay Speed (m/s)');   ylabel('Fraction');
    nexttile;   scatter(RPL_Clus.median_C,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Rank-Correlation'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_C)]);  ylim([0 max(RPL_Clus.median_flight_len)]);
    nexttile;   scatter(RPL_Clus.median_n,RPL_Clus.median_flight_len,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('N. Active Cells'); ylabel('Trajectory lenght (m)'); xlim([0 max(RPL_Clus.median_n)]);  ylim([0 max(RPL_Clus.median_flight_len)]);
    
    %=== Plot a few examples, sorted by trajectory lenght and using the same timescale
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','ascend');
    cRT_good = find(Rpl_table_ascending.n>5 & Rpl_table_ascending.p<0.05 & Rpl_table_ascending.C>0.7);
    n_event2show = min(150,size(cRT_good,1));
    cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
    Rpl_subtable = Rpl_table_ascending(cRT_good,:);
    
    %=== Plot replays
    figure('units','normalized','outerposition',[.2 0 .6 1]);
    tiledlayout(10,15,'TileSpacing','tight');
    for i=1:n_event2show
        s = Rpl_subtable.spikes{cRT_subset(i),1};
        N = size(s,1);              % Total number of cells
        tgt_seq = [1:N]';           % Target sequence of activation
        
        %=== Left and right extremes
        tL = min(vertcat(s{:}));
        tR = max(vertcat(s{:}));
        nexttile;
        for nc= 1:N
            plot(s{nc}, nc*ones(size(s{nc})), 'r|','MarkerSize', 5);  hold on;             % Raster for each cell
        end
        ylim([0 N+1]); xticks([]); xlim([tL, tL+max(Rpl_subtable.dt)]); yticks([]);
        xlabel([num2str(Rpl_subtable.flight_len(cRT_subset(i)),2),' m'])
    end
    
    
    %=== Plot all linearized flight paths
    figure('units','normalized','outerposition',[.2 0.2 .13 0.3]);
    col = cool(size(Rpl_table,1));
    for i=flip(1:size(Rpl_table_ascending,1))
        plot(Rpl_table_ascending.mean_time1D{i,1},Rpl_table_ascending.mean_path1D{i,1},'Color',col(i,:));    hold on;
    end
    xlabel('Time(s)');  ylabel('Trajectory length (m)');
    
    
    
    %=== Plot pairwise timings vs field distances
    figure('units','normalized','outerposition',[.2 0.2 .7 0.3]);
    tiledlayout(1,6,'TileSpacing','tight');
    for i=1:6
        Rpl_table_sst = Rpl_table(Rpl_table.n>4+i & Rpl_table.p<0.05 & Rpl_table.C>0.7,:);
        pw_timing = vertcat(Rpl_table_sst.pw_timing{:});
        pw_distance = vertcat(Rpl_table_sst.pw_distance{:});
        pw_timing(pw_distance<0) = -pw_timing(pw_distance<0);
        %subset = datasample(1:numel(pw_timing),100000,'Replace',false);
        subset = [1:numel(pw_timing)];
        ax(i) = nexttile;
        scatter(abs(pw_distance(subset)),pw_timing(subset),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.05); hold on;
        
        x_data = abs(pw_distance(subset));
        y_data = pw_timing(subset);
        to_keep = x_data<6 & ~isnan(x_data) & ~isnan(y_data);
        y_data = y_data(to_keep);
        x_data = x_data(to_keep);
        p_fit = polyfit(x_data,y_data,1);
        y_fit = polyval(p_fit,[x_data;20]);
        plot([x_data;20],y_fit,'k--');
        ylabel('dt during F replay (s)'); xlabel('Place field distance (m)');
        title(['Min ',num2str(4+i),' active cells'])
    end
    linkaxes(ax);
    
    %=== Plot pairwise timings vs field distances
    figure('units','normalized','outerposition',[.2 0.2 .17 0.3]);
    col = cool(2);
    Rpl_table_sst = Rpl_table(Rpl_table.n>5 & Rpl_table.p<0.05 & Rpl_table.C>0.7,:);
    pw_timing = vertcat(Rpl_table_sst.pw_timing{:});
    pw_distance = vertcat(Rpl_table_sst.pw_distance{:});
    pw_timing(pw_distance<0) = -pw_timing(pw_distance<0);
    [pw_distance,id_tmp] = sort(pw_distance);
    pw_timing = pw_timing(id_tmp);
    subset = datasample(1:numel(pw_timing),10000,'Replace',false);
    %subset = [1:numel(pw_timing)];
    scatter(abs(pw_distance(subset)),pw_timing(subset),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.3); hold on;
    d_thresh = [5 20];
    for i=1:2
        x_data = abs(pw_distance(subset));
        y_data = pw_timing(subset);
        to_keep = x_data<d_thresh(i) & ~isnan(x_data) & ~isnan(y_data);
        y_data = y_data(to_keep);
        x_data = x_data(to_keep);
        p_fit = polyfit(x_data,y_data,1);
        y_fit = polyval(p_fit,[x_data;20]);
        plot([x_data;20],y_fit,'Color',col(i,:));
    end
    colormap('cool');colorbar;
    title(['Min ',num2str(5),' active cells'])
    ylabel('dt during F replay (s)'); xlabel('Place field distance (m)');
    
end

%% LOOK AT FORWARD vs REVERSE REPLAY FEATURES
for hide=1
    
    %=== Params
    min_n = 5;                      % Minimum number of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.7;                    % Minimum Correlation
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                           % Id for each cluster from each session
    SB.GoodFReplay = rowfun(@(x,y,z) (x>min_n) && (y> min_C) && (z<min_p), SB, 'InputVariables', {'n','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.GoodRReplay = rowfun(@(x,y,z) (x>min_n) && (y<-min_C) && (z<min_p), SB, 'InputVariables', {'n','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                    % Derived duration
    
    %=== Select good events
    FRpl_table = SB(SB.GoodFReplay & SB.bat_on_feeder,:);  nFR = size(FRpl_table,1);  % Forward Replays
    RRpl_table = SB(SB.GoodRReplay & SB.bat_on_feeder,:);  nRR = size(RRpl_table,1);  % Reverse Replays
    
    %=== Plot Forward Replays
    n_event2show = min(200,nFR);
    cRT_subset = sort(datasample(1:nFR,n_event2show,'Replace',false));
    figure('units','normalized','outerposition',[0 0 .5 1]);
    tiledlayout(10,20,'TileSpacing','tight');
    for i=1:n_event2show
        s = FRpl_table.spikes{cRT_subset(i),1};
        N = size(s,1);  tgt_seq = [1:N]';
        tL = min(vertcat(s{:}));    tR = max(vertcat(s{:}));
        nexttile;
        for nc= 1:N,plot(s{nc}, nc*ones(size(s{nc})), 'k|','MarkerSize', 5);  hold on;end
        ylim([0 N+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[],[num2str((tR-tL)*1000,3), ' ms']});
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
        nexttile;
        for nc= 1:N,plot(s{nc}, nc*ones(size(s{nc})), 'r|','MarkerSize', 5);  hold on;end
        ylim([0 N+1]); xticks([tL, tR]); xlim([tL, tR]); yticks([]); xticklabels({[],[num2str((tR-tL)*1000,3), ' ms']});
    end
    
    %=== Summary Statistics
    figure('units','normalized','outerposition',[.3 .3 .5 .27]);
    tiledlayout(1,6,'TileSpacing','compact');
    nexttile;
    histogram(FRpl_table.C,'Normalization','count','facealpha',.5,'edgecolor','none','FaceColor','k');  hold on;
    histogram(RRpl_table.C,'Normalization','count','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('Rank-Correlation');   ylabel('Counts');
    nexttile;
    histogram(FRpl_table.bat_on_feeder,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');   hold on;
    histogram(RRpl_table.bat_on_feeder,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('On Feeder');   ylabel('Fraction');
    nexttile;
    histogram(FRpl_table.bat_mdg,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');     hold on;
    histogram(RRpl_table.bat_mdg,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('Accelerometer');   ylabel('Fraction');
    nexttile;
    histogram(FRpl_table.v,[0 :100],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');   hold on;
    histogram(RRpl_table.v,[-100:0],'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('Replay Speed (m/s)');   ylabel('Fraction');
    nexttile;
    histogram(FRpl_table.d1,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');   hold on;
    histogram(RRpl_table.d1,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('Distance from Takeoff (m)');   ylabel('Fraction');
    nexttile;
    histogram(FRpl_table.d2,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','k');   hold on;
    histogram(RRpl_table.d2,'Normalization','probability','facealpha',.5,'edgecolor','none','FaceColor','r');
    xlabel('Distance from Landing (m)');   ylabel('Fraction');
    
end
