%% POOL REPLAYS AND QUANTIFY INTERESTING FEATURES

%=== Load Data
for hide=1
    clr;
    
    %=== Default settings
    set(groot,'defaultAxesFontSize',12);
    set(groot,'defaultHistogramEdgeColor','none','defaultHistogramFaceAlpha',0.5);
    
    %=== Params
    min_n = 5;                      % Minimum number of active cells
    min_f = 0.2;                    % Minimum fraction of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.1;                    % Minimum Correlation
    
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
        'Dataset_2','14611','231209';...
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
        'Dataset_3','14543','240422';...
        'Dataset_3','29959','240402';...
        'Dataset_3','29959','240403';...
        'Dataset_3','29959','240404';...
        'Dataset_3','29959','240405';...
        %'Dataset_3','29959','240406';...
        %'Dataset_3','29959','240407';...
        ...
        };
    
    
    RPSW_CrosCorr = [];
    SWSW_AutoCorr = [];
    
    %=== Load data and aggregate them in the Multi_Day structure
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   SB = []; ids = [];
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(Rpl_table) && any(cellfun(@(row) isequal(row, Rpl_table.unique_ID(1,:)), num2cell(sessions2include, 2)))
            
            %=== Accumulate the replay table
            SB = [SB; Rpl_table];
            
            %=== Load SWR
            session_id = strsplit(FileList(nc).folder,'\'); session_id = session_id{end-1};
            data_directory = ['C:\Users\Angelo\Desktop\Temporary_Work\NP_datasets','\',session_id];
            cd(data_directory);
            load('RPL_probe2.mat');
            
            %=== Extract good replays
            Rpl_table.good =  rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (abs(z)> min_C) && (t<min_p), Rpl_table, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');
            Rpl_stable = Rpl_table(Rpl_table.good,:);
            
            %=== Calculate cross-correlation Replay and SWR and accumulate
            [RS_cross_corr,RS_bin_centers] = cross_correlogram_AF_v0(Rpl_stable.tL+Rpl_stable.dt/2,RPL_out.table.t,0.5,0.025);
            [SS_cross_corr,RS_bin_centers] = cross_correlogram_AF_v0(RPL_out.table.t,RPL_out.table.t,0.5,0.025);
            
            %=== Accumulate cross-correlations
            RPSW_CrosCorr = [RPSW_CrosCorr,RS_cross_corr];
            SWSW_AutoCorr = [SWSW_AutoCorr,SS_cross_corr];
            
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except SB RPSW_CrosCorr SWSW_AutoCorr RS_bin_centers
end

%% LOOK AT RIPPLE-SWR CROSSCORRELATIONS
for hide=1
    
    figure('units','normalized','outerposition',[.3 .3 .1 .35]);
    hold on;
    bar(RS_bin_centers,mean(RPSW_CrosCorr,2,'omitnan'),'FaceColor','k','EdgeColor','none','FaceAlpha',0.5);
    bar(RS_bin_centers,-mean(SWSW_AutoCorr,2,'omitnan'),'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
    xlabel('Time(s)');  ylabel('Fraction');
    
    
end

%% LOOK AT FORWARD VS REVERSE REPLAYS
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 5;                      % Minimum number of active cells
    min_f = 0.2;                    % Minimum fraction of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.1;                    % Minimum Correlation
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                           % Id for each cluster from each session
    SB.GoodFReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z> min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.GoodRReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z<-min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.Not__Replay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z>-inf  ) && (t<=inf ), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');  % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                    % Derived duration
    
    %=== Select good events
    FRpl_table = SB(SB.GoodFReplay,:);   nFR = size(FRpl_table,1);  % Forward Replays
    RRpl_table = SB(SB.GoodRReplay,:);   nRR = size(RRpl_table,1);  % Reverse Replays
    NRpl_table = SB(SB.Not__Replay,:);   nNR = size(NRpl_table,1);  % Reverse Replays
    
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
    
    
end

%% EXTRACT EVENTS OF INTEREST (FORWARD) AND LOOK AT SPEED
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 5;                      % Minimum number of active cells
    min_p = 0.05;                   % Minimum p value
    min_C = 0.6;                    % Minimum Correlation
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                           % Id for each cluster from each session
    SB.GoodReplay = rowfun(@(x,y,z) (x>min_n) && (y>min_C) && (z<min_p), SB, 'InputVariables', {'n','C','p'},'OutputFormat', 'uniform');    % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                    % Derived duration
    
    %=== Select good events
    Rpl_table = SB(SB.GoodReplay,:);
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','ascend');
    
    %=== Summary statistics for replays from a given bat & cluster
    RPL_Clus = groupsummary(Rpl_table,'groupID','median',{'flight_len','dt','v','c_ratio','derived_dt','C','n'});
    
    %=== Plot cumulative distribution for flight lengths and example trajectories
    figure('units','normalized','outerposition',[.2 0.2 .15 0.3]);
    plot(sort(RPL_Clus.median_flight_len)); ylim([0 max(RPL_Clus.median_flight_len)]);  xlabel('Flight Cluster Id');    ylabel('Lenght (m)');
    
    %=== Plot a few trajectory types
    figure('units','normalized','outerposition',[0 0.4 1 0.3]);
    tiledlayout(1,15,'TileSpacing','compact');
    for jj=round(linspace(1,size(Rpl_table_ascending,1),15))
        nexttile;
        clus_path = Rpl_table_ascending.mean_path3D{jj,1};
        plot3(clus_path(:,1),clus_path(:,2),clus_path(:,3),'-','LineWidth',1,'Color', 'k');
        view(0,90);    axis equal;  xlim([-6 6]); ylim([-3 3]);   zlim([0 3]);  title([num2str(Rpl_table_ascending.flight_len(jj),3),' m']);
    end
    
    figure('units','normalized','outerposition',[.3 .3 .25 .3]);
    tiledlayout(1,2,'TileSpacing','compact');
    nexttile;   scatter(med_curv,gaussR2,15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    xlabel('Max curvature');    ylabel('R2 Gaussian');
    idx = dbscan([med_curv,gaussR2],0.3,5);
    nexttile;   gscatter(med_curv,gaussR2,idx);
    
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
    %=== Plot replays
    rng(2);
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','descend');
    cRT_good = find(Rpl_table_ascending.n>6 & Rpl_table_ascending.p<0.05 & Rpl_table_ascending.C>0.6 & Rpl_table_ascending.f>0.1);
    n_event2show = min(150,size(cRT_good,1));
    cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
    Rpl_subtable = Rpl_table_ascending(cRT_good,:);
    figure('units','normalized','outerposition',[.2 0 .6 1]);
    tiledlayout(10,15,'TileSpacing','compact');
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
    
    %=== Plot replays (vertical version)
    Rpl_table_ascending = sortrows(Rpl_table,'flight_len','descend');
    counter =1;
    wR = 0.07;
    for sd = round(linspace(1,14,14))
        rng(sd);
        cRT_good = find(Rpl_table_ascending.n>7 & Rpl_table_ascending.p<0.01 & Rpl_table_ascending.C>0.7 & Rpl_table_ascending.f>0.1);
        n_event2show = min(15,size(cRT_good,1));
        cRT_subset = sort(datasample(cRT_good,n_event2show,'Replace',false));
        Rpl_subtable = Rpl_table_ascending(cRT_good,:);
        figure('units','normalized','outerposition',[wR*(counter-1) 0 wR 1]);
        tiledlayout(15,1,'TileSpacing','compact');
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
            legend([num2str(Rpl_subtable.flight_len(cRT_subset(i)),2),' m'],'Location','eastoutside');
        end
        sgtitle(['Seed: ', num2str(sd)]);
        counter=counter+1;
    end
    
    
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
    scatter(abs(pw_distance(subset)),pw_timing(subset),15,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',.1); hold on;
    d_thresh = [5 20];
    colormap('cool');colorbar;
    title(['Min ',num2str(5),' active cells'])
    ylabel('dt during F replay (s)'); xlabel('Place field distance (m)');
    
    %=== Pairwise timings vs field distances
    Rpl_table_sst = Rpl_table(Rpl_table.n>7 & Rpl_table.p<0.01 & Rpl_table.C>0.6,:);
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
    hold on;plot(fit_result);
    
end

%% LOOK AT AVERAGE DECODED SHAPE OF REPLAYS AND FLIGHT TYPE ANALYSIS
for hide=1
    
    %=== Params
    n_events = size(SB,1);          % Number of candidate replays
    min_n = 6;                      % Minimum number of active cells (6)
    min_f = 0.1;                    % Minimum fraction of active cells (0.1)
    min_p = 0.05;                   % Minimum p value (0.05)
    min_C = 0.6;                    % Minimum Correlation (0.6)
    
    %=== Add a few features
    [~,~,SB.groupID] =  unique([string(SB.unique_ID),string(SB.clus_id)],'rows');                                                                               % Id for each cluster from each session
    SB.GoodReplay = rowfun(@(x,y,z,t) (x>min_n) && (y>min_f) && (z> min_C) && (t<min_p), SB, 'InputVariables', {'n','f','C','p'},'OutputFormat', 'uniform');    % Good Replay or not
    SB.derived_dt = SB.flight_len./SB.v;                                                                                                                        % Derived duration
    
    %=== Select good events
    Rpl_table = SB(SB.GoodReplay,:);
    n_good = size(Rpl_table,1);
    
    %=== Calculate weighted correlation
    for i=1:n_good
        
        posterior = Rpl_table.posterior{i,1}(3:end-2,:);                % Get the posterior probability (cut the edges)
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial');
        [max_P,max_L] = max(pstP);
        invalidP = (max_P==0 | isnan(max_P));                           % Invalidate NaNs or zeros
        max_P(invalidP)=[]; max_L(invalidP)=[];
        corrmat = weightedcorrs([[1:numel(max_L)]',max_L'], max_P');
        Rpl_table.weightcorr(i) = corrmat(1,2);
        
    end
    
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
        mean_v1D_rsz(:,jj) = interp1(mean_v1D./max(mean_v1D),linspace(1,numel(mean_v1D),size(mean_v1D_rsz,1)),'linear')';
        %mean_v1D_rsz(:,jj) = interp1(mean_v1D,linspace(1,numel(mean_v1D),100),'linear')';
        
    end
    Rpl_table.med_curv =  med_curv(tmp_idx1);   % Assign to Rpl table entries
    Rpl_table.gaussR2 =  gaussR2(tmp_idx1);     % Assign to Rpl table entries
    
    %=== Cluster flight types
    %for i=1:30
    rng(10); warning('off');
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
    %end
    
    %plot(mean_v1D_rsz(:,flight_type==1))   'chebychev', i = 10
    
    %=== Plot a few trajectory types
    figure('units','normalized','outerposition',[0 0.4 1 0.3]);
    tiledlayout(1,15,'TileSpacing','compact');
    for jj=round(linspace(1,size(Rpl_table,1),15))
        nexttile;
        clus_path = Rpl_table.mean_path3D{jj,1};
        plot(clus_path(:,1),clus_path(:,2),'-','LineWidth',1,'Color', 'k');
        view(0,90);    axis equal;  xlim([-6 6]); ylim([-3 3]);   title([num2str(Rpl_table.flight_len(jj),3),' m']);
    end
    
    %=== Fit each event with a logistic function and resize the spatial bins
    logisEqn = '1/(1+exp(-b*(x-c)))';
    nS_res = 60;                            % Number of spatial bins for resampling (60)
    nT_res = 40;                            % Number of temporal bins for windowing (40)
    min_wc = 0.6;                           % Min weighted correlation (0.6)
    min_R2_sigm = 0.6;                      % Min R2 for logistic fitting (0.6)
    min_R2_g1v = 0.6;                       % Min R2 for fitting v profile with 1 Gaussian (0.6)
    replay_tbin_size = 0.005;
    replay_sbin_size = 0.15;
    Rpl_table.central_bin = NaN(n_good,1);
    Rpl_table.slope = NaN(n_good,1);
    Rpl_table.lrs = NaN(n_good,1);
    Rpl_table.res_posterior = Rpl_table.posterior;
    Rpl_table.cut_posterior = repmat({NaN(nS_res,nT_res)},n_good,1);
    accum_posterior1 = [];  accum_f1 = [];  accum_s1 = [];  accum_rd1 = []; accum_fd1 = []; accum_v1 = [];
    accum_posterior2 = [];  accum_f2 = [];  accum_s2 = [];  accum_rd2 = []; accum_fd2 = []; accum_v2 = [];
    for i=1:n_good
        
        posterior = Rpl_table.posterior{i,1}(3:end-2,:);                % Get the posterior probability (cut the edges)
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial'); % Smooth the posterior
        [max_P,max_L] = max(pstP);                                      % Calculate location at max posterior
        invalidP = (max_P==0 | isnan(max_P));                           % Invalidate NaNs or zeros
        y1 = normalize(max_L,'range')';                                 % Normalize between 0 and 1 and assign x,y and weights
        x1 = [1:numel(y1)]';  w = max_P;                                 % Assign x and weights
        y1(1) = 0;   y1(end)=1;                                          % Anchor to start and stop (to help fitting)
        x1(invalidP)=[]; w(invalidP)=[]; y1(invalidP)=[];                % Remove invalid points
        [f,gof] = fit(x1,y1,logisEqn,'Start',[3 x1(end)/2],'Weights',w); % Fit with logistic function
        %plot(f,x,y);    legend('off');
        fitcoeff = coeffvalues(f);                                      % Get the fit coefficients
        Rpl_table.central_bin(i) = fitcoeff(2);                         % Get the center of the logistics f
        Rpl_table.slope(i) = fitcoeff(1);                               % Get the slope of the logistics f
        Rpl_table.lrs(i) = gof.rsquare;                                 % Get the R square
        %Rpl_table.res_posterior{i,1} = imresize(pstP,[nS_res nT]);     % Resize the filtered posterior
        Rpl_table.res_posterior{i,1} = imresize(posterior,[nS_res nT]); % Resize the unfiltered posterior
        flight = diff(Rpl_table.mean_path1D{i,1})./diff(Rpl_table.mean_time1D{i,1});                            % Get the velocity
        
        %=== Calculate Replay Speed
        t_replay = [1:nT]'*replay_tbin_size;
        %[max_P,max_L]  = max(imresize(pstP,[nS_res nT]));
        [max_P,max_L]  = max(pstP);
        %         invalidP = (max_P==0 | isnan(max_P));
        %         max_L(invalidP) = NaN;
        replay_vel = [0, diff(max_L*replay_sbin_size)./diff(t_replay')]';
        
        %=== Now keep an interval around the center of the replay
        ct = round(Rpl_table.central_bin(i));
        interval = ct+[-nT_res/2:(nT_res/2-1)];
        if numel(intersect(1:nT,interval))==numel(interval) && Rpl_table.lrs(i)>min_R2_sigm && Rpl_table.weightcorr(i)>min_wc
            Rpl_table.cut_posterior{i,1}  = Rpl_table.res_posterior{i,1}(:,interval);
            if Rpl_table.gaussR2(i)>min_R2_g1v
                accum_posterior1 = cat(3,accum_posterior1,Rpl_table.cut_posterior{i,1});
                accum_f1 = [accum_f1, normalize(interp1(flight,linspace(1,numel(flight),100)),'range')'];
                accum_s1 = [accum_s1,   Rpl_table.slope(i)];
                accum_rd1 = [accum_rd1,   Rpl_table.dt(i)];
                accum_fd1 = [accum_fd1,    Rpl_table.mean_time1D{i,1}(end)];
                accum_v1 = [accum_v1, replay_vel(interval)];
            else
                accum_posterior2 = cat(3,accum_posterior2,Rpl_table.cut_posterior{i,1});
                accum_f2 = [accum_f2, normalize(interp1(flight,linspace(1,numel(flight),100)),'range')'];
                accum_s2 = [accum_s2,   Rpl_table.slope(i)];
                accum_rd2 = [accum_rd2,   Rpl_table.dt(i)];
                accum_fd2 = [accum_fd2,    Rpl_table.mean_time1D{i,1}(end)];
                accum_v2 = [accum_v2, replay_vel(interval)];
                
            end
        end
    end
    
    %=== Average the 'cut and resized' posterior
    avg_dec_event1 = squeeze(mean(accum_posterior1,3,'omitnan'));
    avg_dec_event2 = squeeze(mean(accum_posterior2,3,'omitnan'));
    avg_flight1 = mean(accum_f1,2);
    avg_flight2 = mean(accum_f2,2);
    avg_rpl_v1 = mean(smoothdata(accum_v1,2,'movmedian',1),2,'omitnan');
    avg_rpl_v2 = mean(smoothdata(accum_v2,2,'movmedian',1),2,'omitnan');
    %avg_rpl_v1 = mean(accum_v1,2,'omitnan');
    %avg_rpl_v2 = mean(accum_v2,2,'omitnan');
    
    %=== Fit
    [max_P,max_L] = max(avg_dec_event1);                                % Calculate location at max posterior
    invalidP = (max_P==0 | isnan(max_P));                               % Invalidate NaNs or zeros
    y1 = normalize(max_L,'range')';                                     % Normalize between 0 and 1 and assign x,y and weights
    x1 = [1:numel(y1)]';  w = max_P;                                    % Assign x and weights
    y1(1) = 0;   y1(end)=1;                                             % Anchor to start and stop (to help fitting)
    x1(invalidP)=[]; w(invalidP)=[]; y1(invalidP)=[];                   % Remove invalid points
    f1 = fit(x1,y1,logisEqn,'Start',[3 x1(end)/2],'Weights',w);         % Fit with logistic function
    [max_P,max_L] = max(avg_dec_event2);                                % Calculate location at max posterior
    invalidP = (max_P==0 | isnan(max_P));                               % Invalidate NaNs or zeros
    y2 = normalize(max_L,'range')';                                     % Normalize between 0 and 1 and assign x,y and weights
    x2 = [1:numel(y2)]';  w = max_P;                                    % Assign x and weights
    y2(1) = 0;   y2(end)=1;                                             % Anchor to start and stop (to help fitting)
    x2(invalidP)=[]; w(invalidP)=[]; y2(invalidP)=[];                   % Remove invalid points
    f2 = fit(x2,y2,logisEqn,'Start',[3 x2(end)/2],'Weights',w);         % Fit with logistic function
    
    %=== Plot replays
    smooth_f = [1 1];
    figure('units','normalized','outerposition',[.2 0.2 .12 0.4]);
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile;   imagesc(avg_dec_event1,prctile(avg_dec_event1,[1 99],'all')');
    colormap('hot');  set(gca,'YDir','normal');   title('Straight'); ylabel('Raw'); xticks([]); yticks([]);
    nexttile;   imagesc(avg_dec_event2,prctile(avg_dec_event2,[1 99],'all')');
    colormap('hot');  set(gca,'YDir','normal');   title('Loops');                   xticks([]); yticks([]);
    nexttile;   imagesc(imgaussfilt(avg_dec_event1,smooth_f),prctile(avg_dec_event1,[1 99.5],'all')');
    colormap('hot'); set(gca,'YDir','normal');  ylabel('Filtered');                 xticks([]); yticks([]);
    nexttile;   imagesc(imgaussfilt(avg_dec_event2,smooth_f),prctile(avg_dec_event2,[1 99.5],'all')');
    colormap('hot');  set(gca,'YDir','normal');                                     xticks([]); yticks([]);
    
    %=== Plot velocity profiles
    figure('units','normalized','outerposition',[.4 0.5 .33 0.25]);
    tiledlayout(1,4,'TileSpacing','compact');
    data_tmp = accum_f1';
    nexttile;   plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'r');
    data_tmp = accum_f2';
    hold on;   plotWinterval_AF_v0(1:size(data_tmp,2),mean(data_tmp,'omitnan'),mean(data_tmp,'omitnan')-std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),mean(data_tmp,'omitnan')+std(data_tmp,[],'omitnan')./sqrt(size(data_tmp,1)),'b');
    title('Flight');    xlabel('Flight %'); ylabel('Speed (norm)');    xticks([]); yticks([]);
    nexttile;
    m1 = smoothdata(avg_rpl_v1,'movmean',8);
    sem1 = smoothdata(std(accum_v1,[],2,'omitnan'),'movmean',8)./sqrt(size(accum_v1,2));
    plotWinterval_AF_v0(1:size(m1,1),m1,m1-sem1,m1+sem1,'r');   hold on;
    m2 = smoothdata(avg_rpl_v2,'movmean',8);
    sem2 = smoothdata(std(accum_v2,[],2,'omitnan'),'movmean',8)./sqrt(size(accum_v2,2));
    plotWinterval_AF_v0(1:size(m2,1),m2,m2-sem2,m2+sem2,'b');
    title('Replay');    xlabel('Replay %'); ylabel('Speed (norm)');    xticks([]); yticks([]); xlim('tight');
    nexttile; h1Fit = plot(f1,x1,y1);    legend('off');
    hold on;  h2Fit = plot(f2,x2,y2);    legend('off');
    xlabel('Time'); ylabel('Space');    title('Replay');    xticks([]); yticks([]);
    set(h1Fit(1), 'Marker', 'o', 'MarkerSize', 3, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
    set(h1Fit(2), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
    set(h2Fit(1), 'Marker', 'o', 'MarkerSize', 3, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    set(h2Fit(2), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2);
    nexttile;
    invalid_v1 = m1<0;    invalid_v2 = m2<0;
    t_bin1 = [1:size(m1,1)]'*replay_tbin_size;
    t_bin2 = [1:size(m2,1)]'*replay_tbin_size;
    m1(invalid_v1) = [];    t_bin1(invalid_v1) = [];
    m2(invalid_v2) = [];    t_bin2(invalid_v2) = [];
    m1 = m1./max(m1);
    m2 = m2./max(m2);
    [fg1,gof1] = fit(t_bin1,m1,'gauss1','StartPoint',[1 0.1 0.1]);
    [fg2,gof2] = fit(t_bin2,m2,'gauss2','StartPoint',[0.5 0.03 0.01 1 0.15 0.01]);
    t_bin = [0 : replay_tbin_size :  max([t_bin1;t_bin2])];
    y_1 = feval(fg1,t_bin);
    y_2 = feval(fg2,t_bin);
    hold on;
    plot(t_bin1,m1,'o','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','none');       plot(t_bin2,m2,'o','MarkerSize',3,'MarkerFaceColor','b','MarkerEdgeColor','none');
    plot(t_bin,y_1,'r','LineWidth',2);       plot(t_bin,y_2,'b','LineWidth',2);
    
    %=== Look at slope and duration
    figure('units','normalized','outerposition',[.5 0.2 .3 0.3]);
    tiledlayout(1,4,'TileSpacing','compact');
    nexttile;   plot_Data_AF_v0({accum_s1;accum_s2});
    xticks([1 2]);  xticklabels({'Straight','Loops'});  ylim([0 prctile([accum_s1,accum_s2],90,'all')]);    ylabel('Slope Coefficient');
    title(['p = ',num2str(ranksum(accum_s1,accum_s2),2)]);
    nexttile;   plot_Data_AF_v0({accum_rd1;accum_rd2});
    xticks([1 2]);  xticklabels({'Straight','Loops'}); ylabel('Replay Duration (s)');
    title(['p = ',num2str(ranksum(accum_rd1,accum_rd2),2)]);
    nexttile;   plot_Data_AF_v0({accum_fd1;accum_fd2});
    xticks([1 2]);  xticklabels({'Straight','Loops'}); ylabel('Flight Duration (s)');
    title(['p = ',num2str(ranksum(accum_fd1,accum_fd2),2)]);
    nexttile;   plot_Data_AF_v0({accum_fd1./accum_rd1;accum_fd2./accum_rd2});
    xticks([1 2]);  xticklabels({'Straight','Loops'}); ylabel('Compression');
    title(['p = ',num2str(ranksum(accum_fd1./accum_rd1,accum_fd2./accum_rd2),2)]);
    
    %=== Plot some example replays for straight and loop flights
    n_evt = 50;
    RPL_sst_s = Rpl_table(Rpl_table.gaussR2>0.5 & Rpl_table.n>7 & Rpl_table.weightcorr>0.5,:);
    RPL_sst_l = Rpl_table(Rpl_table.gaussR2<0.5 & Rpl_table.n>7 & Rpl_table.weightcorr>0.5,:);
    ssTs = datasample(1:size(RPL_sst_s,1),n_evt,'Replace',false);
    ssTl = datasample(1:size(RPL_sst_l,1),n_evt,'Replace',false);
    figure('units','normalized','outerposition',[0 0 1 0.3]);
    tiledlayout(2,n_evt,'TileSpacing','compact');
    for i=ssTs
        posterior = RPL_sst_s.res_posterior{i,1};                       % Get the posterior probability (cut the edges)
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial'); % Filter
        nexttile;   imagesc(pstP,prctile(pstP,[1 99],'all')')
        colormap('hot');  set(gca,'YDir','normal');                                     xticks([]); yticks([]);
    end
    for i=ssTl
        posterior = RPL_sst_l.res_posterior{i,1};                       % Get the posterior probability (cut the edges)
        nS = size(posterior,1); nT = size(posterior,2);                 % Number of space and time bins
        pstP = imgaussfilt(posterior,[.01 1],'FilterDomain','spatial'); % Filter
        nexttile;   imagesc(pstP,prctile(pstP,[1 99],'all')')
        colormap('hot');  set(gca,'YDir','normal');                                     xticks([]); yticks([]);
    end
    
end
