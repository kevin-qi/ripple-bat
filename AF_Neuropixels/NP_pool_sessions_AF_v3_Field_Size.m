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
        %'Dataset_2','14445','231210';...
        %'Dataset_2','14445','231211';...
        %'Dataset_2','14445','231212';...
        'Dataset_2','14611','231208';...
        'Dataset_2','14611','231209';...
        'Dataset_2','14611','231210';...
        'Dataset_2','14611','231211';...
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
    Folder = cd;    FileList = dir(fullfile(Folder, '**', 'Analyzed_NPs_*'));   SB = []; 
    for nc = 1:length(FileList)
        cd(FileList(nc).folder);
        load(FileList(nc).name);
        if ~isempty(NP_cluster) && any(cellfun(@(row) isequal(row, NP_cluster(1).unique_ID), num2cell(sessions2include, 2)))
            
            %=== Accumulate the replay table and single units
            SB = [SB; NP_cluster];
            
        end
        disp([num2str(length(FileList)-nc),' remaining sessions to load...']);
    end
    cd(Folder); clearvars -except SB 
    
end

%% CALCULATE A FEW FEATURES FOR EACH CLUSTER
for hide=1
    
    n_clusters = size(SB,1);
    
    for i=1:n_clusters
    
        %=== Extract the table of cells for that cluster
        NP_table = SB(i).table;
        
        %=== Extract place cells
        NP_table.place_cond = NP_table.spkPerflight>1 &...  % Min Spikes per flight (DEF: 1)
            NP_table.peakHz>3 &...        % Min Peak Firing Rate (DEF: 3)
            NP_table.stab_m>0.4 &...      % Min Stability (DEF: .4)
            NP_table.sff<0.5;             % Min Peakyness (DEF: .7)
        NP_subtable = NP_table(NP_table.place_cond,:);
        
        %=== Extract a few features
        SB(i).n_place_cells = size(NP_subtable,1);                                      % Number of place cells
        SB(i).meters_per_cell_1 = range(NP_subtable.field_loc_m)/SB(i).n_place_cells;     % Average space covered by one place cell (useful later)
        coeff_line = polyfit([1:SB(i).n_place_cells],sort(NP_subtable.field_loc_m),1);
        SB(i).meters_per_cell_2 = coeff_line(1);
        SB(i).avg_field_size = mean(NP_subtable.f_width);                               % Average field size
        SB(i).mdn_field_size = median(NP_subtable.f_width);                             % Median field size
        
    end
    
    %=== Extract relevant vectors
    f_l = [SB.length]';             % Flight_lenght
    n_c = [SB.n_place_cells]';      % Number of place cells
    f_s = [SB.avg_field_size]';     % Average field size
    p_d = [SB.meters_per_cell_2]';  % Average place distance
    
    %=== Show field size vs flight lenght
    mkr_size = 30;
    figure('units','normalized','outerposition',[.3 .1 .5 .5]);
    tiledlayout(2,5,'TileSpacing','tight');
    nexttile;   scatter([SB.length],[SB.n_place_cells],mkr_size,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);    xlabel('Flight Lenght (m)');    ylabel('N Place Cells');           xlim('tight'); lsline;   
    [corr_val,p_val] =  corr([SB.length]',[SB.n_place_cells]');    title([corr_val,p_val]); xlim_vals = xlim;   xlim([0 xlim_vals(2)]); ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   scatter([SB.length],[SB.meters_per_cell_1],mkr_size,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);xlabel('Flight Lenght (m)');    ylabel('m/cell (1)');              xlim('tight'); lsline;
    [corr_val,p_val] =  corr([SB.length]',[SB.meters_per_cell_1]');    title([corr_val,p_val]);  xlim_vals = xlim;   xlim([0 xlim_vals(2)]);    ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   scatter([SB.length],[SB.meters_per_cell_2],mkr_size,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);xlabel('Flight Lenght (m)');    ylabel('Place Cell Distance (m)'); xlim('tight'); lsline;
    [corr_val,p_val] =  corr([SB.length]',[SB.meters_per_cell_2]');    title([corr_val,p_val]);  xlim_vals = xlim;   xlim([0 xlim_vals(2)]);    ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   scatter([SB.length],[SB.avg_field_size],mkr_size,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   xlabel('Flight Lenght (m)');    ylabel('Mean Field Size (m)');     xlim('tight'); lsline;
    [corr_val,p_val] =  corr([SB.length]',[SB.avg_field_size]');    title([corr_val,p_val]);     xlim_vals = xlim;   xlim([0 xlim_vals(2)]);    ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   scatter([SB.length],[SB.mdn_field_size],mkr_size,'MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',.5);   xlabel('Flight Lenght (m)');    ylabel('Median Field Size (m)');   xlim('tight'); lsline;
    [corr_val,p_val] =  corr([SB.length]',[SB.mdn_field_size]');    title([corr_val,p_val]);     xlim_vals = xlim;   xlim([0 xlim_vals(2)]);    ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    
    nexttile;   plot_distr_AF_v0(n_c(f_l>3&f_l<7), n_c((f_l>7&f_l<13)), {'Short Flights', 'Long Flights'}, 'SEM', 'Number of place cells'); ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;
    nexttile;   plot_distr_AF_v0(f_s(f_l>3&f_l<7), f_s((f_l>7&f_l<13)), {'Short Flights', 'Long Flights'}, 'SEM', 'Place-Field Size (m)');  ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;   plot_distr_AF_v0(p_d(f_l>3&f_l<7), p_d((f_l>7&f_l<13)), {'Short Flights', 'Long Flights'}, 'SEM', 'Place-Cell Distance (m)');   ylim_vals = ylim;   ylim([0 ylim_vals(2)]);
    nexttile;

    %=== Look at the level of single bat, single sessions
    SB_table = struct2table(SB);    size_corr = []; dist_corr = []; numb_corr = [];
    [~,~,SB_table.sessionID] =  unique(string(SB_table.unique_ID),'rows'); 
    for i=1:max(SB_table.sessionID)                 % Calculate correlation between field size (or distance) and cluster lenght
        SB_sst = SB_table(SB_table.sessionID==i,:);
        if size(SB_sst,1)>1
            dist_corr = [dist_corr; corr(SB_sst.meters_per_cell_2,SB_sst.length)];
            size_corr = [size_corr; corr(SB_sst.avg_field_size,SB_sst.length)];
            numb_corr = [numb_corr; corr(SB_sst.n_place_cells,SB_sst.length)];
        end
    end
    
    figure('units','normalized','outerposition',[.3 .1 .3 .3]);
    tiledlayout(1,3,'TileSpacing','tight');
    nexttile;   histogram(numb_corr,[-1:0.2:1],'Normalization','probability');  ylabel('Fraction'); xlabel('Correlation Place Number to Flight Lenght');
    nexttile;   histogram(dist_corr,[-1:0.2:1],'Normalization','probability');  ylabel('Fraction'); xlabel('Correlation Place Distance to Flight Lenght');
    nexttile;   histogram(size_corr,[-1:0.2:1],'Normalization','probability');  ylabel('Fraction'); xlabel('Correlation Place Size to Flight Lenght');
    

end

